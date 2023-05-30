# Import Libraries
import pandas as pd
import numpy as np
import yaml
import math
import os
import sys
import geopandas as gpd


def calculate_route_vmt(transit_routes_df: pd.DataFrame, transit_flows_df: pd.DataFrame, periods_length: dict, daily_factor: float):
    """
    calculates daily vmt of each route using period headway and route length

    args:
        transit_routes_df: dataframe with headway (by AM, PM, OP period) for routes
        transit_flows_df: dataframe with TOMP (to mile post) used to calculate the length of route
        periods_length: dictionary of period lengths in hours
        daily_factor: factor to account for running hours of transit in a day

    returns:
        dataframe with route id and its daily vmt

    """

    work_df = transit_routes_df.copy()
    work_df["Total_Freq"] = 0

    transit_route_length_df = transit_flows_df.groupby(["ROUTE"], as_index=False)["TOMP"].max()
    transit_route_length_df.rename(columns={"TOMP": "Route_Length"}, inplace=True)

    for period in periods_length.keys():
        work_df.loc[:, period + "_Freq"] = np.where(
            work_df[period + "_Headway"] > 0,
            (60/work_df[period + "_Headway"]) * periods_length[period] * daily_factor,
            0
        )
        work_df["Total_Freq"] = work_df["Total_Freq"] + work_df[period + "_Freq"]

    work_df["NT_Freq"] = np.where(
        work_df["Night_Headway"] > 0,
        (60/work_df["Night_Headway"]) * work_df["Night_Hours"] * daily_factor,
        0
    )

    work_df["Total_Freq"] = work_df["Total_Freq"] + work_df["NT_Freq"]

    work_df = pd.merge(
        work_df,
        transit_route_length_df,
        left_on="Route_ID",
        right_on="ROUTE",
        how="left"
    )

    work_df["VMT"] = work_df["Total_Freq"] * work_df["Route_Length"]

    out_df = work_df[["Route_ID", "Route_Name", "Mode", "VMT"]]

    return out_df


def calculate_vmt_reduction(transit_ev_strategy_df: pd.DataFrame, route_vmt_df: pd.DataFrame):
    """
    calculates reduction in vmt by route based on change in conventional fuel percent b/w base year and scenario year

    args:
        transit_ev_strategy_df: dataframe with conventional fuel type percent by route in base year and scenario year
        route_vmt_df: dataframe with daily vmt by route

    returns:
        dataframe with route id and vmt reductions (by fuel type)

    """

    work_df = pd.merge(
        transit_ev_strategy_df,
        route_vmt_df,
        left_on="route",
        right_on="Route_ID",
        how="left"
    )

    # keep only the routes by bus mode
    # in case there is a non-bus route specified in the strategy file
    work_df = work_df[work_df["Mode"].isin([6, 7, 8, 9, 10])]

    work_df["total_fuel_vmt_reduction"] = ((work_df["conventional_fuel_pct_base"]
                                            - work_df["conventional_fuel_pct_scen"])
                                           * work_df["VMT"])

    work_df["gasoline_vmt_reduction"] = (work_df["total_fuel_vmt_reduction"]
                                         * work_df["gasoline_fleet_proportion"])

    work_df["cng_vmt_reduction"] = (work_df["total_fuel_vmt_reduction"]
                                    * work_df["cng_fleet_proportion"])

    work_df.rename(columns={"Route_Name": "route_name"}, inplace=True)
    work_df.rename(columns={"route": "route_id"}, inplace=True)

    work_df["route"] = (work_df["route_name"]/1000).astype(int)

    out_df = work_df[["route_id", "route", "route_name",
                      "total_fuel_vmt_reduction", "gasoline_vmt_reduction", "cng_vmt_reduction"]]

    return out_df


def get_transit_emission_factors(input_emission_data: pd.DataFrame, scen_year: int):
    """
    prepare a dataframe of transit emission factors for the scen year

    args:
        input_emission_data: dataframe with all emission rates
        scen_year: scen_year

    returns:
        dataframe with co2 runex emission rates for transit buses

    """

    emission_df = input_emission_data[input_emission_data["Year"] == scen_year]
    emission_df = emission_df[emission_df["Vehicle Type"].isin(["Bus - Gas", "Bus - NG"])]
    emission_df.reset_index(inplace=True, drop=True)

    emission_df["Type"] = np.where(emission_df["Vehicle Type"] == "Bus - Gas", "gasoline", "cng")

    out_df = emission_df[["Year", "Type", "CO2 RunEx Emission Factor (gr/mile)"]]

    return out_df


def calculate_ghg_reduction(vmt_reduction: pd.DataFrame, emission_factors: pd.DataFrame):
    """
    calculates ghg reductions (by fuel type)

    args:
        vmt_reduction: dataframe with vmt reduction by route and by fuel type
        emission_factors: emission factor (ton/mile) by fuel type

    returns:
        dataframe with ghg reduction by fuel type

    """
    GRAMS_TO_SHORT_TONS = 0.0000011

    work_df = vmt_reduction.copy()

    work_df = pd.melt(
        vmt_reduction_df,
        id_vars=["route_id"],
        value_vars=["gasoline_vmt_reduction", "cng_vmt_reduction"],
        var_name="fuel_type",
        value_name="vmt_reduction"
    )

    work_df["fuel_type"] = work_df["fuel_type"].str.replace("_vmt_reduction", "")

    work_df = pd.merge(work_df, emission_factors, left_on="fuel_type", right_on="Type", how="left")

    work_df["ghg_reduction"] = (work_df["vmt_reduction"]
                                * work_df["CO2 RunEx Emission Factor (gr/mile)"]
                                * GRAMS_TO_SHORT_TONS)

    out_df = work_df.groupby(["fuel_type"])["ghg_reduction"].sum().to_frame()
    out_df.reset_index(inplace=True)
    out_df = out_df.append(
        {'fuel_type': 'TOTAL', 'ghg_reduction': out_df["ghg_reduction"].sum()}, ignore_index=True)
    out_df["fuel_type"] = out_df["fuel_type"].str.upper()
    out_df.rename(columns={"ghg_reduction": "ghg_reduction (short tons)"}, inplace=True)

    return out_df


def write_results(results_dict, out_file_name, out_dir):
    with pd.ExcelWriter(os.path.join(out_dir, out_file_name)) as writer:
        for key, value in results_dict.items():
            value.to_excel(writer, sheet_name=key, index=False)


if __name__ == "__main__":
    args = sys.argv
    config_filename = args[1]

    if not os.path.exists(config_filename):
        msg = "Configuration file doesn't exist at: {}".format(config_filename)
        raise ValueError(msg)

    with open(config_filename, "r") as yml_file:
        config = yaml.safe_load(yml_file)

    # inputs
    transit_routes_file = config['inputs']['transit_routes_file']
    transit_flow_file = config['inputs']['transit_flow_file']
    emission_factors_file = config['inputs']['emission_factors_file']
    transit_ev_strategy_file = config['inputs']['transit_ev_strategy_file']

    # parameters
    scen_year = config['parameters']['scen_year']
    transit_working_hours_factor = config['parameters']['transit_working_hours_factor']
    periods_length = {'AM': 3, "PM": 3.5, "OP": 6.5}

    # read config: outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # read data
    transit_routes_df = pd.read_csv(transit_routes_file)
    transit_flow_df = pd.read_csv(transit_flow_file)
    emission_df = pd.read_excel(emission_factors_file)
    transit_ev_strategy_df = pd.read_csv(transit_ev_strategy_file)

    # raise error if there are transit routes in strategy file that are not part of scenario
    scenario_routes = transit_routes_df["Route_ID"].to_list()
    strategy_routes = transit_ev_strategy_df["route"].to_list()
    missing_route_list = list(set(strategy_routes) - set(scenario_routes))

    if missing_route_list:
        msg = "ERROR: These transit routes {} are in the EV strategy input file but are not part of the scenario.".format(
            missing_route_list)
        raise ValueError(msg)

    # calculate total daily vmt for each route
    route_vmt_df = calculate_route_vmt(
        transit_routes_df,
        transit_flow_df,
        periods_length,
        transit_working_hours_factor
    )

    # calculate vmt reduction by strategy routes
    vmt_reduction_df = calculate_vmt_reduction(transit_ev_strategy_df, route_vmt_df)

    # get emission factors
    emission_factors_df = get_transit_emission_factors(emission_df, scen_year)

    # calculate ghg reductions
    ghg_reduction_df = calculate_ghg_reduction(vmt_reduction_df, emission_factors_df)

    # writing the result
    results_dict = {"GHG_Reduction": ghg_reduction_df,
                    "VMT_Reduction": vmt_reduction_df,
                    "Emission_Factors": emission_factors_df}

    write_results(results_dict, output_results_filename, output_dir)
