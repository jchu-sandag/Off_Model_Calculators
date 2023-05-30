# Import Libraries
import pandas as pd
import numpy as np
import yaml
import math
import os
import sys
import geopandas as gpd


def calculate_truck_vmt(highway_load_data_dict: dict, corridor_links_df: pd.DataFrame):
    """
    calculates daily truck vmt of each route for each truck class

    args:
        highway_load_data_dict: dictionary for highway load dataframe for different time periods
        corridor_links_df: dataframe with highway links for the corridor

    returns:
        dataframe with vmt by truck class

    """
    first = True

    for period in highway_load_data_dict.keys():
        highway_load_df = highway_load_data_dict[period]

        vmt_period_df = pd.merge(
            corridor_links_df,
            highway_load_df,
            how="left",
            left_on="hwycov_id",
            right_on="ID1"
        )
        vmt_period_df = vmt_period_df.fillna(0)

        light_truck_vol_cols = [x for x in vmt_period_df.columns if "Flow_lhd" in x]
        medium_truck_vol_cols = [x for x in vmt_period_df.columns if "Flow_mhd" in x]
        heavy_truck_vol_cols = [x for x in vmt_period_df.columns if "Flow_hhd" in x]

        vmt_period_df["light_flow"] = vmt_period_df.loc[:, light_truck_vol_cols].sum(axis=1)
        vmt_period_df["medium_flow"] = vmt_period_df.loc[:, medium_truck_vol_cols].sum(axis=1)
        vmt_period_df["heavy_flow"] = vmt_period_df.loc[:, heavy_truck_vol_cols].sum(axis=1)

        vmt_period_df["light"] = vmt_period_df["light_flow"] * vmt_period_df["len_mile"]
        vmt_period_df["medium"] = vmt_period_df["medium_flow"] * vmt_period_df["len_mile"]
        vmt_period_df["heavy"] = vmt_period_df["heavy_flow"] * vmt_period_df["len_mile"]

        vmt_period_df = vmt_period_df[["corridor_name", "light", "medium", "heavy"]]

        if first:
            vmt_all_period_df = vmt_period_df
            first = False
        else:
            vmt_all_period_df = pd.concat([vmt_all_period_df, vmt_period_df], axis=0)

    vmt_df = vmt_all_period_df.groupby(['corridor_name']).sum()
    vmt_df.reset_index(inplace=True)

    out_df = pd.melt(
        vmt_df,
        id_vars=["corridor_name"],
        value_vars=["light", "medium", "heavy"],
        var_name="truck_type",
        value_name="vmt"
    )

    return out_df


def get_truck_emission_factors(input_emission_data: pd.DataFrame, scen_year: int):
    """
    prepare a dataframe of truck emission factors for the scen year

    args:
        input_emission_data: dataframe with all emission rates
        scen_year: scen_year

    returns:
        dataframe with co2 runex emission rates for trucks

    """

    emission_df = input_emission_data[input_emission_data["Year"] == scen_year]
    emission_df = emission_df[emission_df["Vehicle Type"].isin(
        ["Light HDT", "Medium HDT", "Heavy HDT"])]
    emission_df.reset_index(inplace=True, drop=True)
    emission_df["Vehicle Type"] = emission_df["Vehicle Type"].str.replace(" HDT", "")
    emission_df["Vehicle Type"] = emission_df["Vehicle Type"].str.lower()

    out_df = emission_df[["Year", "Vehicle Type", "CO2 RunEx Emission Factor (gr/mile)"]]

    return out_df


def calculate_ghg_reduction(truck_vmt_df: pd.DataFrame, emission_factors: pd.DataFrame, freight_ev_strategy_df: pd.DataFrame):
    """
    calculates ghg reductions by truck type

    args:
        truck_vmt_df: dataframe with daily vmt by truck type
        emission_factors: truck emission factor (gr/mile) by truck type
        freight_ev_strategy_df: dataframe with ev conversion factor by truck type for the corridor

    returns:
        dataframe with ghg reduction by truck type

    """
    GRAMS_TO_SHORT_TONS = 0.0000011

    work_df = pd.merge(
        truck_vmt_df,
        emission_factors,
        left_on="truck_type",
        right_on="Vehicle Type",
        how="left"
    )

    work_df = pd.merge(
        work_df,
        freight_ev_strategy_df,
        on="truck_type",
        how="left"
    )

    work_df["ghg_reduction"] = (work_df["vmt"]
                                * work_df["CO2 RunEx Emission Factor (gr/mile)"]
                                * work_df["ev_conversion_factor"]
                                * GRAMS_TO_SHORT_TONS)

    corridor_name = work_df.loc[1]["corridor_name"]
    year = work_df.loc[1]["Year"]

    out_df = work_df[["Year", "corridor_name", "truck_type", "ghg_reduction"]]

    out_df = out_df.append({
        "Year": year,
        "corridor_name": corridor_name,
        "truck_type": 'TOTAL',
        "ghg_reduction": out_df["ghg_reduction"].sum()},
        ignore_index=True)

    out_df["truck_type"] = out_df["truck_type"].str.upper()
    out_df.rename(
        columns={"ghg_reduction": "ghg_reduction (short tons)", "Year": "year"},
        inplace=True
    )

    return out_df


def reformat_vmt_data(truck_vmt_df: pd.DataFrame, scen_year: int):
    """
    reformat the vmt by truck type dataframe

    args:
        truck_vmt_df: dataframe with daily vmt by truck type
        scen_year: scenario year
    returns:
        dataframe with vmt by truck type
    """
    out_df = truck_vmt_df.copy()

    corridor_name = out_df.loc[1]["corridor_name"]

    out_df["year"] = scen_year

    out_df = out_df.append({
        "year": scen_year,
        "corridor_name": corridor_name,
        "truck_type": 'TOTAL',
        "vmt": out_df["vmt"].sum()},
        ignore_index=True)

    out_df["truck_type"] = out_df["truck_type"].str.upper()
    out_df = out_df[["year", "corridor_name", "truck_type", "vmt"]]

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
    highway_load_am_file = config['inputs']['highway_load_am_file']
    highway_load_md_file = config['inputs']['highway_load_md_file']
    highway_load_pm_file = config['inputs']['highway_load_pm_file']
    highway_load_ea_file = config['inputs']['highway_load_ea_file']
    highway_load_ev_file = config['inputs']['highway_load_ev_file']
    highway_network_shapefile = config['inputs']['highway_network_shapefile']
    cmcp_corridors_shapefile = config['inputs']['cmcp_corridors_shapefile']
    emission_factors_file = config['inputs']['emission_factors_file']
    freight_ev_strategy_file = config['inputs']['freight_ev_strategy_file']

    # outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # parameters
    scen_year = config['parameters']['scen_year']
    corridor_name = config['parameters']['corridor_name']

    # read data
    highway_network_df = gpd.read_file(highway_network_shapefile)
    cmcp_corridors_df = gpd.read_file(cmcp_corridors_shapefile)

    highway_load_am_df = pd.read_csv(highway_load_am_file)
    highway_load_pm_df = pd.read_csv(highway_load_pm_file)
    highway_load_md_df = pd.read_csv(highway_load_md_file)
    highway_load_ea_df = pd.read_csv(highway_load_ea_file)
    highway_load_ev_df = pd.read_csv(highway_load_ev_file)

    emission_df = pd.read_excel(emission_factors_file)
    freight_ev_strategy_df = pd.read_csv(freight_ev_strategy_file)

    highway_load_data_dict = {"am": highway_load_am_df, "pm": highway_load_pm_df,
                              "md": highway_load_md_df, "ea": highway_load_ea_df, "ev": highway_load_ev_df}

    # get the corridor name
    corridor_list = cmcp_corridors_df["Name"].drop_duplicates().to_list()
    corridor_list.sort()
    if corridor_name not in corridor_list:
        msg = "Specified corridor_name: '{}' in config file is not correct. Corridor name must be one of these values: [{}].".format(
            corridor_name, "; ".join(corridor_list))
        raise ValueError(msg)

    # join corridor and highway shapefiles
    corridor_links_df = gpd.sjoin(
        highway_network_df,
        cmcp_corridors_df,
        how='inner',
        predicate='intersects'
    )
    corridor_links_df = corridor_links_df[corridor_links_df["Name"] == corridor_name]
    corridor_links_df = corridor_links_df[["Name", "hwycov_id", "len_mile"]]
    corridor_links_df.rename(columns={"Name": "corridor_name"}, inplace=True)

    # calculate vmt by truck class
    truck_vmt_df = calculate_truck_vmt(highway_load_data_dict, corridor_links_df)

    # get emission factors
    emission_factors_df = get_truck_emission_factors(emission_df, scen_year)

    # calculate ghg reduction by truck class and fuel type
    ghg_reduction_df = calculate_ghg_reduction(
        truck_vmt_df,
        emission_factors_df,
        freight_ev_strategy_df
    )

    # reformat data before writing it out
    truck_vmt_df = reformat_vmt_data(truck_vmt_df, scen_year)

    # writing results
    results_dict = {"GHG_Reduction": ghg_reduction_df,
                    "VMT_Truck": truck_vmt_df,
                    "Emission_Factors": emission_factors_df}

    write_results(results_dict, output_results_filename, output_dir)
