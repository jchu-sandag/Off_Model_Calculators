import pandas as pd
import numpy as np
import openmatrix as omx
import yaml
import os
import sys


def get_msa_travel_times(base_skim_file, scen_skim_file, military_taz_list, geo_xwalk, msa_names, external_gateways, sov_core_name, hov_core_name):
    # read travel time skims and calculate time savings by msa-msa pair
    travel_time_base_df = get_omx_cores_as_df(
        base_skim_file,
        [sov_core_name],
        ["sov_time_base"]
    )

    travel_time_scen_df = get_omx_cores_as_df(
        scen_skim_file,
        [sov_core_name, hov_core_name],
        ["sov_time_scen", "hov_time_scen"]
    )

    travel_time_df = pd.merge(
        travel_time_base_df,
        travel_time_scen_df,
        on=["orig_taz", "dest_taz"],
        how="left"
    )

    gateways_df = external_gateways.copy()
    gateways_df = gateways_df[["Home_County", "TAZ"]]
    gateways_df.loc[:, "Home_County"] = gateways_df["Home_County"].str.upper()

    taz_msa_df = geo_xwalk[["taz", "msa"]].drop_duplicates(subset=["taz"])

    travel_time_df = pd.merge(
        travel_time_df,
        taz_msa_df,
        left_on="orig_taz",
        right_on="taz",
        how="left"
    )
    travel_time_df = pd.merge(travel_time_df, msa_names_df, on="msa", how="left")
    travel_time_df.rename(columns={"msa_name": "orig_msa"}, inplace=True)
    travel_time_df.drop(["taz", "msa"], axis=1, inplace=True)

    travel_time_df = pd.merge(
        travel_time_df,
        taz_msa_df,
        left_on="dest_taz",
        right_on="taz",
        how="left"
    )
    travel_time_df = pd.merge(travel_time_df, msa_names_df, on="msa", how="left")
    travel_time_df.rename(columns={"msa_name": "dest_msa"}, inplace=True)
    travel_time_df.drop(["taz", "msa"], axis=1, inplace=True)

    travel_time_df = pd.merge(
        travel_time_df,
        gateways_df,
        left_on="orig_taz",
        right_on="TAZ",
        how="left"
    )
    travel_time_df["orig_msa"] = np.where(travel_time_df["Home_County"].isna(),
                                          travel_time_df["orig_msa"],
                                          travel_time_df["Home_County"])
    travel_time_df.rename(columns={"TAZ": "orig_ext_taz"}, inplace=True)
    travel_time_df.drop(["Home_County"], axis=1, inplace=True)

    travel_time_df = pd.merge(
        travel_time_df,
        gateways_df,
        left_on="dest_taz",
        right_on="TAZ",
        how="left"
    )
    travel_time_df["dest_msa"] = np.where(travel_time_df["Home_County"].isna(),
                                          travel_time_df["dest_msa"],
                                          travel_time_df["Home_County"])
    travel_time_df.rename(columns={"TAZ": "dest_ext_taz"}, inplace=True)
    travel_time_df.drop(["Home_County"], axis=1, inplace=True)

    travel_time_df = travel_time_df[travel_time_df.orig_msa.notna()]
    travel_time_df = travel_time_df[travel_time_df.dest_msa.notna()]
    travel_time_df = travel_time_df[
        ~(travel_time_df.orig_ext_taz.notna() & travel_time_df.dest_ext_taz.notna())
    ]

    travel_time_msa_df = travel_time_df.copy()
    travel_time_msa_df["time_savings_scen"] = (travel_time_msa_df["hov_time_scen"]
                                               - travel_time_msa_df["sov_time_scen"])
    travel_time_msa_df = travel_time_msa_df.groupby(["orig_msa", "dest_msa"]).agg(
        {"sov_time_base": 'mean', "time_savings_scen": 'mean'}
    )
    travel_time_msa_df.reset_index(inplace=True)

    travel_time_msa_mil_df = travel_time_df.copy()
    travel_time_msa_mil_df = travel_time_msa_mil_df[
        travel_time_msa_mil_df["dest_taz"].isin(military_taz_list)
    ]
    travel_time_msa_mil_df["time_savings_scen"] = (travel_time_msa_mil_df["hov_time_scen"]
                                                   - travel_time_msa_mil_df["sov_time_scen"])
    travel_time_msa_mil_df = travel_time_msa_mil_df.groupby(["orig_msa", "dest_msa"]).agg(
        {"sov_time_base": 'mean', "time_savings_scen": 'mean'}
    )
    travel_time_msa_mil_df.reset_index(inplace=True)
    travel_time_msa_mil_df.rename(
        columns={'sov_time_base': 'sov_time_base_mil', 'time_savings_scen': 'time_savings_scen_mil'}, inplace=True
    )

    out_df = pd.merge(
        travel_time_msa_df,
        travel_time_msa_mil_df,
        on=["orig_msa", "dest_msa"],
        how="left"
    )
    out_df["sov_time_base"] = out_df["sov_time_base"].fillna(0)
    out_df["time_savings_scen"] = out_df["time_savings_scen"].fillna(0)
    out_df["sov_time_base_mil"] = out_df["sov_time_base_mil"].fillna(0)
    out_df["time_savings_scen_mil"] = out_df["time_savings_scen_mil"].fillna(0)

    return out_df


def get_msa_skims(omx_file_name, core_names, out_cols, geo_xwalk, msa_names):
    time_df = get_omx_cores_as_df(omx_file_name, core_names, out_cols)

    out_df = pd.merge(time_df, geo_xwalk, left_on="orig_taz", right_on="taz", how="left")
    out_df = pd.merge(out_df, msa_names, on="msa", how="left")
    out_df.rename(columns={"msa_name": "orig_msa"}, inplace=True)
    out_df.drop(["mgra", "taz", "msa"], axis=1, inplace=True)

    out_df = pd.merge(out_df, geo_xwalk, left_on="dest_taz", right_on="taz", how="left")
    out_df = pd.merge(out_df, msa_names, on="msa", how="left")
    out_df.rename(columns={"msa_name": "dest_msa"}, inplace=True)
    out_df.drop(["mgra", "taz", "msa"], axis=1, inplace=True)

    return out_df


def get_omx_cores_as_df(omx_file_name, core_names, out_cols):
    omx_matrix = omx.open_file(omx_file_name)
    dim = omx_matrix.shape()[0]

    out_df = pd.DataFrame(
        {"orig_taz": np.repeat(1 + np.arange(dim), dim), "dest_taz": np.tile(1 + np.arange(dim), dim)})

    for i, core in enumerate(core_names):
        core_mtx = omx_matrix[core]
        col = out_cols[i]
        values_mat = np.array(core_mtx)
        values_df = pd.DataFrame({col: np.reshape(values_mat, (dim ** 2))})

        out_df = pd.concat([out_df, values_df], axis=1)

    return out_df


def get_msa_employment(mgra_input, emp_external_gateways, geo_xwalk, msa_names, year):
    input_df = mgra_input[["mgra", "emp_fed_non_mil", "emp_fed_mil", "emp_total"]]
    emp_df = input_df.copy()
    emp_df.loc[:, "emp_non_fed"] = (emp_df["emp_total"]
                                    - emp_df["emp_fed_non_mil"]
                                    - emp_df["emp_fed_mil"])

    emp_df = pd.merge(emp_df, geo_xwalk, how="left", on="mgra")
    emp_df = pd.merge(emp_df, msa_names, how="left", on="msa")

    emp_external_df = emp_external_gateways[emp_external_gateways["Year"] == year]

    external_msa_df = pd.DataFrame(
        {
            "County": ["IM", "LA", "OR", "RV", "SB"],
            "msa_name": ["IMPERIAL COUNTY", "LOS ANGELES COUNTY", "ORANGE COUNTY", "RIVERSIDE COUNTY", "SAN BERNARDINO COUNTY"]
        }
    )

    emp_external_df = pd.merge(emp_external_df, external_msa_df, how="left", on="County")
    emp_external_df["msa_name"] = emp_external_df["msa_name"].fillna("OTHER")

    emp_external_df.loc[:, "emp_non_fed"] = emp_external_df["Employment"]
    emp_external_df.loc[:, "emp_fed_non_mil"] = emp_external_df["Employment"]
    emp_external_df.loc[:, "emp_fed_mil"] = emp_external_df["Employment"]

    emp_external_df = emp_external_df[["msa_name", "emp_non_fed", "emp_fed_non_mil", "emp_fed_mil"]]

    out_df = emp_df.groupby(["msa_name"]).agg(
        {"emp_non_fed": 'sum', "emp_fed_non_mil": 'sum', "emp_fed_mil": 'sum'}
    )
    out_df.reset_index(inplace=True)

    out_df = pd.concat([out_df, emp_external_df], axis=0)

    return out_df


def get_emp_and_pop_by_corridor(mgra_input, corridor_mgra):
    emp_pop_df = mgra_input[["mgra", "emp_total", "pop"]]

    work_df = pd.merge(corridor_mgra, emp_pop_df, on="mgra", how="left")

    work_df["emp_corridor"] = work_df["emp_total"] * work_df["weight"]
    work_df["pop_corridor"] = work_df["pop"] * work_df["weight"]

    out_df = work_df.groupby(["corridor"]).agg({"emp_corridor": 'sum', "pop_corridor": 'sum'})
    out_df.reset_index(inplace=True)
    out_df["emp_corridor"] = out_df["emp_corridor"].astype(int)
    out_df["emp_pct_corridor"] = out_df["emp_corridor"] / emp_pop_df["emp_total"].sum()
    out_df["pop_corridor"] = out_df["pop_corridor"].astype(int)
    out_df["pop_pct_corridor"] = out_df["pop_corridor"] / emp_pop_df["pop"].sum()

    return out_df


def get_emission_factors(input_emission_data, scen_year):
    # prepare a dataframe of emission factors for the scen year

    emission_df = input_emission_data[input_emission_data["Year"] == scen_year]
    emission_df = emission_df[emission_df["Vehicle Type"] == "Passenger Car"]
    emission_df.reset_index(inplace=True, drop=True)

    co2_runex_emission_factor = emission_df["CO2 RunEx Emission Factor (gr/mile)"][0]
    co2_strex_emission_factor = emission_df["CO2 StrEx Emission Factor (gr/trip)"][0]

    out_data = {
        "Variable": ["Year",
                     "CO2 RunEx Emission Factor (gr/mile)",
                     "CO2 StrEx Emission Factor (gr/trip)"],
        "Value": [scen_year,
                  co2_runex_emission_factor,
                  co2_strex_emission_factor]
    }

    out_df = pd.DataFrame(out_data)

    return out_df


def compute_emission_results(vanpool_demand, vanpool_mean_stats, emission_factors, regional_population):
    # compute emission results
    co2_runex_emission_factor = emission_factors.loc[
        emission_factors.Variable == "CO2 RunEx Emission Factor (gr/mile)", "Value"].values[0]
    co2_strex_emission_factor = emission_factors.loc[
        emission_factors.Variable == "CO2 StrEx Emission Factor (gr/trip)", "Value"].values[0]

    daily_trips_reduction = 2 * np.sum(
        vanpool_mean_stats["Avg_Passenger_ExDriver"] * vanpool_demand["Num_Vans_Total_Scen"])
    daily_vmt_reduction = np.sum(
        vanpool_mean_stats["Avg_Passenger_ExDriver"] * vanpool_demand["Num_Vans_Total_Scen"] * vanpool_mean_stats[
            "Avg_RoundTrip_Mileage"])
    incounty_vmt_reduction = np.sum(
        vanpool_mean_stats["Avg_Passenger_ExDriver"] * vanpool_demand["Num_Vans_Total_Scen"] * vanpool_mean_stats[
            "Avg_InCounty_RoundTrip_Mileage"])

    GRAMS_TO_SHORT_TONS = 0.0000011

    coldstart_ghg_reduction_tons = (daily_trips_reduction
                                    * co2_strex_emission_factor
                                    * GRAMS_TO_SHORT_TONS)
    vmt_ghg_reduction_tons = (incounty_vmt_reduction
                              * co2_runex_emission_factor
                              * GRAMS_TO_SHORT_TONS)
    daily_total_ghg_reduction_tons = coldstart_ghg_reduction_tons + vmt_ghg_reduction_tons
    daily_ghg_reduction_percapita_lbs = daily_total_ghg_reduction_tons * 2000 / regional_population

    out_data = {
        "Variable": ["Regional Population",
                     "Total daily vehicle trip reduction",
                     "Total daily VMT reduction by vanpooling",
                     "VMT reduced in San Diego County by vanpooling",
                     "GHG reduction due to cold starts (short tons)",
                     "GHG reduction due to VMT (short tons)",
                     "Daily Total GHG reduction (short tons)",
                     "Daily Per capita GHG reduction (lbs/person)"],
        "Value": [regional_population,
                  daily_trips_reduction,
                  daily_vmt_reduction,
                  incounty_vmt_reduction,
                  coldstart_ghg_reduction_tons,
                  vmt_ghg_reduction_tons,
                  daily_total_ghg_reduction_tons,
                  daily_ghg_reduction_percapita_lbs]
    }

    out_df = pd.DataFrame(out_data)

    return out_df


def compute_corridor_reductions(regional_results, corridor_emp_pop_data, emission_factors, scenario_year):
    # apportion regional vmt reduction to corridors based on employment percent

    co2_runex_emission_factor = emission_factors.loc[
        emission_factors.Variable == "CO2 RunEx Emission Factor (gr/mile)", "Value"].values[0]

    co2_strex_emission_factor = emission_factors.loc[
        emission_factors.Variable == "CO2 StrEx Emission Factor (gr/trip)", "Value"].values[0]

    GRAMS_TO_SHORT_TONS = 0.0000011

    total_trip_reduction = regional_results.loc[
        regional_results.Variable == "Total daily vehicle trip reduction", "Value"].values[0]

    total_incounty_vmt_reduction = regional_results.loc[
        regional_results.Variable == "VMT reduced in San Diego County by vanpooling", "Value"].values[0]

    work_df = corridor_emp_pop_data.copy()

    work_df["Scenario Year"] = scenario_year

    work_df["Total daily vehicle trip reduction"] = work_df["emp_pct_corridor"] * total_trip_reduction

    work_df["Total daily VMT reduction"] = work_df["emp_pct_corridor"] * total_incounty_vmt_reduction

    work_df["GHG reduction due to cold starts (short tons)"] = (work_df["Total daily vehicle trip reduction"]
                                                                * co2_strex_emission_factor
                                                                * GRAMS_TO_SHORT_TONS)

    work_df["GHG reduction due to VMT (short tons)"] = (work_df["Total daily VMT reduction"]
                                                        * co2_runex_emission_factor
                                                        * GRAMS_TO_SHORT_TONS)

    work_df["Daily Total GHG reduction (short tons)"] = (work_df["GHG reduction due to cold starts (short tons)"]
                                                         + work_df["GHG reduction due to VMT (short tons)"])

    work_df.drop(
        ["emp_corridor", "pop_corridor", "emp_pct_corridor", "pop_pct_corridor"],
        axis=1,
        inplace=True
    )
    work_df.set_index('corridor', inplace=True)

    out_df = work_df.transpose()
    out_df.reset_index(inplace=True)
    out_df.rename(columns={"index": "Variable"}, inplace=True)

    return out_df


def write_results(results_dict, out_file_name, out_dir):
    with pd.ExcelWriter(os.path.join(out_dir, out_file_name)) as writer:
        for key, value in results_dict.items():
            value.to_excel(writer, sheet_name=key, index=False)


if __name__ == "__main__":
    args = sys.argv
    config_filename = args[1]
    # print(config_filename)

    if not os.path.exists(config_filename):
        msg = "Configuration file doesn't exist at: {}".format(config_filename)
        raise ValueError(msg)

    with open(config_filename, "r") as yml_file:
        config = yaml.safe_load(yml_file)

    # read config: inputs
    vanpool_od_file = config['inputs']['vanpool_od_file']
    mgra_scen_input_file = config['inputs']['mgra_scen_input_file']
    mgra_base_input_file = config['inputs']['mgra_base_input_file']
    geography_xwalk_file = config['inputs']['geography_xwalk_file']
    msa_names_file = config['inputs']['msa_names_file']
    employment_forecast_scag_file = config['inputs']['employment_forecast_scag_file']
    skim_base_file = config['inputs']['skim_base_file']
    skim_scen_file = config['inputs']['skim_scen_file']
    emission_factors_file = config['inputs']['emission_factors_file']
    zipcode_coordinates_file = config['inputs']['zipcode_coordinates_file']
    external_gateways_file = config['inputs']['external_gateways_file']
    corridor_mgra_xwalk_file = config['inputs']['corridors_mgra_xwalk_file']
    individual_tours_file = config['inputs']['individual_tours_output_file']

    # read config: outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # read config: parameters
    base_year = config['parameters']['base_year']
    scen_year = config['parameters']['scen_year']
    c_ivt = config['parameters']['c_ivt']
    avg_vanpool_occupancy = config['parameters']['avg_vanpool_occupancy']
    military_base_taz = config['parameters']['military_base_taz']
    pct_work_trips_over_50mi = config['parameters']['pct_work_trips_over_50mi']
    sov_time_core_name = config['parameters']['sov_am_time_core_name']
    hov_time_core_name = config['parameters']['hov_am_time_core_name']
    abm_version = config['parameters']['abm_version']

    # read data
    vanpool_od_df = pd.read_csv(vanpool_od_file)
    mgra_scen_input_df = pd.read_csv(mgra_scen_input_file)
    mgra_base_input_df = pd.read_csv(mgra_base_input_file)
    geo_xwalk_df = pd.read_csv(geography_xwalk_file)
    msa_names_df = pd.read_csv(msa_names_file)
    emp_forecast_scag_df = pd.read_csv(employment_forecast_scag_file)
    emission_df = pd.read_excel(emission_factors_file)
    zipcode_coordinates_df = pd.read_csv(zipcode_coordinates_file)
    external_gateways_df = pd.read_csv(external_gateways_file)
    corridor_mgra_df = pd.read_csv(corridor_mgra_xwalk_file)
    indiv_tours_df = pd.read_csv(individual_tours_file)

    # format data columns
    zipcode_coordinates_df.columns = zipcode_coordinates_df.columns.str.replace(' ', '_')
    external_gateways_df.columns = external_gateways_df.columns.str.replace(' ', '_')
    vanpool_od_df.columns = vanpool_od_df.columns.str.replace(' ', '_')
    emp_forecast_scag_df.columns = emp_forecast_scag_df.columns.str.replace(' ', '_')

    vanpool_od_df["Origin"] = vanpool_od_df["Origin"].str.lower()
    vanpool_od_df["Destination"] = vanpool_od_df["Destination"].str.lower()
    external_gateways_df["Home_County"] = external_gateways_df["Home_County"].str.lower()

    # join origin zip code
    data_df = pd.merge(vanpool_od_df, zipcode_coordinates_df,
                       left_on="Origin_Zip", right_on="Zip_Code", how="left")
    data_df.rename(
        columns={"Lat": "Origin_Lat",
                 "Long": "Origin_Long",
                 "Lat_Feet": "Origin_Lat_Feet",
                 "Long_Feet": "Origin_Long_Feet"
                 },
        inplace=True
    )
    data_df.drop(columns=["Zip_Code"], inplace=True)

    # join destination zip code
    data_df = pd.merge(data_df, zipcode_coordinates_df,
                       left_on="Destination_Zip", right_on="Zip_Code", how="left")
    data_df.rename(
        columns={"Lat": "Destination_Lat",
                 "Long": "Destination_Long",
                 "Lat_Feet": "Destination_Lat_Feet",
                 "Long_Feet": "Destination_Long_Feet"
                 },
        inplace=True
    )
    data_df.drop(columns=["Zip_Code"], inplace=True)

    # join external gateways
    data_df = pd.merge(
        data_df,
        external_gateways_df,
        left_on="Origin",
        right_on="Home_County",
        how="left"
    )
    data_df.Gateway_Latitude = data_df.Gateway_Latitude.fillna(data_df.Origin_Lat_Feet)
    data_df.Gateway_Longitude = data_df.Gateway_Longitude.fillna(data_df.Origin_Long_Feet)
    data_df.drop(columns=["Home_County"], inplace=True)

    # calculate van miles traveled
    data_df["OutofCounty_Oneway_Mileage"] = np.sqrt(
        (data_df["Origin_Lat_Feet"] - data_df["Gateway_Latitude"]) ** 2 +
        (data_df["Origin_Long_Feet"] - data_df["Gateway_Longitude"]) ** 2
    ) / 5280

    data_df["InCounty_Roundtrip_Mileage"] = np.maximum(
        data_df["Daily_Round_Trip_Mileage"] - 2 * data_df["OutofCounty_Oneway_Mileage"], 0)

    # calculate total vans by industry type, origin and destination
    vanpool_df = data_df.groupby(["Industry_Type", "Origin", "Destination"])[
        "Van_ID"].count().to_frame()
    vanpool_df.reset_index(inplace=True)
    vanpool_df.rename(columns={"Van_ID": "Num_Vans_Base"}, inplace=True)

    # get total vans by destination
    vanpool_total_by_dest_df = data_df.groupby(["Industry_Type", "Destination"])[
        "Van_ID"].count().to_frame()
    vanpool_total_by_dest_df.reset_index(inplace=True)
    vanpool_total_by_dest_df.rename(columns={"Van_ID": "Num_Vans_By_Dest_Base"}, inplace=True)

    # scen year employment
    emp_scen_by_msa_df = get_msa_employment(
        mgra_scen_input_df,
        emp_forecast_scag_df,
        geo_xwalk_df,
        msa_names_df,
        scen_year
    )
    emp_scen_by_msa_df.columns = ["Name", "Emp_Non_Fed_Scen",
                                  "Emp_Fed_Non_Mil_Scen", "Emp_Fed_Mil_Scen"]

    # base year employment
    emp_base_by_msa_df = get_msa_employment(
        mgra_base_input_df,
        emp_forecast_scag_df,
        geo_xwalk_df,
        msa_names_df,
        base_year
    )
    emp_base_by_msa_df.columns = ["Name", "Emp_Non_Fed_Base",
                                  "Emp_Fed_Non_Mil_Base", "Emp_Fed_Mil_Base"]

    # calculate employment growth factor by msa and type
    emp_by_msa_df = pd.merge(emp_scen_by_msa_df, emp_base_by_msa_df, on="Name", how="left")
    emp_by_msa_df["Emp_Non_Fed_Growth"] = (
        emp_by_msa_df["Emp_Non_Fed_Scen"] / emp_by_msa_df["Emp_Non_Fed_Base"]
    )
    emp_by_msa_df["Emp_Fed_Non_Mil_Growth"] = (
        emp_by_msa_df["Emp_Fed_Non_Mil_Scen"] / emp_by_msa_df["Emp_Fed_Non_Mil_Base"]
    )
    emp_by_msa_df["Emp_Fed_Mil_Growth"] = (
        emp_by_msa_df["Emp_Fed_Mil_Scen"] / emp_by_msa_df["Emp_Fed_Mil_Base"]
    )

    emp_by_msa_df["Emp_Fed_Non_Mil_Growth"] = emp_by_msa_df["Emp_Fed_Non_Mil_Growth"].fillna(1)
    emp_by_msa_df["Emp_Fed_Mil_Growth"] = emp_by_msa_df["Emp_Fed_Mil_Growth"].fillna(1)
    emp_by_msa_df["Emp_Non_Fed_Growth"] = emp_by_msa_df["Emp_Non_Fed_Growth"].fillna(1)
    emp_by_msa_df["Name"] = emp_by_msa_df["Name"].str.lower()

    emp_growth_by_msa_df = emp_by_msa_df[[
        "Name", "Emp_Non_Fed_Growth", "Emp_Fed_Non_Mil_Growth", "Emp_Fed_Mil_Growth"]]
    emp_growth_by_msa_df.columns = ["Name", "Non-Federal", "Federal", "Military"]
    emp_growth_by_msa_df = pd.melt(emp_growth_by_msa_df,
                                   id_vars=["Name"],
                                   value_vars=["Non-Federal", "Federal", "Military"],
                                   var_name="Type",
                                   value_name="Emp_Growth"
                                   )

    # calculate vanpool in scenario year based on employment growth
    vanpool_total_by_dest_df = pd.merge(
        vanpool_total_by_dest_df,
        emp_growth_by_msa_df,
        left_on=["Destination", "Industry_Type"],
        right_on=["Name", "Type"],
        how="left"
    )
    vanpool_total_by_dest_df["Emp_Growth"] = vanpool_total_by_dest_df["Emp_Growth"].fillna(0)
    vanpool_total_by_dest_df["Num_Vans_By_Dest_Scen"] = (vanpool_total_by_dest_df["Num_Vans_By_Dest_Base"]
                                                         * vanpool_total_by_dest_df["Emp_Growth"])
    vanpool_total_by_dest_df["Num_Vans_By_Dest_Scen"] = vanpool_total_by_dest_df["Num_Vans_By_Dest_Scen"].astype(
        int)

    vanpool_df = pd.merge(
        vanpool_df,
        vanpool_total_by_dest_df,
        on=["Industry_Type", "Destination"],
        how="left"
    )

    vanpool_df["Num_Vans_Emp_Scen"] = vanpool_df["Num_Vans_Base"] * vanpool_df["Num_Vans_By_Dest_Scen"] / vanpool_df[
        "Num_Vans_By_Dest_Base"]
    vanpool_df["Num_Vans_Emp_Scen"] = vanpool_df["Num_Vans_Emp_Scen"].astype(int)
    vanpool_df = vanpool_df[["Type", "Origin", "Destination", "Num_Vans_Base", "Num_Vans_Emp_Scen"]]

    # parameter from vanpool demand
    weekday_vanpool_demand = 2 * vanpool_od_df["Vehicle_Capacity"].sum() * avg_vanpool_occupancy

    # work trips from the model
    if abm_version == 'ABM2':
        work_tours_query = "tour_purpose == 'Work' and (tour_mode == 1 or tour_mode == 2)"
    elif abm_version == 'ABM2+':
        work_tours_query = "tour_purpose == 'Work' and tour_mode == 1"
    else:
        msg = "ABM version is not correct in the config file. Specify it as either 'ABM2' or 'ABM2+'."
        raise ValueError(msg)

    work_tours_df = indiv_tours_df.query(work_tours_query, engine="python")
    num_work_trips = work_tours_df.shape[0]

    # pct_work_trips_over_50mi is percent of daily work trips with one-way distance of 50 miles or more
    potential_weekday_vanpool_demand = num_work_trips * (pct_work_trips_over_50mi/100) * 2

    prob_vanpool = round(weekday_vanpool_demand / potential_weekday_vanpool_demand, 4)

    # get msa to msa (average) travel times for base and scenario years
    travel_times_msa_df = get_msa_travel_times(
        skim_base_file,
        skim_scen_file,
        military_base_taz,
        geo_xwalk_df,
        msa_names_df,
        external_gateways_df,
        sov_time_core_name,
        hov_time_core_name
    )

    # calculate elasticity and growth in vanpool demand because of travel time savings
    elasticity_df = travel_times_msa_df.copy()
    elasticity_df["Mil_Elasticity"] = (c_ivt
                                       * elasticity_df["sov_time_base_mil"]
                                       * (1 - prob_vanpool))
    elasticity_df["NonMil_Elasticity"] = (c_ivt
                                          * elasticity_df["sov_time_base"]
                                          * (1 - prob_vanpool))

    elasticity_df["Mil_Growth"] = np.where(
        elasticity_df["sov_time_base_mil"] > 0,
        1 + elasticity_df["Mil_Elasticity"] * elasticity_df["time_savings_scen_mil"] / elasticity_df[
            "sov_time_base_mil"],
        1
    )

    elasticity_df["NonMil_Growth"] = np.where(
        elasticity_df["sov_time_base"] > 0,
        1 + elasticity_df["NonMil_Elasticity"] * elasticity_df["time_savings_scen"] / elasticity_df[
            "sov_time_base"],
        1
    )

    elasticity_df["orig_msa"] = elasticity_df["orig_msa"].str.lower()
    elasticity_df["dest_msa"] = elasticity_df["dest_msa"].str.lower()

    # calculate vanpool in scenario year based on travel time savings
    vanpool_df = pd.merge(
        vanpool_df,
        elasticity_df[["orig_msa", "dest_msa", "Mil_Growth", "NonMil_Growth"]],
        left_on=["Origin", "Destination"],
        right_on=["orig_msa", "dest_msa"],
        how="left"
    )

    vanpool_df["Num_Vans_ML_Scen"] = np.where(
        vanpool_df["Type"] == "Military",
        vanpool_df["Num_Vans_Emp_Scen"] * (vanpool_df["Mil_Growth"] - 1),
        vanpool_df["Num_Vans_Emp_Scen"] * (vanpool_df["NonMil_Growth"] - 1)
    )
    vanpool_df["Num_Vans_ML_Scen"].fillna(0, inplace=True)
    vanpool_df["Num_Vans_ML_Scen"] = np.where(
        vanpool_df["Num_Vans_ML_Scen"] < 0,
        0,
        vanpool_df["Num_Vans_ML_Scen"])

    # total vanpool for the scenario
    vanpool_df["Num_Vans_Total_Scen"] = (vanpool_df["Num_Vans_Emp_Scen"]
                                         + vanpool_df["Num_Vans_ML_Scen"])

    vanpool_df["Num_Vans_Emp_Scen"] = vanpool_df["Num_Vans_Emp_Scen"].astype(int)
    vanpool_df["Num_Vans_ML_Scen"] = vanpool_df["Num_Vans_ML_Scen"].astype(int)
    vanpool_df["Num_Vans_Total_Scen"] = vanpool_df["Num_Vans_Total_Scen"].astype(int)

    # vanpool demand by type
    vanpool_demand_df = vanpool_df.groupby(["Type"]).agg(
        {'Num_Vans_Base': 'sum',
         'Num_Vans_Emp_Scen': 'sum',
         'Num_Vans_ML_Scen': 'sum',
         'Num_Vans_Total_Scen': 'sum'
         }
    )
    vanpool_demand_df.reset_index(inplace=True)

    # Calculations for required mean stats for vanpool od data
    vanpool_mean_stat_df = data_df.groupby(["Industry_Type"])[
        ["Vehicle_Capacity", "Daily_Round_Trip_Mileage", "InCounty_Roundtrip_Mileage"]].mean()
    vanpool_mean_stat_df.reset_index(inplace=True)
    vanpool_mean_stat_df.rename(
        columns={"Vehicle_Capacity": "Avg_Vanpool_Capacity",
                 "Daily_Round_Trip_Mileage": "Avg_RoundTrip_Mileage",
                 "InCounty_Roundtrip_Mileage": "Avg_InCounty_RoundTrip_Mileage"
                 },
        inplace=True
    )
    vanpool_mean_stat_df["Avg_Vanpool_Capacity_ExDriver"] = vanpool_mean_stat_df["Avg_Vanpool_Capacity"] - 1
    vanpool_mean_stat_df["Avg_Passenger_ExDriver"] = (vanpool_mean_stat_df["Avg_Vanpool_Capacity_ExDriver"]
                                                      * avg_vanpool_occupancy)

    emission_factors_df = get_emission_factors(emission_df, scen_year)

    regional_population_scen = mgra_scen_input_df["pop"].sum()

    regional_results_df = compute_emission_results(
        vanpool_demand_df,
        vanpool_mean_stat_df,
        emission_factors_df,
        regional_population_scen
    )

    corridor_emp_pop_df = get_emp_and_pop_by_corridor(mgra_scen_input_df, corridor_mgra_df)

    corridor_results_df = compute_corridor_reductions(
        regional_results_df,
        corridor_emp_pop_df,
        emission_factors_df,
        scen_year
    )

    vanpool_demand_df = vanpool_demand_df.append(
        {
            'Type': 'TOTAL',
            'Num_Vans_Base': vanpool_demand_df["Num_Vans_Base"].sum(),
            'Num_Vans_Emp_Scen': vanpool_demand_df["Num_Vans_Emp_Scen"].sum(),
            'Num_Vans_ML_Scen': vanpool_demand_df["Num_Vans_ML_Scen"].sum(),
            'Num_Vans_Total_Scen': vanpool_demand_df["Num_Vans_Total_Scen"].sum()
        },
        ignore_index=True
    )

    vanpool_demand_df.rename(
        columns={"Num_Vans_Base": "Num_Vans_" + str(base_year),
                 "Num_Vans_Emp_Scen": "Num_Vans_Emp_" + str(scen_year),
                 "Num_Vans_ML_Scen": "Num_Vans_ML_" + str(scen_year),
                 "Num_Vans_Total_Scen": "Num_Vans_Total_" + str(scen_year)
                 },
        inplace=True
    )

    results_dict = {"Regional_Results": regional_results_df,
                    "Corridor_Results": corridor_results_df,
                    "Vanpool_Growth": vanpool_demand_df,
                    "Emission_Factors": emission_factors_df}

    write_results(results_dict, output_results_filename, output_dir)
