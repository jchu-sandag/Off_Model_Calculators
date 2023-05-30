# Import Libraries
import geopandas as gpd
import pandas as pd
import numpy as np
import yaml
import sys
import openmatrix as omx
import os


def get_target_population_trips(trips_df, income_threshold, person_type_file, trip_purpose_file, origin_geography_file=None, destination_geography_file=None):
    """
    returns the target population trips based on different attributes
    i.e., person type, trip type, income threshold, orig MGRAs, dest MGRAs

    Args:
        trip_df: dataframe with all the individual trips
        income_threshold: income threshold value to identify low income population
        person_type_file: input file name that specifies the person types to include for target population
        trip_purpose_file: input file name that specifies the trip purpose types to include for target population
        origin_geography_file: input file name that specifies the origin mgra of the trips for target population
        destination_geography_file: input file name that specifies the destination mgra of the trips for target population

    Returns:
        dataframe with trips for target population

    """

    work_df = trips_df.copy()

    if income_threshold != None:
        work_df = work_df[work_df["income"] <= income_threshold]

    if person_type_file != None:
        person_type_target_df = pd.read_csv(person_type_file)
        person_type_target_df["target_population"] = person_type_target_df["target_population"].str.lower()
        target_person_types = person_type_target_df[
            person_type_target_df["target_population"] == "yes"]["person_type"].to_list()

        work_df = work_df[work_df["type"].isin(target_person_types)]

    if trip_purpose_file != None:
        trip_purpose_target_df = pd.read_csv(trip_purpose_file)
        trip_purpose_target_df["target_population"] = trip_purpose_target_df[
            "target_population"].str.lower()
        target_trip_purposes = trip_purpose_target_df[
            trip_purpose_target_df["target_population"] == "yes"]["tour_purpose"].to_list()

        work_df = work_df[work_df["tour_purpose"].isin(target_trip_purposes)]

    if origin_geography_file != None:
        orig_geo_df = pd.read_csv(origin_geography_file)
        target_orig_mgra = orig_geo_df["orig_mgra"].to_list()

        work_df = work_df[work_df["orig_mgra"].isin(target_orig_mgra)]

    if destination_geography_file != None:
        dest_geo_df = pd.read_csv(destination_geography_file)
        target_dest_mgra = dest_geo_df["dest_mgra"].to_list()

        work_df = work_df[work_df["dest_mgra"].isin(target_dest_mgra)]

    return work_df


def get_omx_cores_as_df(omx_file_name, core_names, out_cols):
    """
    get the skim cores from omx files as dataframe

    Args:
        omx_file_name: omx file name
        core_names: matrix core names to read from the omx file
        out_cols: desired output column names in the output dataframe

    Returns:
        dataframe with the desired matrix cores as columns

    """

    omx_matrix = omx.open_file(omx_file_name)
    dim = omx_matrix.shape()[0]

    out_df = pd.DataFrame({
        "orig_taz": np.repeat(1 + np.arange(dim), dim),
        "dest_taz": np.tile(1 + np.arange(dim), dim)
    })

    for i, core in enumerate(core_names):
        core_mtx = omx_matrix[core]
        col = out_cols[i]
        values_mat = np.array(core_mtx)
        values_df = pd.DataFrame({col: np.reshape(values_mat, (dim ** 2))})

        out_df = pd.concat([out_df, values_df], axis=1)

    return out_df


def get_trip_distance(trips_df, geography_xwalk_df, skim_file, sov_core_name, hov_core_name):
    """
    get the orig-dest trip distance for the trips

    Args:
        trips_df: dataframe with person trips
        geography_xwalk_df: dataframe for mgra to taz cross-walk
        skim_file: skim matrix file
        sov_core_name: core name for sov distance
        hov_core_name: core name for hov distance
    Returns:
        dataframe with trip distance appened to the input dataframe

    """

    trip_dist_df = get_omx_cores_as_df(
        skim_file,
        [sov_core_name, hov_core_name],
        ["sov_distance", "hov_distance"]
    )

    mgra_taz_df = geography_xwalk_df[["mgra", "taz"]]

    work_df = trips_df.copy()

    work_df = pd.merge(work_df, mgra_taz_df, left_on="orig_mgra", right_on="mgra", how="left")
    work_df.rename(columns={"taz": "orig_taz"}, inplace=True)
    work_df.drop(["mgra"], axis=1, inplace=True)

    work_df = pd.merge(work_df, mgra_taz_df, left_on="dest_mgra", right_on="mgra", how="left")
    work_df.rename(columns={"taz": "dest_taz"}, inplace=True)
    work_df.drop(["mgra"], axis=1, inplace=True)

    out_df = pd.merge(work_df, trip_dist_df, on=["orig_taz", "dest_taz"], how="left")

    return out_df


def calculate_vmt_reduction(trips_df, transit_modes, auto_modes, fare_ridership_elasticity, percent_fare_subsidy, auto_to_transit_switch_factor):
    """
    calculate the vmt reduction due to increase in transit ridership

    Args:
        trips_df: dataframe with trips for target population
        transit_modes: list of transit mode numbers
        auto_modes: list of auto mode numbers
        fare_ridership_elasticity: fare elasticity of transit ridership
        percent_fare_subsidy: percent of fare subsidy
        auto_to_transit_switch_factor: factor representing how many vehicle trips will be replaced by transit

    Returns:
        daily vmt reduction for the target population

    """

    transit_trips_df = trips_df[trips_df['trip_mode'].isin(transit_modes)]
    total_transit_trips = transit_trips_df.shape[0]

    auto_trips_df = trips_df.copy()
    auto_trips_df = auto_trips_df[auto_trips_df['trip_mode'].isin(auto_modes)]
    auto_trips_df['trip_distance'] = np.where(
        auto_trips_df['trip_mode'] == 1, auto_trips_df['sov_distance'], auto_trips_df['hov_distance'])
    average_auto_trips_length = np.mean(auto_trips_df['trip_distance'])

    vmt_reduction = (total_transit_trips
                     * average_auto_trips_length
                     * fare_ridership_elasticity
                     * (percent_fare_subsidy/100.0)
                     * auto_to_transit_switch_factor)

    vmt_reduction = np.absolute(vmt_reduction)

    return vmt_reduction


def get_emission_factors(input_emission_data, scen_year):
    """
    get a dataframe of emission factors for the scen year

    Args:
        input_emission_data: dataframe with all emission rates
        scen_year: scenario year
    Returns:
        dataframe with emission rates for the scen year

    """
    emission_df = input_emission_data.copy()
    emission_df = emission_df[emission_df["Year"] == scen_year]
    emission_df.reset_index(inplace=True, drop=True)

    return emission_df


def calculate_ghg_reduction(vmt_reduction, emission_factors, scen_year):
    """
    calculate total ghg reduction because of vmt reduction

    Args:
        vmt_reduction: daily vmt reduction for the target population
        emission_factors: emission factor dataframe
        scen_year: scenario year

    Returns:
        total ghg reductions

    """
    GRAMS_TO_SHORT_TONS = 0.0000011

    auto_runex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                  "Passenger Car"].reset_index().at[0, "CO2 RunEx Emission Factor (gr/mile)"]

    passenger_car_daily_ghg_reduction = (vmt_reduction
                                         * auto_runex_emission_factor
                                         * GRAMS_TO_SHORT_TONS)

    out_data = {
        "Variable": ["GHG reduction (short tons) due to vmt reduction for passenger cars"],
        "Value": [passenger_car_daily_ghg_reduction]
    }

    out_df = pd.DataFrame(out_data)

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
    individual_trip_file = config['inputs']['individual_trips_output_file']
    person_data_file = config['inputs']['person_data_file']
    household_data_file = config['inputs']['household_data_file']
    emission_factors_file = config['inputs']['emission_factors_file']
    skim_matrix_file = config['inputs']['skim_matrix_file']
    geography_xwalk_file = config['inputs']['geography_xwalk_file']

    origin_geography_file = config['inputs']['origin_geography_file']
    destination_geography_file = config['inputs']['destination_geography_file']
    person_type_target_population_file = config['inputs']['person_type_target_population_file']
    trip_purpose_target_population_file = config['inputs']['trip_purpose_target_population_file']

    # outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # parameters
    scen_year = config['parameters']['scen_year']
    income_threshold = config['parameters']['low_income_threshold']
    fare_ridership_elasticity = config['parameters']['fare_ridership_elasticity']
    percent_fare_subsidy = config['parameters']['percent_fare_subsidy']
    auto_to_transit_switch_factor = config['parameters']['auto_to_transit_switch_factor']
    auto_modes = config['parameters']['auto_modes']
    transit_modes = config['parameters']['transit_modes']
    sov_dist_core_name = config['parameters']['sov_ev_dist_core_name']
    hov_dist_core_name = config['parameters']['hov_ev_dist_core_name']

    # read inputs
    individual_trip_df = pd.read_csv(individual_trip_file, encoding='latin1')
    person_df = pd.read_csv(person_data_file)
    household_df = pd.read_csv(household_data_file)
    emission_df = pd.read_excel(emission_factors_file)
    geography_mapping_df = pd.read_csv(geography_xwalk_file)

    # join trip list, person data and household data
    individual_trip_df = individual_trip_df[[
        "hh_id", "person_id", "tour_purpose", "orig_mgra", "dest_mgra", "trip_mode"]]
    person_df = person_df[["hh_id", "person_id", "type"]]
    household_df = household_df[["hh_id", "income"]]

    data_df = pd.merge(
        individual_trip_df,
        person_df,
        how="left",
        on=["hh_id", "person_id"]
    )

    data_df = pd.merge(
        data_df,
        household_df,
        how="left",
        on=["hh_id"]
    )

    # get the target population trips
    target_population_trips_df = get_target_population_trips(
        data_df,
        income_threshold,
        person_type_target_population_file,
        trip_purpose_target_population_file,
        origin_geography_file,
        destination_geography_file
    )

    # add trip distance
    target_population_trips_df = get_trip_distance(
        target_population_trips_df,
        geography_mapping_df,
        skim_matrix_file,
        sov_dist_core_name,
        hov_dist_core_name
    )

    # calculate VMT change because of transit fare subsidy
    vmt_reduction = calculate_vmt_reduction(
        target_population_trips_df,
        transit_modes,
        auto_modes,
        fare_ridership_elasticity,
        percent_fare_subsidy,
        auto_to_transit_switch_factor
    )

    vmt_reduction_df = pd.DataFrame(
        {
            "Variable": ["Year", "VMT reduction from new transit trips due to transit subsidy"],
            "Value": [scen_year, vmt_reduction]
        }
    )

    # get emission factors
    emission_factors_df = get_emission_factors(emission_df, scen_year)

    # calculate ghg reductions
    ghg_reduction_df = calculate_ghg_reduction(vmt_reduction, emission_factors_df, scen_year)

    results_df = pd.concat([vmt_reduction_df, ghg_reduction_df])

    # writing results
    results_dict = {"GHG_Reduction": results_df, "Emission_Factors": emission_factors_df}
    write_results(results_dict, output_results_filename, output_dir)
