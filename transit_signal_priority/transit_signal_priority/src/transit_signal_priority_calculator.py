import geopandas as gpd
import pandas as pd
import numpy as np
import yaml
import sys
from utils import *


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
    highway_load_pm_file = config['inputs']['highway_load_pm_file']
    highway_network_shapefile = config['inputs']['highway_network_shapefile']
    transit_flow_file = config['inputs']['transit_flow_file']
    transit_aggflow_file = config['inputs']['transit_aggflow_file']
    transit_routes_file = config['inputs']['transit_routes_file']
    transit_stops_file = config['inputs']['transit_stops_file']
    transit_links_file = config['inputs']['transit_links_file']
    transit_onoff_file = config['inputs']['transit_onoff_file']
    tsp_strategy_file = config['inputs']['tsp_strategy_file']
    emission_factors_file = config['inputs']['emission_factors_file']

    # outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # parameters
    scen_year = config['parameters']['scen_year']
    periods_length = config['parameters']['periods_length_in_hr']
    AVG_DWELL_TIME_MIN = config['parameters']['avg_dwell_time_mins']
    GREEN_TO_CYCLE_WITHOUT_TSP = config['parameters']['green_to_cycle_without_tsp']
    GREEN_TO_CYCLE_PRIMARY_WITH_TSP = config['parameters']['green_to_cycle_primary_with_tsp']
    GREEN_TO_CYCLE_SECONDARY_WITH_TSP = config['parameters']['green_to_cycle_secondary_with_tsp']
    AVG_CYCLE_LENGTH = config['parameters']['avg_cycle_length_seconds']
    TT_ELASTICITY_TRANSIT_RIDERSHIP = config['parameters']['elasticity_transit_ridership']
    AVG_AUTO_OCCUPANCY = config['parameters']['avg_auto_occupancy']

    if not os.path.exists(tsp_strategy_file):
        msg = "TSP strategy file doesn't exist at: {}".format(tsp_strategy_file)
        raise ValueError(msg)

    with open(tsp_strategy_file, "r") as yml_file:
        tsp_config = yaml.safe_load(yml_file)
        tsp_strategy_dict = tsp_config['intersections']

    number_of_tsp_intersections = len(tsp_strategy_dict["intersection_node_id"])

    # read inputs
    links_df = gpd.read_file(highway_network_shapefile)
    highway_load_am_df = pd.read_csv(highway_load_am_file)
    highway_load_pm_df = pd.read_csv(highway_load_pm_file)
    transit_flow_df = pd.read_csv(transit_flow_file)
    transit_aggflow_df = pd.read_csv(transit_aggflow_file)
    transit_routes_df = pd.read_csv(transit_routes_file)
    transit_stops_df = pd.read_csv(transit_stops_file)
    transit_links_df = pd.read_csv(transit_links_file)
    transit_onoff_df = pd.read_csv(transit_onoff_file)
    emission_df = pd.read_excel(emission_factors_file)

    highway_load_data_dict = {"AM": highway_load_am_df, "PM": highway_load_pm_df}

    # get corridor links
    corridor_links_df = get_corridor_links(links_df, tsp_strategy_dict)

    query = "corridor == 'Primary' and intersection > 0"
    primary_links_df = corridor_links_df.query(query, engine="python")

    # calculate average transit headway for the intersection
    am_average_transit_headway = get_average_transit_headway(
        primary_links_df,
        transit_routes_df,
        transit_links_df,
        AVG_CYCLE_LENGTH,
        "AM"
    )
    pm_average_transit_headway = get_average_transit_headway(
        primary_links_df,
        transit_routes_df,
        transit_links_df,
        AVG_CYCLE_LENGTH,
        "PM"
    )

    # calculate average trip length of transit riders in the corridor
    am_avg_transit_trip_length = get_average_transit_trip_length(
        primary_links_df,
        transit_links_df,
        transit_flow_df,
        transit_onoff_df,
        "AM"
    )
    pm_avg_transit_trip_length = get_average_transit_trip_length(
        primary_links_df,
        transit_links_df,
        transit_flow_df,
        transit_onoff_df,
        "PM"
    )

    # calculate average transit flow at intersection
    am_avg_transit_flow = get_average_transit_flow(primary_links_df, transit_aggflow_df, "AM")
    pm_avg_transit_flow = get_average_transit_flow(primary_links_df, transit_aggflow_df, "PM")

    # calculate total stops along the primary corridor
    ab_transit_stops, ba_transit_stops = get_corridor_transit_stops(
        corridor_links_df,
        transit_routes_df,
        transit_stops_df
    )

    # get link attributes for all links (primary and secondary) in the corridor
    corridor_links_df = get_link_attributes(
        corridor_links_df,
        highway_load_data_dict,
        transit_aggflow_df
    )

    # calculate total travel time for the corridor by direction
    ab_am_time = round(corridor_links_df[corridor_links_df["corridor"] == "Primary"]["ab_am_time"].sum(
    ), 2) + (ab_transit_stops * AVG_DWELL_TIME_MIN)
    ba_am_time = round(corridor_links_df[corridor_links_df["corridor"] == "Primary"]["ba_am_time"].sum(
    ), 2) + (ba_transit_stops * AVG_DWELL_TIME_MIN)
    ab_pm_time = round(corridor_links_df[corridor_links_df["corridor"] == "Primary"]["ab_pm_time"].sum(
    ), 2) + (ab_transit_stops * AVG_DWELL_TIME_MIN)
    ba_pm_time = round(corridor_links_df[corridor_links_df["corridor"] == "Primary"]["ba_pm_time"].sum(
    ), 2) + (ba_transit_stops * AVG_DWELL_TIME_MIN)

    # prepare average intersection data
    data_df = corridor_links_df.copy()
    data_df = data_df.query("intersection > 0", engine="python")

    data_df["am_flow"] = np.where(data_df["approach"] == "AB",
                                  data_df["ab_am_flow"],
                                  data_df["ba_am_flow"])
    data_df["pm_flow"] = np.where(data_df["approach"] == "AB",
                                  data_df["ab_pm_flow"],
                                  data_df["ba_pm_flow"])

    data_df["am_auto_flow"] = np.where(data_df["approach"] == "AB",
                                       data_df["ab_am_auto_flow"],
                                       data_df["ba_am_auto_flow"])
    data_df["pm_auto_flow"] = np.where(data_df["approach"] == "AB",
                                       data_df["ab_pm_auto_flow"],
                                       data_df["ba_pm_auto_flow"])

    data_df["am_truck_flow"] = np.where(data_df["approach"] == "AB",
                                        data_df["ab_am_truck_flow"],
                                        data_df["ba_am_truck_flow"])
    data_df["pm_truck_flow"] = np.where(data_df["approach"] == "AB",
                                        data_df["ab_pm_truck_flow"],
                                        data_df["ba_pm_truck_flow"])

    data_df["am_capacity"] = np.where(data_df["approach"] == "AB",
                                      data_df["ab_am_capacity"],
                                      data_df["ba_am_capacity"])
    data_df["pm_capacity"] = np.where(data_df["approach"] == "AB",
                                      data_df["ab_pm_capacity"],
                                      data_df["ba_pm_capacity"])

    data_df["am_transit_flow"] = np.where(data_df["approach"] == "AB",
                                          data_df["ab_am_transit_flow"],
                                          data_df["ba_am_transit_flow"])
    data_df["pm_transit_flow"] = np.where(data_df["approach"] == "AB",
                                          data_df["ab_pm_transit_flow"],
                                          data_df["ba_pm_transit_flow"])

    voc_df = data_df.groupby(["intersection", "corridor"], as_index=False).agg(
        {"am_flow": 'sum', "pm_flow": 'sum', "am_capacity": 'sum', "pm_capacity": 'sum'})

    voc_df["am_voc"] = voc_df["am_flow"]/voc_df["am_capacity"]
    voc_df["pm_voc"] = voc_df["pm_flow"]/voc_df["pm_capacity"]
    voc_df = voc_df[["intersection", "corridor", "am_voc", "pm_voc"]]

    data_df = pd.merge(data_df, voc_df, how="left", on=["intersection", "corridor"])

    data_df = data_df[["hwycov_id", "from_node", "to_node", "link_name",
                       "corridor", "intersection", "approach", "am_auto_flow", "pm_auto_flow",
                       "am_truck_flow", "pm_truck_flow", "am_voc", "pm_voc", "am_transit_flow", "pm_transit_flow"]]

    data_intersection_df = data_df.groupby(["intersection", "corridor"], as_index=False).agg({
        "am_auto_flow": 'sum',
        "am_truck_flow": 'sum',
        "am_transit_flow": 'sum',
        "am_voc": 'mean',
        "pm_auto_flow": 'sum',
        "pm_truck_flow": 'sum',
        "pm_transit_flow": 'sum',
        "pm_voc": 'mean'
    })

    data_avg_intersection_df = data_intersection_df.groupby(["corridor"], as_index=False).agg({
        "am_auto_flow": 'mean',
        "am_truck_flow": 'mean',
        "am_transit_flow": 'mean',
        "am_voc": 'mean',
        "pm_auto_flow": 'mean',
        "pm_truck_flow": 'mean',
        "pm_transit_flow": 'mean',
        "pm_voc": 'mean'
    })

    data_avg_intersection_df.set_index('corridor', inplace=True)

    work_df = pd.DataFrame.from_dict({
        "type": ["auto", "truck", "bus"],
        "am_flow_primary": [data_avg_intersection_df.loc["Primary"].at["am_auto_flow"], data_avg_intersection_df.loc["Primary"].at["am_truck_flow"], data_avg_intersection_df.loc["Primary"].at["am_transit_flow"]],
        "am_flow_secondary": [data_avg_intersection_df.loc["Secondary"].at["am_auto_flow"], data_avg_intersection_df.loc["Secondary"].at["am_truck_flow"], 0],
        "am_voc_primary": [data_avg_intersection_df.loc["Primary"].at["am_voc"], data_avg_intersection_df.loc["Primary"].at["am_voc"], data_avg_intersection_df.loc["Primary"].at["am_voc"]],
        "am_voc_secondary": [data_avg_intersection_df.loc["Secondary"].at["am_voc"], data_avg_intersection_df.loc["Secondary"].at["am_voc"], data_avg_intersection_df.loc["Secondary"].at["am_voc"]],
        "pm_flow_primary": [data_avg_intersection_df.loc["Primary"].at["pm_auto_flow"], data_avg_intersection_df.loc["Primary"].at["pm_truck_flow"], data_avg_intersection_df.loc["Primary"].at["pm_transit_flow"]],
        "pm_flow_secondary": [data_avg_intersection_df.loc["Secondary"].at["pm_auto_flow"], data_avg_intersection_df.loc["Secondary"].at["pm_truck_flow"], 0],
        "pm_voc_primary": [data_avg_intersection_df.loc["Primary"].at["pm_voc"], data_avg_intersection_df.loc["Primary"].at["pm_voc"], data_avg_intersection_df.loc["Primary"].at["pm_voc"]],
        "pm_voc_secondary": [data_avg_intersection_df.loc["Secondary"].at["pm_voc"], data_avg_intersection_df.loc["Secondary"].at["pm_voc"], data_avg_intersection_df.loc["Secondary"].at["pm_voc"]]
    })

    work_df["am_voc_primary"] = work_df["am_voc_primary"].apply(lambda x: min(1, x))
    work_df["am_voc_secondary"] = work_df["am_voc_secondary"].apply(lambda x: min(1, x))
    work_df["pm_voc_primary"] = work_df["pm_voc_primary"].apply(lambda x: min(1, x))
    work_df["pm_voc_secondary"] = work_df["pm_voc_secondary"].apply(lambda x: min(1, x))

    am_prob_bus_in_cycle = AVG_CYCLE_LENGTH / (60 * am_average_transit_headway)
    pm_prob_bus_in_cycle = AVG_CYCLE_LENGTH / (60 * pm_average_transit_headway)

    work_df["am_primary_delay_without_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_WITHOUT_TSP)**2) / (1 - work_df["am_voc_primary"] * GREEN_TO_CYCLE_WITHOUT_TSP)
    work_df["am_secondary_delay_without_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_WITHOUT_TSP)**2) / (1 - work_df["am_voc_secondary"] * GREEN_TO_CYCLE_WITHOUT_TSP)
    work_df["am_primary_delay_with_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_PRIMARY_WITH_TSP)**2) / (1 - work_df["am_voc_primary"] * GREEN_TO_CYCLE_PRIMARY_WITH_TSP)
    work_df["am_secondary_delay_with_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_SECONDARY_WITH_TSP)**2) / (1 - work_df["am_voc_secondary"] * GREEN_TO_CYCLE_SECONDARY_WITH_TSP)

    work_df["pm_primary_delay_without_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_WITHOUT_TSP)**2) / (1 - work_df["pm_voc_primary"] * GREEN_TO_CYCLE_WITHOUT_TSP)
    work_df["pm_secondary_delay_without_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_WITHOUT_TSP)**2) / (1 - work_df["pm_voc_secondary"] * GREEN_TO_CYCLE_WITHOUT_TSP)
    work_df["pm_primary_delay_with_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_PRIMARY_WITH_TSP)**2) / (1 - work_df["pm_voc_primary"] * GREEN_TO_CYCLE_PRIMARY_WITH_TSP)
    work_df["pm_secondary_delay_with_tsp"] = 0.5 * (AVG_CYCLE_LENGTH / 60) * (
        (1 - GREEN_TO_CYCLE_SECONDARY_WITH_TSP)**2) / (1 - work_df["pm_voc_secondary"] * GREEN_TO_CYCLE_SECONDARY_WITH_TSP)

    work_df["am_delay_savings"] = (
        (work_df["am_flow_primary"]
         * (work_df["am_primary_delay_without_tsp"] - work_df["am_primary_delay_with_tsp"])
         * am_prob_bus_in_cycle)
        +
        (work_df["am_flow_secondary"]
         * (work_df["am_secondary_delay_without_tsp"] - work_df["am_primary_delay_with_tsp"])
         * am_prob_bus_in_cycle)
    ) / 3600

    work_df["pm_delay_savings"] = (
        (work_df["pm_flow_primary"]
         * (work_df["pm_primary_delay_without_tsp"] - work_df["pm_primary_delay_with_tsp"])
         * pm_prob_bus_in_cycle) +
        (work_df["pm_flow_secondary"]
         * (work_df["pm_secondary_delay_without_tsp"] - work_df["pm_primary_delay_with_tsp"])
         * pm_prob_bus_in_cycle)
    ) / 3600

    work_df["total_delay_savings"] = (work_df["am_delay_savings"] * periods_length["AM"]
                                      + work_df["pm_delay_savings"] * periods_length["PM"]) * number_of_tsp_intersections

    work_df["Variable"] = ["Delay Savings (hours) for Autos due to TSP",
                           "Delay Savings (hours) for Trucks due to TSP",
                           "Delay Savings (hours) for Buses due to TSP"]

    work_df["Value"] = work_df["total_delay_savings"]

    delay_savings_df = work_df[["type", "total_delay_savings", "Variable", "Value"]]

    # calculate reduction in VMT due to marginal mode shift to transit because of travel time reduction
    am_primary_delay_without_tsp = work_df.loc[0].at["am_primary_delay_without_tsp"]
    am_primary_delay_with_tsp = work_df.loc[0].at["am_primary_delay_with_tsp"]
    am_pct_change_bus_travel_time = (
        am_primary_delay_with_tsp - am_primary_delay_without_tsp) * number_of_tsp_intersections / ab_am_time

    pm_primary_delay_without_tsp = work_df.loc[0].at["pm_primary_delay_without_tsp"]
    pm_primary_delay_with_tsp = work_df.loc[0].at["pm_primary_delay_with_tsp"]
    pm_pct_change_bus_travel_time = (
        pm_primary_delay_with_tsp - pm_primary_delay_without_tsp) * number_of_tsp_intersections / ab_pm_time

    am_new_transit_riders = (am_pct_change_bus_travel_time
                             * TT_ELASTICITY_TRANSIT_RIDERSHIP
                             * am_avg_transit_flow
                             * periods_length["AM"])
    pm_new_transit_riders = (pm_pct_change_bus_travel_time
                             * TT_ELASTICITY_TRANSIT_RIDERSHIP
                             * pm_avg_transit_flow
                             * periods_length["PM"])

    am_vmt_savings = am_avg_transit_trip_length * (am_new_transit_riders / AVG_AUTO_OCCUPANCY)
    pm_vmt_savings = pm_avg_transit_trip_length * (pm_new_transit_riders / AVG_AUTO_OCCUPANCY)

    total_vmt_savings = am_vmt_savings + pm_vmt_savings

    vmt_savings_df = pd.DataFrame.from_dict(
        {
            'Variable': ["VMT AM savings due to mode shift to transit", "VMT PM savings due to mode shift to transit", "VMT TOTAL savings due to mode shift to transit"],
            'Value': [am_vmt_savings, pm_vmt_savings, total_vmt_savings]
        }
    )

    # get emission factors
    emission_factors_df = get_emission_factors(emission_df, scen_year)

    # calculate ghg reductions
    ghg_reduction_df = calculate_ghg_reduction(
        delay_savings_df,
        total_vmt_savings,
        emission_factors_df
    )

    delay_savings_df = delay_savings_df[["Variable", "Value"]]

    # writing results
    results_dict = {"GHG_Reduction": ghg_reduction_df,
                    "VMT_Reduction": vmt_savings_df,
                    "Delay_Reduction": delay_savings_df,
                    "Emission_Factors": emission_factors_df}

    write_results(results_dict, output_results_filename, output_dir)
