import geopandas as gpd
import pandas as pd
import numpy as np
import yaml
from utils import *
import sys

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
    highway_load_md_file = config['inputs']['highway_load_md_file']
    highway_load_ea_file = config['inputs']['highway_load_ea_file']
    highway_load_ev_file = config['inputs']['highway_load_ev_file']
    highway_network_shapefile = config['inputs']['highway_network_shapefile']
    atms_strategy_file = config['inputs']['atms_strategy_file']
    emission_factors_file = config['inputs']['emission_factors_file']
    intersection_delays_file = config['inputs']['intersection_delays_file']

    # outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # parameters
    scen_year = config['parameters']['scen_year']
    percent_atms_delay_reduction = config['parameters']['percent_atms_delay_reduction']
    apply_adaptive_signal_reduction = config['parameters']['apply_adaptive_signal_reduction']

    if apply_adaptive_signal_reduction:
        percent_adaptive_signal_reduction = config['parameters']['percent_adaptive_signal_reduction']
        total_percent_delay_reduction = percent_atms_delay_reduction + percent_adaptive_signal_reduction
    else:
        total_percent_delay_reduction = percent_atms_delay_reduction

    # other paramteters
    MAX_DELAY = 0.8

    if not os.path.exists(atms_strategy_file):
        msg = "ATMS strategy file doesn't exist at: {}".format(atms_strategy_file)
        raise ValueError(msg)

    with open(atms_strategy_file, "r") as yml_file:
        atms_config = yaml.safe_load(yml_file)
        atms_strategy_dict = atms_config['intersections']

    number_of_atms_intersections = len(atms_strategy_dict["intersection_node_id"])

    # read inputs
    links_df = gpd.read_file(highway_network_shapefile)
    highway_load_am_df = pd.read_csv(highway_load_am_file)
    highway_load_pm_df = pd.read_csv(highway_load_pm_file)
    highway_load_md_df = pd.read_csv(highway_load_md_file)
    highway_load_ea_df = pd.read_csv(highway_load_ea_file)
    highway_load_ev_df = pd.read_csv(highway_load_ev_file)
    emission_df = pd.read_excel(emission_factors_file)

    highway_load_data_dict = {"AM": highway_load_am_df, "PM": highway_load_pm_df,
                              "MD": highway_load_md_df, "EA": highway_load_ea_df, "EV": highway_load_ev_df}

    # intersection delays
    if (atms_strategy_dict["intersection_delays_in_minutes"]) == None:
        USER_SPECIFIED_DELAY = 0

        intersection_delays_df = pd.read_csv(intersection_delays_file)
        intersection_delays_df = intersection_delays_df[[
            "From", "To", "Inter_Time_EA",  "Inter_Time_AM",  "Inter_Time_MD", "Inter_Time_PM", "Inter_Time_EV"]]
        intersection_delays_df = intersection_delays_df.drop_duplicates()
        intersection_delays_df.rename(
            columns={
                "Inter_Time_EA": "ea_delay",
                "Inter_Time_AM": "am_delay",
                "Inter_Time_MD": "md_delay",
                "Inter_Time_PM": "pm_delay",
                "Inter_Time_EV": "ev_delay"
            },
            inplace=True
        )
    else:
        USER_SPECIFIED_DELAY = 1

        intersection_delays_df = pd.DataFrame(
            {'node_id': atms_strategy_dict["intersection_node_id"], 'delay': atms_strategy_dict["intersection_delays_in_minutes"]})

        if max(intersection_delays_df["delay"]) > MAX_DELAY:
            sys.exit(
                'The user specified delay is more than the maximum allowed delay for one or more intersections')

        cols = ['ea_delay', 'am_delay', 'md_delay', 'pm_delay', 'ev_delay']

        for c in cols:
            intersection_delays_df[c] = intersection_delays_df['delay']

        intersection_delays_df = intersection_delays_df.drop(columns=["delay"])

    # get two approaching links for each atms intersection
    corridor_links_df = get_corridor_links(links_df, atms_strategy_dict)
    intersection_links_df = corridor_links_df[corridor_links_df.intersection > 0]

    # get link attributes for approaching links for all intersections
    intersection_links_df = get_link_attributes(intersection_links_df, highway_load_data_dict)

    # prepare data
    data_df = intersection_links_df.copy()

    data_df["am_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_am_flow"],
        data_df["ba_am_flow"]
    )

    data_df["pm_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_pm_flow"],
        data_df["ba_pm_flow"]
    )

    data_df["md_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_md_flow"],
        data_df["ba_md_flow"]
    )

    data_df["ea_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_ea_flow"],
        data_df["ba_ea_flow"]
    )

    data_df["ev_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_ev_flow"],
        data_df["ba_ev_flow"]
    )

    data_df["am_auto_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_am_auto_flow"],
        data_df["ba_am_auto_flow"]
    )

    data_df["pm_auto_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_pm_auto_flow"],
        data_df["ba_pm_auto_flow"]
    )

    data_df["md_auto_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_md_auto_flow"],
        data_df["ba_md_auto_flow"]
    )

    data_df["ea_auto_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_ea_auto_flow"],
        data_df["ba_ea_auto_flow"]
    )

    data_df["ev_auto_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_ev_auto_flow"],
        data_df["ba_ev_auto_flow"]
    )

    data_df["am_truck_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_am_truck_flow"],
        data_df["ba_am_truck_flow"]
    )

    data_df["pm_truck_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_pm_truck_flow"],
        data_df["ba_pm_truck_flow"]
    )

    data_df["md_truck_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_md_truck_flow"],
        data_df["ba_md_truck_flow"]
    )

    data_df["ea_truck_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_ea_truck_flow"],
        data_df["ba_ea_truck_flow"]
    )

    data_df["ev_truck_flow"] = np.where(
        data_df["approach"] == "AB",
        data_df["ab_ev_truck_flow"],
        data_df["ba_ev_truck_flow"]
    )

    data_df["from_node"] = np.where(
        data_df["approach"] == "AB",
        data_df["from_node"],
        data_df["to_node"]
    )

    data_df["to_node"] = np.where(
        data_df["approach"] == "AB",
        data_df["to_node"],
        data_df["from_node"]
    )

    # join intersection link flows and delays into one dataframe
    flows_df = data_df.copy()

    if USER_SPECIFIED_DELAY == 0:
        flows_df = pd.merge(
            flows_df,
            intersection_delays_df,
            how="left",
            left_on=["from_node", "to_node"],
            right_on=["From", "To"]
        )
    else:
        flows_df = pd.merge(
            flows_df,
            intersection_delays_df,
            how="left",
            left_on=["to_node"],
            right_on=["node_id"]
        )

    # calculate total delay savings for passenger cars and trucks
    delay_saving_df = flows_df.copy()

    delay_saving_df["ea_auto_delay_savings"] = (delay_saving_df["ea_auto_flow"]
                                                * (delay_saving_df["ea_delay"] / 60.0)
                                                * (total_percent_delay_reduction / 100.0))

    delay_saving_df["am_auto_delay_savings"] = (delay_saving_df["am_auto_flow"]
                                                * (delay_saving_df["am_delay"] / 60.0)
                                                * (total_percent_delay_reduction / 100.0))

    delay_saving_df["md_auto_delay_savings"] = (delay_saving_df["md_auto_flow"]
                                                * (delay_saving_df["md_delay"] / 60.0)
                                                * (total_percent_delay_reduction / 100.0))
    delay_saving_df["pm_auto_delay_savings"] = (delay_saving_df["pm_auto_flow"]
                                                * (delay_saving_df["pm_delay"] / 60.0)
                                                * (total_percent_delay_reduction / 100.0))

    delay_saving_df["ev_auto_delay_savings"] = (delay_saving_df["ev_auto_flow"]
                                                * (delay_saving_df["ev_delay"] / 60.0)
                                                * (total_percent_delay_reduction / 100.0))

    delay_saving_df["ea_truck_delay_savings"] = (delay_saving_df["ea_truck_flow"]
                                                 * (delay_saving_df["ea_delay"] / 60.0)
                                                 * (total_percent_delay_reduction / 100.0))

    delay_saving_df["am_truck_delay_savings"] = (delay_saving_df["am_truck_flow"]
                                                 * (delay_saving_df["am_delay"] / 60.0)
                                                 * (total_percent_delay_reduction / 100.0))

    delay_saving_df["md_truck_delay_savings"] = (delay_saving_df["md_truck_flow"]
                                                 * (delay_saving_df["md_delay"] / 60.0)
                                                 * (total_percent_delay_reduction / 100.0))

    delay_saving_df["pm_truck_delay_savings"] = (delay_saving_df["pm_truck_flow"]
                                                 * (delay_saving_df["pm_delay"] / 60.0)
                                                 * (total_percent_delay_reduction / 100.0))

    delay_saving_df["ev_truck_delay_savings"] = (delay_saving_df["ev_truck_flow"]
                                                 * (delay_saving_df["ev_delay"] / 60.0)
                                                 * (total_percent_delay_reduction / 100.0))

    delay_saving_df["daily_auto_delay_savings"] = (delay_saving_df["ea_auto_delay_savings"]
                                                   + delay_saving_df["am_auto_delay_savings"]
                                                   + delay_saving_df["md_auto_delay_savings"]
                                                   + delay_saving_df["pm_auto_delay_savings"]
                                                   + delay_saving_df["ev_auto_delay_savings"])

    delay_saving_df["daily_truck_delay_savings"] = (delay_saving_df["ea_truck_delay_savings"]
                                                    + delay_saving_df["am_truck_delay_savings"]
                                                    + delay_saving_df["md_truck_delay_savings"]
                                                    + delay_saving_df["pm_truck_delay_savings"]
                                                    + delay_saving_df["ev_truck_delay_savings"])

    auto_daily_delay_savings = delay_saving_df["daily_auto_delay_savings"].sum()
    truck_daily_delay_savings = delay_saving_df["daily_truck_delay_savings"].sum()

    delay_savings_df = pd.DataFrame.from_dict(
        {
            'Variable': ["Total delays savings for passenger cars (hours)", "Total delay savings for trucks (hours)"],
            'Value': [auto_daily_delay_savings, truck_daily_delay_savings]
        }
    )

    # get emission factors
    emission_factors_df = get_emission_factors(emission_df, scen_year)

    # calculate ghg reductions
    ghg_reduction_df = calculate_ghg_reduction(
        auto_daily_delay_savings,
        truck_daily_delay_savings,
        emission_factors_df
    )

    results_df = pd.concat([delay_savings_df, ghg_reduction_df])

    # writing results
    results_dict = {"GHG_Reduction": results_df, "Emission_Factors": emission_factors_df}
    write_results(results_dict, output_results_filename, output_dir)
