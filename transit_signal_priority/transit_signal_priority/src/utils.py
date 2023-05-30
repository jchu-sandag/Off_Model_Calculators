import os
import pandas as pd
import numpy as np
import geopandas as gpd
import osmnx as ox
import networkx as nx
from shapely.geometry import Point, LineString


def create_nodes_gdf(edges_gdf):
    """
    create nodes geodataframe from the links geodataframe using from/to nodes of links

    Args:
        edges_gdf : dataframe of links

    Returns: dataframe of nodes
    """
    edges_dict = edges_gdf.to_dict('records')

    nodes_list = []
    for row in edges_dict:
        from_node = row["from_node"]
        to_node = row["to_node"]
        from_x, from_y = list(row["geometry"].coords)[0]
        to_x, to_y = list(row["geometry"].coords)[-1]

        from_node_dict = {"id": from_node, "x": from_x, "y": from_y}
        to_node_dict = {"id": to_node, "x": to_x, "y": to_y}

        nodes_list.append(from_node_dict)
        nodes_list.append(to_node_dict)

    nodes_df = pd.DataFrame.from_dict(nodes_list, orient='columns')
    nodes_df.drop_duplicates(subset=["id"], inplace=True)

    nodes_geometries = nodes_df.apply(lambda row: Point(row["x"], row["y"]), axis=1)
    nodes_gdf = gpd.GeoDataFrame(nodes_df, geometry=nodes_geometries)
    nodes_gdf.crs = edges_gdf.crs

    return nodes_gdf


def ox_graph(nodes_df, links_df):
    """
    create an osmnx-flavored network graph

    Args:
        nodes_df : GeoDataFrame of nodes
        link_df : GeoDataFrame of links

    Returns: a networkx multidigraph
    """

    graph_nodes = nodes_df.copy()
    graph_nodes.gdf_name = "network_nodes"

    graph_links = links_df.copy()

    # have to change this over into u,v b/c this is what osmnx is expecting
    graph_links["u"] = graph_links["from_node"]
    graph_links["v"] = graph_links["to_node"]
    graph_links["id"] = graph_links["hwycov_id"]
    graph_links["key"] = graph_links["hwycov_id"]

    graph_links.set_index(['u', 'v', 'key'], inplace=True)

    G = ox.graph_from_gdfs(graph_nodes, graph_links)

    return G


def get_primary_corridor_links(links_df, nodes_df, tsp_intersections_dict):
    """
    return the link ids for the links along the mainline corridor from the
    first tsp intersection to last tsp intersection.

    args:
        links_df: network links dataframe
        nodes_df: network nodes dataframe
        tsp_strategy_dict: dictionary of node id and mainline corridor name for tsp intersections

    returns:
        list of primary link ids

    """
    first_node_id = tsp_intersections_dict["intersection_node_id"][0]
    first_corridor_name = tsp_intersections_dict["mainline_corridor_name"][0]

    last_node_id = tsp_intersections_dict["intersection_node_id"][-1]
    last_corridor_name = tsp_intersections_dict["mainline_corridor_name"][-1]

    mainline_corridor_names = tsp_intersections_dict["mainline_corridor_name"]

    candidate_links_df = links_df.copy()
    candidate_links_df = candidate_links_df[candidate_links_df["link_name"].isin(
        mainline_corridor_names)]

    G = ox_graph(nodes_df, candidate_links_df)

    sp_route = nx.shortest_path(G, first_node_id, last_node_id)

    sp_links = links_df[links_df["from_node"].isin(
        sp_route) & links_df["to_node"].isin(sp_route)]
    sp_link_ids = sp_links["hwycov_id"].tolist()

    first_intersection_links = candidate_links_df[
        (candidate_links_df["from_node"] == first_node_id) | (
            candidate_links_df["to_node"] == first_node_id)
    ]

    last_intersection_links = candidate_links_df[
        (candidate_links_df["from_node"] == last_node_id) | (
            candidate_links_df["to_node"] == last_node_id)
    ]

    other_first_intersection_links = first_intersection_links[~first_intersection_links["hwycov_id"].isin(
        sp_link_ids)]
    other_last_intersection_links = last_intersection_links[~last_intersection_links["hwycov_id"].isin(
        sp_link_ids)]

    add_first_link_ids = other_first_intersection_links["hwycov_id"].tolist()
    add_last_link_ids = other_last_intersection_links["hwycov_id"].tolist()

    selected_link_ids = add_first_link_ids + sp_link_ids + add_last_link_ids

    return selected_link_ids


def get_secondary_links(links_df, tsp_intersections_dict):
    """
    return the link ids for the secondary links at the tsp intersections

    args:
        links_df: network links dataframe
        tsp_strategy_dict: dictionary of node id and mainline corridor name for tsp intersections

    returns:
        list of secondary link ids

    """
    selected_link_ids = []

    for node_id, corridor_name in zip(tsp_intersections_dict["intersection_node_id"], tsp_intersections_dict["mainline_corridor_name"]):
        intersection_links_df = links_df[
            (links_df["from_node"] == node_id) | (links_df["to_node"] == node_id)
        ]

        secondary_links_df = intersection_links_df[
            intersection_links_df["link_name"] != corridor_name
        ]

        selected_link_ids = selected_link_ids + secondary_links_df["hwycov_id"].tolist()

    return selected_link_ids


def get_corridor_links(links_df, tsp_strategy_dict):
    """
    get the corridor links_df. this includes :
    primary and secondary links for each tsp intersection, and
    intermediate links (b/w the tsp intersections) along the mainline corridor

    args:
        links_df: network links dataframe
        tsp_strategy_dict: dictionary of node id and mainline corridor name for tsp intersections

    returns:
        dataframe of identified corridor links

    """
    nodes_df = create_nodes_gdf(links_df)

    primary_corridor_link_ids = get_primary_corridor_links(links_df, nodes_df, tsp_strategy_dict)

    secondary_link_ids = get_secondary_links(links_df, tsp_strategy_dict)

    corridor_link_ids = primary_corridor_link_ids + secondary_link_ids

    corridor_links_df = links_df[links_df["hwycov_id"].isin(corridor_link_ids)]
    corridor_links_df = corridor_links_df[["hwycov_id", "link_name", "from_node", "to_node"]]
    corridor_links_df["corridor"] = np.where(corridor_links_df["hwycov_id"].isin(
        primary_corridor_link_ids), "Primary", "Secondary")

    corridor_links_df["intersection"] = np.NAN
    for i, node in enumerate(tsp_strategy_dict["intersection_node_id"]):
        corridor_links_df["intersection"] = np.where(((corridor_links_df["from_node"] == node) | (
            corridor_links_df["to_node"] == node)), i+1, corridor_links_df["intersection"])

    corridor_links_df["approach"] = None
    corridor_links_df["approach"] = np.where(corridor_links_df["to_node"].isin(
        tsp_strategy_dict["intersection_node_id"]), "AB", corridor_links_df["approach"])
    corridor_links_df["approach"] = np.where(corridor_links_df["from_node"].isin(
        tsp_strategy_dict["intersection_node_id"]), "BA", corridor_links_df["approach"])

    corridor_links_df.sort_values(by=["intersection", "corridor", "approach"], inplace=True)
    corridor_links_df.reset_index(inplace=True, drop=True)

    return corridor_links_df


def get_average_transit_headway(primary_links_df, transit_routes_df, transit_links_df, cycle_length, period):
    """
    get all the period specific transit routes passing through a link and
    calculate the average headway based on total arriving transit vehicles in an hour and
    then average over the all the intersections to get the average intersection headway

    Args:
        primary_links_df : dataframe of primary links for TSP intersections
        transit_routes_df: dataframe of transit routes
        transit_links_df : dataframe of links for transit routes
        period: time period (AM or PM)
        cycle_length: average cycle length

    Returns: average headway for the intersection
    """
    headway_field = period + "_Headway"

    work_df = pd.merge(primary_links_df, transit_links_df,
                       left_on="hwycov_id", right_on="Link_ID", how="left")

    work_df = pd.merge(work_df, transit_routes_df, on="Route_ID", how="left")

    work_df["Direction"] = (
        (work_df["Config"] - 1000 * (work_df["Config"]/1000).astype(int))/100).astype(int)

    work_df["Transit_Veh_Per_Hour"] = 60.0/work_df[headway_field]

    veh_by_dir_df = work_df.groupby(["intersection", "Link_ID", "Direction"], as_index=False).agg({
        "Transit_Veh_Per_Hour": 'sum'})

    veh_by_link_df = veh_by_dir_df.groupby(["intersection", "Link_ID"], as_index=False).agg({
        "Transit_Veh_Per_Hour": 'mean'})

    veh_by_intersection_df = veh_by_link_df.groupby(
        ["intersection"], as_index=False).agg({"Transit_Veh_Per_Hour": 'sum'})

    veh_by_intersection_df["average_headway"] = 60.0/veh_by_intersection_df["Transit_Veh_Per_Hour"]

    average_intersection_headway = round(veh_by_intersection_df["average_headway"].mean(), 2)

    if average_intersection_headway < (cycle_length/60):
        average_intersection_headway = round(cycle_length/60, 2)

    return average_intersection_headway


def get_average_transit_trip_length(primary_links_df, transit_links_df, transit_flow_df, transit_onoff_df, period):
    """
    get all the transit routes touching the TSP intersections primary links and
    compute the average transit trip length

    Args:
        primary_links_df : dataframe of primary links for TSP intersections
        transit_links_df: dataframe of links for transit routes
        transit_flow_df: dataframe of transit flow (stop to stop) by route
        transit_onoff_df : dataframe of transit boardings by route
        period: time period

    Returns: average transit trip length
    """
    primary_links_ids = primary_links_df["hwycov_id"].to_list()

    routes_df = transit_links_df[transit_links_df["Link_ID"].isin(primary_links_ids)]

    route_list = routes_df["Route_ID"].to_list()
    route_list = list(set(route_list))  # to get the unique routes

    flow_df = transit_flow_df[transit_flow_df["ROUTE"].isin(route_list)]
    flow_df = flow_df[flow_df["TOD"] == period]
    flow_df["PassengerMiles"] = flow_df["TRANSITFLOW"] * (flow_df["TOMP"] - flow_df["FROMMP"])

    total_passenger_miles = flow_df["PassengerMiles"].sum()

    boardings_df = transit_onoff_df[transit_onoff_df["ROUTE"].isin(route_list)]
    boardings_df = boardings_df[boardings_df["TOD"] == period]

    total_boardings = boardings_df["BOARDINGS"].sum()

    average_trip_length = round(total_passenger_miles/total_boardings, 2)

    return average_trip_length


def get_average_transit_flow(primary_links_df, transit_aggflow_df, period):
    """
    get all the transit routes touching the TSP intersections primary links and
    compute the average transit flow

    Args:
        primary_links_df : dataframe of primary links for TSP intersections
        transit_aggflow_df: dataframe of transit flow by link
        period: time period

    Returns: average transit flow
    """
    flow_df = transit_aggflow_df[transit_aggflow_df["TOD"] == period]

    flow_df = flow_df.groupby(["LINK_ID"], as_index=False).agg(
        {"AB_TransitFlow": 'sum', "BA_TransitFlow": 'sum'})

    work_df = pd.merge(primary_links_df, flow_df, left_on="hwycov_id",
                       right_on="LINK_ID", how="left")
    work_df["Transit_Flow"] = np.where(
        work_df["approach"] == "AB", work_df["AB_TransitFlow"], work_df["BA_TransitFlow"])

    out_df = work_df.groupby(["intersection"], as_index=False).agg({"Transit_Flow": 'sum'})

    average_transit_flow = round(out_df["Transit_Flow"].mean(), 2)

    return average_transit_flow


def get_corridor_transit_stops(corridor_links_df, transit_routes_df, transit_stops_df):
    """
    get the average (by direction) transit stops across all transit routes within the corridor

    Args:
        corridor_links_df : dataframe of corridor all (primary and secondary) links
        transit_routes_df: dataframe of transit routes
        transit_stops_df: dataframe of transit stops

    Returns: ab_avg_stops and ba_avg_stops
    """
    work_df = corridor_links_df[corridor_links_df["corridor"] == "Primary"]
    work_df = pd.merge(work_df, transit_stops_df, left_on="hwycov_id",
                       right_on="Link_ID", how="left")
    work_df = work_df[~work_df["Stop_ID"].isna()]
    work_df = pd.merge(work_df, transit_routes_df, how="left", on="Route_ID")
    work_df["Direction"] = (
        (work_df["Config"] - 1000 * (work_df["Config"]/1000).astype(int))/100).astype(int)

    ab_stops_df = work_df[work_df["Direction"] == 1]
    ab_stops = ab_stops_df.shape[0]/len(pd.unique(ab_stops_df['Route_ID']))

    ba_stops_df = work_df[work_df["Direction"] == 2]
    ba_stops = ba_stops_df.shape[0]/len(pd.unique(ba_stops_df['Route_ID']))

    return ab_stops, ba_stops


def get_link_attributes(corridor_links_df, highway_load_data_dict, transit_aggflow_df):
    """
    get the link attributes for the corridor links

    Args:
        corridor_links_df : dataframe of corridor all (primary and secondary) links
        highway_load_data_dict: dictionary for highway load dataframe for different time periods
        transit_aggflow_df: dataframe of transit flow by link

    Returns: dataframe with link attributes appended to input dataframe
    """
    out_df = corridor_links_df.copy()

    for period in highway_load_data_dict.keys():
        highway_load_df = highway_load_data_dict[period]

        light_truck_ab_vol_cols = [x for x in highway_load_df.columns if "AB_Flow_lhd" in x]
        medium_truck_ab_vol_cols = [x for x in highway_load_df.columns if "AB_Flow_mhd" in x]
        heavy_truck_ab_vol_cols = [x for x in highway_load_df.columns if "AB_Flow_hhd" in x]
        truck_ab_vol_cols = light_truck_ab_vol_cols + medium_truck_ab_vol_cols + heavy_truck_ab_vol_cols

        light_truck_ba_vol_cols = [x for x in highway_load_df.columns if "BA_Flow_lhd" in x]
        medium_truck_ba_vol_cols = [x for x in highway_load_df.columns if "BA_Flow_mhd" in x]
        heavy_truck_ba_vol_cols = [x for x in highway_load_df.columns if "BA_Flow_hhd" in x]
        truck_ba_vol_cols = light_truck_ba_vol_cols + medium_truck_ba_vol_cols + heavy_truck_ba_vol_cols

        highway_load_df["AB_Truck_Flow"] = highway_load_df.loc[:, truck_ab_vol_cols].sum(axis=1)

        highway_load_df["BA_Truck_Flow"] = highway_load_df.loc[:, truck_ba_vol_cols].sum(axis=1)

        highway_load_df["AB_Auto_Flow"] = (highway_load_df["AB_Flow"]
                                           - highway_load_df["AB_Truck_Flow"])
        highway_load_df["BA_Auto_Flow"] = (highway_load_df["BA_Flow"]
                                           - highway_load_df["BA_Truck_Flow"])

        highway_load_df["AB_Capacity"] = np.where(
            highway_load_df["AB_VOC"] > 0, highway_load_df["AB_Flow"]/highway_load_df["AB_VOC"], 0)
        highway_load_df["BA_Capacity"] = np.where(
            highway_load_df["BA_VOC"] > 0, highway_load_df["BA_Flow"]/highway_load_df["BA_VOC"], 0)

        highway_load_df = highway_load_df[["ID1", "AB_Time", "BA_Time", "AB_Flow", "BA_Flow",
                                           "AB_Auto_Flow", "BA_Auto_Flow", "AB_Truck_Flow", "BA_Truck_Flow", "AB_Capacity", "BA_Capacity"]]

        transit_flow_df = transit_aggflow_df[transit_aggflow_df["TOD"] == period]
        transit_flow_df = transit_flow_df.groupby(["LINK_ID"], as_index=False).agg(
            {"AB_TransitFlow": 'sum', "BA_TransitFlow": 'sum'})

        out_df = pd.merge(out_df, highway_load_df, how="left",
                          left_on="hwycov_id", right_on="ID1")

        out_df = pd.merge(out_df, transit_flow_df, how="left",
                          left_on="hwycov_id", right_on="LINK_ID")

        out_df.rename(
            columns={
                "AB_Time": "ab_" + period.lower() + "_time",
                "BA_Time": "ba_" + period.lower() + "_time",
                "AB_Flow": "ab_" + period.lower() + "_flow",
                "BA_Flow": "ba_" + period.lower() + "_flow",
                "AB_Auto_Flow": "ab_" + period.lower() + "_auto_flow",
                "BA_Auto_Flow": "ba_" + period.lower() + "_auto_flow",
                "AB_Truck_Flow": "ab_" + period.lower() + "_truck_flow",
                "BA_Truck_Flow": "ba_" + period.lower() + "_truck_flow",
                "AB_Capacity": "ab_" + period.lower() + "_capacity",
                "BA_Capacity": "ba_" + period.lower() + "_capacity",
                "AB_TransitFlow": "ab_" + period.lower() + "_transit_flow",
                "BA_TransitFlow": "ba_" + period.lower() + "_transit_flow"
            },
            inplace=True
        )

        out_df.drop(columns=["ID1", "LINK_ID"], inplace=True)

    return out_df


def get_emission_factors(input_emission_data, scen_year):
    """
    get a dataframe of emission factors for the scen year

    args:
        input_emission_data: dataframe with all emission rates
        scen_year: scen_year
    returns:
        dataframe with emission rates for the scen year

    """
    emission_df = input_emission_data.copy()
    emission_df = emission_df[emission_df["Year"] == scen_year]
    emission_df.reset_index(inplace=True, drop=True)

    return emission_df


def calculate_ghg_reduction(delay_savings_df, vmt_savings, emission_factors):
    """
    calculates total ghg reductions because of delay savings and vmt savings

    args:
        delay_savings_df: dataframe with delay saving because of TSP and by vehicle type
        vmt_savings: vmt reduction due to marginal mode shift to transit
        emission_factors: emission factor dataframe
    returns:
        total ghg reductions

    """
    GRAMS_TO_SHORT_TONS = 0.0000011

    auto_runex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                  "Passenger Car"].reset_index().at[0, "CO2 RunEx Emission Factor (gr/mile)"]

    auto_idlex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                  "Passenger Car"].reset_index().at[0, "CO2 IdlEx Emission Factor (gr/hour)"]

    truck_idlex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                   "HDT - All Weight Types"].reset_index().at[0, "CO2 IdlEx Emission Factor (gr/hour)"]

    bus_idlex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                 "Bus - All Fuel Types"].reset_index().at[0, "CO2 IdlEx Emission Factor (gr/hour)"]

    ghg_reductions_vmt_savings = (vmt_savings
                                  * auto_runex_emission_factor
                                  * GRAMS_TO_SHORT_TONS)

    idling_factors = {
        'type': ["auto", "truck", "bus"],
        'idlex_factor': [auto_idlex_emission_factor, truck_idlex_emission_factor, bus_idlex_emission_factor]
    }

    idling_emission_factor_df = pd.DataFrame(idling_factors)

    work_df = delay_savings_df.copy()
    work_df = pd.merge(work_df, idling_emission_factor_df, how='left', on='type')

    work_df["ghg_reductions"] = (work_df["total_delay_savings"]
                                 * work_df["idlex_factor"]
                                 * GRAMS_TO_SHORT_TONS)

    ghg_reductions_delay_savings = work_df["ghg_reductions"].sum()

    total_ghg_reductions = ghg_reductions_vmt_savings + ghg_reductions_delay_savings

    out_data = {
        "Variable": ["GHG reduction due to delay savings (short tons)",
                     "GHG reduction due to vmt savings (short tons)",
                     "Total GHG reduction (short tons)"],
        "Value": [ghg_reductions_delay_savings,
                  ghg_reductions_vmt_savings,
                  total_ghg_reductions]
    }

    out_df = pd.DataFrame(out_data)

    return out_df


def write_results(results_dict, out_file_name, out_dir):
    with pd.ExcelWriter(os.path.join(out_dir, out_file_name)) as writer:
        for key, value in results_dict.items():
            value.to_excel(writer, sheet_name=key, index=False)
