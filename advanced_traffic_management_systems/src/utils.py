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


def get_primary_corridor_links(links_df, nodes_df, atms_intersections_dict):
    """
    return the link ids for the links along the mainline corridor from the
    first tsp intersection to last tsp intersection.

    Args:
        links_df: network links dataframe
        nodes_df: network nodes dataframe
        atms_intersections_dict: dictionary of node id and mainline corridor name for atms intersections

    Returns:
        list of primary link ids

    """
    first_node_id = atms_intersections_dict["intersection_node_id"][0]
    first_corridor_name = atms_intersections_dict["mainline_corridor_name"][0]

    last_node_id = atms_intersections_dict["intersection_node_id"][-1]
    last_corridor_name = atms_intersections_dict["mainline_corridor_name"][-1]

    mainline_corridor_names = atms_intersections_dict["mainline_corridor_name"]

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


def get_corridor_links(links_df, atms_strategy_dict):
    """
    get the primary corridor links for all atms intersections

    Args:
        links_df: network links dataframe
        atms_strategy_dict: dictionary of node id and mainline corridor name for atms intersections

    Returns:
        dataframe of identified corridor links

    """
    nodes_df = create_nodes_gdf(links_df)

    primary_corridor_link_ids = get_primary_corridor_links(links_df, nodes_df, atms_strategy_dict)

    corridor_links_df = links_df[links_df["hwycov_id"].isin(primary_corridor_link_ids)]
    corridor_links_df = corridor_links_df[["hwycov_id", "link_name", "from_node", "to_node"]]
    corridor_links_df["corridor"] = np.where(corridor_links_df["hwycov_id"].isin(
        primary_corridor_link_ids), "Primary", "Secondary")

    corridor_links_df["intersection"] = np.NAN
    for i, node in enumerate(atms_strategy_dict["intersection_node_id"]):
        corridor_links_df["intersection"] = np.where(((corridor_links_df["from_node"] == node) | (
            corridor_links_df["to_node"] == node)), i+1, corridor_links_df["intersection"])

    corridor_links_df["approach"] = None
    corridor_links_df["approach"] = np.where(corridor_links_df["to_node"].isin(
        atms_strategy_dict["intersection_node_id"]), "AB", corridor_links_df["approach"])
    corridor_links_df["approach"] = np.where(corridor_links_df["from_node"].isin(
        atms_strategy_dict["intersection_node_id"]), "BA", corridor_links_df["approach"])

    corridor_links_df.sort_values(by=["intersection", "corridor", "approach"], inplace=True)
    corridor_links_df.reset_index(inplace=True, drop=True)

    return corridor_links_df


def get_link_attributes(corridor_links_df, highway_load_data_dict):
    """
    get the link attributes for the corridor links

    Args:
        corridor_links_df : dataframe of corridor links
        highway_load_data_dict: dictionary for highway load dataframe for different time periods

    Returns:
        dataframe with link attributes appended to input dataframe

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

        highway_load_df = highway_load_df[["ID1", "AB_Flow", "BA_Flow",
                                           "AB_Auto_Flow", "BA_Auto_Flow", "AB_Truck_Flow", "BA_Truck_Flow"]]

        out_df = pd.merge(
            out_df,
            highway_load_df,
            how="left",
            left_on="hwycov_id",
            right_on="ID1"
        )

        out_df.rename(
            columns={
                "AB_Flow": "ab_" + period.lower() + "_flow",
                "BA_Flow": "ba_" + period.lower() + "_flow",
                "AB_Auto_Flow": "ab_" + period.lower() + "_auto_flow",
                "BA_Auto_Flow": "ba_" + period.lower() + "_auto_flow",
                "AB_Truck_Flow": "ab_" + period.lower() + "_truck_flow",
                "BA_Truck_Flow": "ba_" + period.lower() + "_truck_flow",
            },
            inplace=True
        )

        out_df.drop(columns=["ID1"], inplace=True)

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


def calculate_ghg_reduction(daily_auto_delay_savings, daily_truck_delay_savings, emission_factors):
    """
    calculates total ghg reductions because of delay savings

    Args:
        daily_auto_delay_savings: total delays savings (hours) for passenger cars
        daily_truck_delay_savings: total delay savings (hours) for trucks
        emission_factors: emission factor dataframe

    Returns:
        total ghg reductions

    """
    GRAMS_TO_SHORT_TONS = 0.0000011

    auto_idlex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                  "Passenger Car"].reset_index().at[0, "CO2 IdlEx Emission Factor (gr/hour)"]

    truck_idlex_emission_factor = emission_factors[emission_factors["Vehicle Type"] ==
                                                   "HDT - All Weight Types"].reset_index().at[0, "CO2 IdlEx Emission Factor (gr/hour)"]

    passenger_car_daily_ghg_reduction = (daily_auto_delay_savings
                                         * auto_idlex_emission_factor
                                         * GRAMS_TO_SHORT_TONS)

    truck_daily_ghg_reduction = (daily_truck_delay_savings
                                 * truck_idlex_emission_factor
                                 * GRAMS_TO_SHORT_TONS)

    total_daily_ghg_reduction = passenger_car_daily_ghg_reduction + truck_daily_ghg_reduction

    out_data = {
        "Variable": ["GHG reduction due to delay savings for passenger cars (short tons)",
                     "GHG reduction due to delay savings for trucks (short tons)",
                     "Total GHG reduction (short tons)"],
        "Value": [passenger_car_daily_ghg_reduction,
                  truck_daily_ghg_reduction,
                  total_daily_ghg_reduction]
    }

    out_df = pd.DataFrame(out_data)

    return out_df


def write_results(results_dict, out_file_name, out_dir):
    with pd.ExcelWriter(os.path.join(out_dir, out_file_name)) as writer:
        for key, value in results_dict.items():
            value.to_excel(writer, sheet_name=key, index=False)
