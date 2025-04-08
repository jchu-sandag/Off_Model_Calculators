# Python 3.11.7

# Import libraries
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

import numpy as np
import openmatrix as omx
import yaml
import math
import sys
import os
import time
from datetime import datetime
import psutil

# Function: record start
def print_start_time():
    
    print("============================================")
    # Record start time
    start_time = time.time()
    
    # Print current time
    current_time = datetime.now()
    formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    print("Processing Start Time: %s" % formatted_time)
    
    # Print memory space in use
    used_memory = psutil.Process().memory_info().rss / 1024 ** 2 # Convert to MB
    total_memory = psutil.virtual_memory().total / 1024 ** 2 # Convert to MB
    print("Memory space in use {:,.0f} / {:,.0f} MB".format(used_memory, total_memory))
    print("============================================")
    return start_time

# Function: record end time
def print_end_time(start_time):
    
    print("============================================")
    # Record end time
    end_time = time.time()
    
    # Print the end time
    current_time = datetime.now()
    formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    print("Processing End Time: %s" % formatted_time)

    # Calculate elapsed time
    elapsed_time = end_time - start_time
    elapsed_minutes = int(elapsed_time // 60)
    elapsed_seconds = int(elapsed_time % 60)
    print("Elapsed time: %s minutes and %s seconds" % (elapsed_minutes, elapsed_seconds))
    
    # Print memory space in use
    used_memory = psutil.Process().memory_info().rss / 1024 ** 2 # Convert to MB
    total_memory = psutil.virtual_memory().total / 1024 ** 2 # Convert to MB
    print("Memory space in use {:,.0f} / {:,.0f} MB".format(used_memory, total_memory))
    print("============================================")
    return end_time

#### READ INPUTS & OUTPUTS FROM CONFIG FILE

# Open config file & read inputs
start_time = print_start_time()

###########################################################################
########   READ INPUT/OUTPUT FILE NAMES FROM CONFIG FILE   ################
###########################################################################

if __name__ == "__main__":
    args = sys.argv
    config_filename = args[1]
    # print(config_filename)

    if not os.path.exists(config_filename):
        msg = "Configuration file doesn't exist at: {}".format(config_filename)
        raise ValueError(msg)

    with open(config_filename, "r") as yml_file:
        config = yaml.safe_load(yml_file)

    # Read config: inputs

    # ABM3 Data
    mgra_base_input_file = config['inputs']['mgra_base_input_file']
    mgra_scen_input_file = config['inputs']['mgra_scen_input_file']
    skim_base_file = config['inputs']['skim_base_file']
    skim_scen_file = config['inputs']['skim_scen_file']
    individual_tours_file = config['inputs']['individual_tours_output_file']

    # External Data
    vanpool_od_file = config['inputs']['vanpool_od_file'] 
    employment_forecast_scag_file = config['inputs']['employment_forecast_scag_file'] 
    emission_factors_file = config['inputs']['emission_factors_file']
    zipcode_coordinates_file = config['inputs']['zipcode_coordinates_file']
    external_gateways_file = config['inputs']['external_gateways_file']
    geography_xwalk_file = config['inputs']['geography_xwalk_file'] 
    msa_names_file = config['inputs']['msa_names_file'] 

    # Read config: outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # Read config: parameters
    base_year = config['parameters']['base_year']
    scen_year = config['parameters']['scen_year']
    c_ivt = config['parameters']['c_ivt'] 
    avg_vanpool_occupancy = config['parameters']['avg_vanpool_occupancy'] 
    pct_work_trips_over_50mi = config['parameters']['pct_work_trips_over_50mi'] 
    sov_time_core_name = config['parameters']['sov_am_time_core_name'] 
    hov_time_core_name = config['parameters']['hov_am_time_core_name'] 
    military_base_taz = config['parameters']['military_base_taz'] 
    abm_version = config['parameters']['abm_version']

    ##############################################
    ########   READ INPUT FILES   ################
    ##############################################

    vanpool_od_df = pd.read_csv(vanpool_od_file)
    mgra_scen_input_df = pd.read_csv(mgra_scen_input_file)
    mgra_base_input_df = pd.read_csv(mgra_base_input_file)
    geo_xwalk_df = pd.read_csv(geography_xwalk_file)
    msa_names_df = pd.read_csv(msa_names_file)
    emp_forecast_scag_df = pd.read_csv(employment_forecast_scag_file)
    emission_df = pd.read_excel(emission_factors_file)
    zipcode_coordinates_df = pd.read_csv(zipcode_coordinates_file)
    external_gateways_df = pd.read_csv(external_gateways_file)
    indiv_tours_df = pd.read_csv(individual_tours_file)

    #### VANPOOL DEMAND DUE TO REGIONAL EMPLOYMENT GROWTH

    # calculate distance traveled inside/outside of county

    # drop unused fields from vanpool df
    data_df = vanpool_od_df.drop(columns = ["Employer", "Van_Type", "Industry"])

    # join origin/destination zip code coordinates to vanpool df
    od_type = ["O", "D"]

    for od in od_type:
        data_df = data_df.merge(zipcode_coordinates_df, left_on = od + "_Zip", right_on = "Zip_Code", how = "left")
        data_df.rename(
            columns={"Lat" : od + "_Lat",
                    "Long" : od + "_Long",
                    "Lat_ft" : od + "_Lat_ft",
                    "Long_ft" : od + "_Long_ft"
                    },
            inplace = True
        )
        data_df.drop(columns = ["Zip_Code"], inplace=True)

    # join external gateways
    data_df = data_df.merge(external_gateways_df, left_on = "O_MSA", right_on = "Home_County", how = "left")

    # calculate van miles traveled outside of SD county (external distance) based on distance between external gateway lat/long and outside county (any county outside of SD)
    data_df["OutofCounty_Oneway_Mileage"] = round(np.sqrt(
        (data_df["O_Lat_ft"] - data_df["Gateway_Lat_ft"]) ** 2 +
        (data_df["O_Long_ft"] - data_df["Gateway_Long_ft"]) ** 2
        ) / 5280, 2)

    # for all interal trips, replace out of county oneway mileage (NA) value with Zero
    data_df["OutofCounty_Oneway_Mileage"].fillna(0, inplace=True)

    # subtract (2 X inside county distance traveled) from total daily round trip mileage to estimate outside travel distance
    data_df["InCounty_Round_Trip_Mileage"] = np.maximum(
        data_df["Daily_Round_Trip_Mileage"] - 2 * (data_df["OutofCounty_Oneway_Mileage"]), 0)

    # remove unused fields
    data_df.drop(columns = ["O_Lat_ft", "O_Long_ft", "D_Lat_ft", "D_Long_ft", "TAZ", "Home_County", "Gateway_Lat_ft", "Gateway_Long_ft"], inplace=True)

    data_df.head()

    # estimate In County emp growth rates

    ### Base Year
    # write input land use data into temporary df
    emp_base_by_msa_internal = mgra_base_input_df[["mgra", "emp_mil", "emp_total"]]
    emp_base_by_msa_internal["emp_non_mil"] = emp_base_by_msa_internal["emp_total"] - emp_base_by_msa_internal["emp_mil"]

    # merge taz, msa name, and external county name to temp df
    emp_base_by_msa_internal = emp_base_by_msa_internal.merge(geo_xwalk_df, how="left", on="mgra")
    emp_base_by_msa_internal = emp_base_by_msa_internal.merge(msa_names_df, how="left", on="msa")

    # groupby msa and aggregate number of employees in each group
    emp_base_by_msa_internal = emp_base_by_msa_internal.groupby(["msa_name"]).agg({"emp_mil": 'sum', "emp_non_mil": 'sum'})
    emp_base_by_msa_internal.reset_index(inplace=True)
    emp_base_by_msa_internal.columns = ["MSA_Name", "emp_mil_base", "emp_non_mil_base"]

    emp_base_by_msa_internal.head()

    ### Scenario Year
    # write input land use data into temporary df
    emp_scen_by_msa_internal = mgra_scen_input_df[["mgra", "emp_mil", "emp_total"]]
    emp_scen_by_msa_internal["emp_non_mil"] = emp_scen_by_msa_internal["emp_total"] - emp_scen_by_msa_internal["emp_mil"]

    # merge taz, msa name, and external county name to temp df
    emp_scen_by_msa_internal = emp_scen_by_msa_internal.merge(geo_xwalk_df, how="left", on="mgra")
    emp_scen_by_msa_internal = emp_scen_by_msa_internal.merge(msa_names_df, how="left", on="msa")

    # groupby msa and aggregate number of employees in each group
    emp_scen_by_msa_internal = emp_scen_by_msa_internal.groupby(["msa_name"]).agg({"emp_mil": 'sum', "emp_non_mil": 'sum'})
    emp_scen_by_msa_internal.reset_index(inplace=True)
    emp_scen_by_msa_internal.columns = ["MSA_Name", "emp_mil_scen", "emp_non_mil_scen"]

    # combine base year and scenario year data into one df
    emp_growth_internal = pd.merge(emp_base_by_msa_internal, emp_scen_by_msa_internal, how="left", on="MSA_Name")

    # calculate growth rates
    emp_growth_internal["emp_mil_growth"] = np.where(emp_growth_internal["emp_mil_base"] != 0, round(emp_growth_internal["emp_mil_scen"] / emp_growth_internal["emp_mil_base"] ,2), 1)
    emp_growth_internal["emp_non_mil_growth"] = np.where(emp_growth_internal["emp_non_mil_base"] != 0, round(emp_growth_internal["emp_non_mil_scen"] / emp_growth_internal["emp_non_mil_base"], 2), 1)

    emp_growth_internal

    # estimate Out of County emp growth rates

    ### Base Year
    # read emp popualtion from external msa by year
    emp_base_by_msa_external = emp_forecast_scag_df[emp_forecast_scag_df["Year"] == base_year]
    emp_base_by_msa_external.drop(columns = ["Year"], inplace = True)
    emp_base_by_msa_external.rename(columns = {"Employment" : "emp_total_base"}, inplace = True)

    ### Scenario Year
    # read emp popualtion from external msa by year
    emp_scen_by_msa_external = emp_forecast_scag_df[emp_forecast_scag_df["Year"] == scen_year]
    emp_scen_by_msa_external.drop(columns = ["Year"], inplace = True)
    emp_scen_by_msa_external.rename(columns = {"Employment" : "emp_total_scen"}, inplace = True)

    # combine base year and scenario year data into one df
    emp_growth_external = pd.merge(emp_base_by_msa_external, emp_scen_by_msa_external, how="left", on="County")

    # calculate growth rates
    # IMPORTANT! since military/non-military data is not availabe from the SCAG outputs, growth rate of both groups is assumed same as total emp growth
    emp_growth_external["emp_mil_growth"] = np.where(emp_growth_external["emp_total_base"] != 0, round(emp_growth_external["emp_total_scen"] / emp_growth_external["emp_total_base"] ,2), 1)
    emp_growth_external["emp_non_mil_growth"] = emp_growth_external["emp_mil_growth"]
    emp_growth_external.rename(columns = {"County" : "MSA_Name"}, inplace = True)

    emp_growth_external

    # combine emp growth into one table
    emp_growth = pd.concat([emp_growth_internal[["MSA_Name", "emp_mil_growth", "emp_non_mil_growth"]], emp_growth_external[["MSA_Name", "emp_mil_growth", "emp_non_mil_growth"]]])
    emp_growth.columns = ["MSA_Name", "Military", "Non_Military"]

    # write growth factors to new df
    emp_growth = pd.melt(emp_growth,
                        id_vars = ["MSA_Name"],
                        value_vars = ["Military", "Non_Military"],
                        var_name = "Industry_Type",
                        value_name = "Emp_Growth"
                        )

    emp_growth

    # calculate van numbers in scenario year
    # count total vans by industry type and destination
    data_df_vans_by_dest = data_df.groupby(by = ["Industry_Type", "D_MSA"])["VAN_ID"].count().to_frame()
    data_df_vans_by_dest.reset_index(inplace = True)
    data_df_vans_by_dest.rename(columns = {"VAN_ID": "Num_Vans_By_Dest_Base"}, inplace = True)

    # calculate vanpool in scenario year based on employment growth
    data_df_vans_by_dest = pd.merge(
        data_df_vans_by_dest,
        emp_growth,
        left_on=["D_MSA", "Industry_Type"],
        right_on=["MSA_Name", "Industry_Type"],
        how="left"
        )

    # drop unused fields
    data_df_vans_by_dest.drop(columns = ["MSA_Name"], inplace = True)

    # apply growth percentage
    data_df_vans_by_dest["Num_Vans_By_Dest_Scen"] = data_df_vans_by_dest["Num_Vans_By_Dest_Base"] * data_df_vans_by_dest["Emp_Growth"]

    # round numbers up
    data_df_vans_by_dest["Num_Vans_By_Dest_Scen"] = data_df_vans_by_dest["Num_Vans_By_Dest_Scen"].apply(lambda x: math.ceil(x))

    data_df_vans_by_dest

    # count total vans by industry type, origin and destination
    data_df_vans_by_od = data_df.groupby(["Industry_Type", "O_MSA", "D_MSA"])["VAN_ID"].count().to_frame()
    data_df_vans_by_od.reset_index(inplace = True)
    data_df_vans_by_od.rename(columns = {"VAN_ID": "Num_Vans_Base"}, inplace = True)

    # merge Vanpool OD Base and Vanpool Growth by destination. The objective is to estimate Vanpool OD for scenario year
    data_df_vans_by_od = pd.merge(
        data_df_vans_by_od,
        data_df_vans_by_dest,
        on = ["Industry_Type", "D_MSA"],
        how="left"
        )

    # The calculation apply vanpool growth numbers to the OD numbers and rounds the numbers up to nearest integer. Therefore, the total vanpool from previous table does not exactly match the vanpool od sum total numbers.
    # data_df_vans_by_od["Num_Vans_Scen"] = data_df_vans_by_od["Num_Vans_Base"] * (data_df_vans_by_od["Num_Vans_By_Dest_Scen"] / data_df_vans_by_od["Num_Vans_By_Dest_Base"])
    data_df_vans_by_od["Num_Vans_Scen"] = data_df_vans_by_od["Num_Vans_Base"] * data_df_vans_by_od["Emp_Growth"]
    data_df_vans_by_od["Num_Vans_Scen"] = data_df_vans_by_od["Num_Vans_Scen"].apply(lambda x: math.ceil(x))

    data_df_vans_by_od.drop(columns = ["Num_Vans_By_Dest_Base", "Emp_Growth", "Num_Vans_By_Dest_Scen"], inplace = True)
    data_df_vans_by_od

    # ### VANPOOL DEMAND DUE TO MANAGED LANE INFRASTRUCTURE INVESTMENTS

    # estimate probability of vanpool demand

    # estimate vanpool demand (two trips per day)
    weekday_vanpool_demand = vanpool_od_df["Vehicle_Capacity"].sum() * avg_vanpool_occupancy * 2

    # calculate number of work tours based on the individual tour file
    work_tours_df = indiv_tours_df[(indiv_tours_df["primary_purpose"] == 'work') & (indiv_tours_df["tour_mode"] == 'DRIVEALONE')]
    num_work_trips = work_tours_df.shape[0]

    # pct_work_trips_over_50mi is percent of daily work trips with one-way distance of 50 miles or more
    potential_weekday_vanpool_demand = num_work_trips * (pct_work_trips_over_50mi/100) * 2

    # probability of vanpool trips
    prob_vanpool = round(weekday_vanpool_demand / potential_weekday_vanpool_demand, 4)

    print(
        "Weekday Vanpool Demand (Vanpool Report OD) = %s employee trips with vanpool \n" % f"{int(weekday_vanpool_demand):,}" + 
        "Number of work trips (ABM3) = %s emplyee trips across region \n" % f"{num_work_trips:,}" + 
        "Potential Vanpool Demand = %s potential employee trips with vanpool \n" % f"{int(potential_weekday_vanpool_demand):,}" + 
        "Probability of Vanpooling = %s percent" % (round(prob_vanpool * 100, 2))
        )

    # prepare a travel time df including travel times for SOV in BASE year

    ######################################
    ###########   BASE YEAR   ############
    ######################################

    # open the skim_data
    omx_matrix = omx.open_file(skim_base_file)

    # read dimension of the matrix (ROW X COL) Square
    dim = omx_matrix.shape()[0]

    # create a template df with o/d taz. The objective is to convert the skim matrix into a long format (2 columns). Col 1 is Oigin, Col 2 is Destination. Shape of the template df is dim ** 2
    travel_time_base_df = pd.DataFrame({"orig_taz": np.repeat(1 + np.arange(dim), dim), "dest_taz": np.tile(1 + np.arange(dim), dim)})

    # create a df with the Same shape as the template that includes travel times 
    core_names = [sov_time_core_name]
    out_cols = ["sov_time_base"]

    for ii, core in enumerate(core_names):
            core_mtx = omx_matrix[core]
            col = out_cols[ii]
            values_mat = np.array(core_mtx)
            values_df = pd.DataFrame({col: np.reshape(values_mat, (dim ** 2))})
            travel_time_base_df = pd.concat([travel_time_base_df, values_df], axis = 1)

    travel_time_base_df.head(25)

    # prepare a travel time df including travel times for SOV, HOV2, and HOV3 in SCEN year

    ######################################
    ###########   SCEN YEAR   ############
    ######################################

    # open the skim_data
    omx_matrix = omx.open_file(skim_scen_file)

    # read dimension of the matrix (ROW X COL) Square
    dim = omx_matrix.shape()[0]

    # create a template df with o/d taz. The objective is to convert the skim matrix into a long format (2 columns). Col 1 is Oigin, Col 2 is Destination. Shape of the template df is dim ** 2
    travel_time_scen_df = pd.DataFrame({"orig_taz": np.repeat(1 + np.arange(dim), dim), "dest_taz": np.tile(1 + np.arange(dim), dim)})

    # create a df with the Same shape as the template that includes travel times

    core_names = [sov_time_core_name, hov_time_core_name]
    out_cols = ["sov_time_scen", "hov_time_scen"]

    for ii, core in enumerate(core_names):
            core_mtx = omx_matrix[core]
            col = out_cols[ii]
            values_mat = np.array(core_mtx)
            values_df = pd.DataFrame({col: np.reshape(values_mat, (dim ** 2))})
            travel_time_scen_df = pd.concat([travel_time_scen_df, values_df], axis = 1)

    travel_time_scen_df.head(25)

    # merge base and scenario year travel times into a single df
    travel_time_df = pd.merge(
            travel_time_base_df,
            travel_time_scen_df,
            on = ["orig_taz", "dest_taz"],
            how = "left"
        )

    # write Orig/Dest MSA names to the travel time df
    # the task is done in 2 steps (step 1 is internal MSA, step 2 is external MSA)

    # create df that includes taz and associated msa
    taz_msa_df = geo_xwalk_df[["taz", "msa"]].drop_duplicates(subset = ["taz"])

    ##################################################
    ############   INTERNAL MSA   ####################
    ##################################################

    # add msa origin and destination associated with each taz
    ### ORIGIN MSA (Internal MSA)
    travel_time_df = pd.merge(travel_time_df, taz_msa_df, left_on = "orig_taz", right_on = "taz", how = "left")
    travel_time_df = pd.merge(travel_time_df, msa_names_df, on = "msa", how = "left")
    travel_time_df.rename(columns={"msa_name": "O_MSA"}, inplace = True)
    travel_time_df.drop(["taz", "msa"], axis = 1, inplace = True)

    ### DESTINATION MSA (Internal MSA)
    travel_time_df = pd.merge(travel_time_df, taz_msa_df, left_on = "dest_taz", right_on = "taz", how="left")
    travel_time_df = pd.merge(travel_time_df, msa_names_df, on = "msa", how = "left")
    travel_time_df.rename(columns={"msa_name": "D_MSA"}, inplace = True)
    travel_time_df.drop(columns = ["taz", "msa"], axis = 1, inplace = True)

    travel_time_df.head()

    ##################################################
    ############   EXTERNAL MSA   ####################
    ##################################################

    # for each external MSA, there is a TAZ associated with the inbound/outbound traffic. (example, trips to/from Orange County enter/exit from TAZ 12 (in ABM2+ data))
    # in other words, Origin/Destination MSA is determined based on the TAZ where trip enter/exit SD County
    gateways_df = external_gateways_df[["Home_County", "TAZ"]]
    travel_time_df = pd.merge(travel_time_df, gateways_df, left_on = "orig_taz", right_on = "TAZ", how = "left")

    ### ORIGIN MSA (External MSA)
    # for rows where the Home County field has data, replace the O_MSA with Home County value
    travel_time_df["O_MSA"] = np.where(travel_time_df["Home_County"].isna(),
                                            travel_time_df["O_MSA"],
                                            travel_time_df["Home_County"])
    travel_time_df.rename(columns={"TAZ": "orig_ext_taz"}, inplace = True)
    travel_time_df.drop(["Home_County"], axis = 1, inplace = True)

    ### DESTINATION MSA (External MSA)
    # for rows where the Home County field has data, replace the D_MSA with Home County value
    travel_time_df = pd.merge(travel_time_df, gateways_df, left_on="dest_taz", right_on="TAZ", how="left")
    travel_time_df["D_MSA"] = np.where(travel_time_df["Home_County"].isna(),
                                            travel_time_df["D_MSA"],
                                            travel_time_df["Home_County"])
    travel_time_df.rename(columns = {"TAZ": "dest_ext_taz"}, inplace = True)
    travel_time_df.drop(["Home_County"], axis=1, inplace=True)

    # for some of the OD pairs, the msa value is None
    # the reason is that some of the external TAZ (entry/exit gates around the SD county) does not exist today (they are conceptual gates) and are dropped from the calculations
    travel_time_df = travel_time_df[travel_time_df.O_MSA.notna()]
    travel_time_df = travel_time_df[travel_time_df.D_MSA.notna()]

    # drop rows where both the origin AND destination is entry/exit gates (to or from) San Diego county.
    # this means the trip enter from one gate and exit from another
    # there should be very few of such examples compared to all data
    travel_time_df = travel_time_df[~(travel_time_df.orig_ext_taz.notna() & travel_time_df.dest_ext_taz.notna())]

    travel_time_df.head()

    # calculate average travel time savings by comparing tavel times for SOV & HOV by MSA for ALL tazs

    # calculate travel time saveings of each taz
    travel_time_msa_df = travel_time_df.copy()
    travel_time_msa_df["time_savings_scen"] = travel_time_msa_df["hov_time_scen"] - travel_time_msa_df["sov_time_scen"]

    # group travel time savings based on the orig/dest MSA and calculate an average travel time saving between MSA
    travel_time_msa_df = travel_time_msa_df.groupby(["O_MSA", "D_MSA"]).agg({"sov_time_base": 'mean', "time_savings_scen": 'mean'})
    travel_time_msa_df.reset_index(inplace = True)

    travel_time_msa_df.head()

    # calculate average travel time savings by comparing tavel times for SOV & HOV by MSA for MILITARY tazs

    # make a copy of the travel time df
    travel_time_msa_mil_df = travel_time_df.copy()

    # filter data to all trips with destination in military taz
    travel_time_msa_mil_df = travel_time_msa_mil_df[travel_time_msa_mil_df["dest_taz"].isin(military_base_taz)]

    # calculate travel time savings for trips with military taz as destination
    travel_time_msa_mil_df["time_savings_scen"] = travel_time_msa_mil_df["hov_time_scen"] - travel_time_msa_mil_df["sov_time_scen"]

    # group travel time savings based on the orig/dest MSA and calculate an average travel time saving between MSA
    travel_time_msa_mil_df = travel_time_msa_mil_df.groupby(["O_MSA", "D_MSA"]).agg({"sov_time_base": 'mean', "time_savings_scen": 'mean'})
    travel_time_msa_mil_df.reset_index(inplace=True)

    # add mil label to column names
    travel_time_msa_mil_df.rename(columns={'sov_time_base': 'sov_time_base_mil', 'time_savings_scen': 'time_savings_scen_mil'}, inplace=True)

    travel_time_msa_mil_df.head()

    # merge results for military taz to all taz results
    travel_time_msa_df = pd.merge(travel_time_msa_df, travel_time_msa_mil_df, on=["O_MSA", "D_MSA"], how="left")

    # replace all NA values with 0
    travel_time_msa_df["sov_time_base"] = travel_time_msa_df["sov_time_base"].fillna(0)
    travel_time_msa_df["time_savings_scen"] = travel_time_msa_df["time_savings_scen"].fillna(0)
    travel_time_msa_df["sov_time_base_mil"] = travel_time_msa_df["sov_time_base_mil"].fillna(0)
    travel_time_msa_df["time_savings_scen_mil"] = travel_time_msa_df["time_savings_scen_mil"].fillna(0)
        
    travel_time_msa_df.head()

    # calculate demand elasticity with respect to travel time on MSA by MSA matrix
    # The elasticity of demand for vanpooling with respect to travel time was approximated using the formula for point elasticity derived from a logit model (Train, 1993)
    # Elasticity = (coefficient of in-vehicle time) * average travel time * (1 – probability of vanpooling)
    # travel_time_msa_elas_df = travel_time_msa_df.copy()
    travel_time_msa_df["Elasticity_mil"] = c_ivt * travel_time_msa_df["sov_time_base_mil"] * (1 - prob_vanpool)
    travel_time_msa_df["Elasticity_non_mil"] = c_ivt * travel_time_msa_df["sov_time_base"] * (1 - prob_vanpool)


    # calculate percentage growth for mil and non_mil sector due to induced demand from travel time saving achieved from implementing managed lanes
    travel_time_msa_df["Growth_mil"] = np.where(travel_time_msa_df["sov_time_base_mil"] > 0,
                                        1 + (travel_time_msa_df["Elasticity_mil"] * (travel_time_msa_df["time_savings_scen_mil"] / travel_time_msa_df["sov_time_base_mil"])), 1
                                        )

    travel_time_msa_df["Growth_non_mil"] = np.where(travel_time_msa_df["sov_time_base"] > 0,
                                            1 + travel_time_msa_df["Elasticity_non_mil"] * (travel_time_msa_df["time_savings_scen"] / travel_time_msa_df["sov_time_base"]), 1
                                            )

    travel_time_msa_df.head()

    # calculate number of vans in the Scenario year based on calculated growth from travel time savings
    # write calculate growth factors to the vanpool OD data and write into new df
    vanpool_df = pd.merge(
        data_df_vans_by_od,
        travel_time_msa_df[["O_MSA", "D_MSA", "Growth_mil", "Growth_non_mil"]],
        on = ["O_MSA", "D_MSA"],    
        how = "left"
        )

    # estimate number of additional vans based on the calculated growth factor for military & non-military
    vanpool_df["Num_Vans_Scen_ML"] = np.where(
        vanpool_df["Industry_Type"] == "Military",
        vanpool_df["Num_Vans_Scen"] * (vanpool_df["Growth_mil"] - 1),
        vanpool_df["Num_Vans_Scen"] * (vanpool_df["Growth_non_mil"] - 1)
    )

    # round up the value
    vanpool_df["Num_Vans_Scen_ML"] = vanpool_df["Num_Vans_Scen_ML"].apply(lambda x: math.ceil(x))

    # total vanpool for the scenario year = number of vans in scenario year + number of new vans induced due to travel time saving
    vanpool_df["Num_Vans_Scen_Total"] = vanpool_df["Num_Vans_Scen"] + vanpool_df["Num_Vans_Scen_ML"]

    vanpool_df

    # summarize vanpool demand by industry type
    # vanpool demand by type
    vanpool_demand_df = vanpool_df.groupby(["Industry_Type"]).agg(
        {'Num_Vans_Base': 'sum',
        'Num_Vans_Scen': 'sum',
        'Num_Vans_Scen_ML': 'sum',
        'Num_Vans_Scen_Total': 'sum'
        }
    )

    vanpool_demand_df.reset_index(inplace=True)
    vanpool_demand_df

    #### VANPOOL VMT & GHG REDUCTIONS

    # read emission factors for the scenario year from the EMFAC outputs
    emission_df_scen = emission_df[(emission_df["Year"] == scen_year) & (emission_df["Vehicle Type"] == "Passenger Car")]
    emission_df_scen.reset_index(inplace = True, drop = True)
    emission_df_scen

    # calculate average vanpool capacity and average passengers per vanpool
    vanpool_mean_stat_df = data_df.groupby(["Industry_Type"])[["Vehicle_Capacity", "Daily_Round_Trip_Mileage", "InCounty_Round_Trip_Mileage"]].mean()
    vanpool_mean_stat_df.reset_index(inplace=True)

    # rename columns
    vanpool_mean_stat_df.rename(
        columns={"Vehicle_Capacity": "Avg_Vanpool_Capacity",
                "Daily_Round_Trip_Mileage": "Avg_RoundTrip_Mileage",
                "InCounty_Round_Trip_Mileage": "Avg_InCounty_RoundTrip_Mileage"
                },
        inplace=True
    )

    # calculate required mean stats for vanpool od data
    # The calculation removes driver from the capacity. confirm driver is included in the capacity in the first place. "ExDriver" field indicated driver excluded
    vanpool_mean_stat_df["Avg_Vanpool_Capacity_ExDriver"] = vanpool_mean_stat_df["Avg_Vanpool_Capacity"] - 1
    vanpool_mean_stat_df["Avg_Passenger_ExDriver"] = vanpool_mean_stat_df["Avg_Vanpool_Capacity_ExDriver"] * avg_vanpool_occupancy

    vanpool_mean_stat_df.head()

    # calculate VMT & GHG reduction
    co2_runex_emission_factor = emission_df_scen["CO2 RunEx Emission Factor (tons/mile)"].values[0]
    co2_strex_emission_factor = emission_df_scen["CO2 StrEx Emission Factor (tons/trip)"].values[0]

    # number of daily trips reduced = 2 trips X Averge Number of Passengers per Vanpool X Number of Vanpools [seperately for Military & Non-Military]
    # the calculation assumes that each passenger would have made a single trip if not using the vanpool service
    daily_trips_reduction = 2 * np.sum(vanpool_mean_stat_df["Avg_Passenger_ExDriver"] * vanpool_demand_df["Num_Vans_Scen_Total"])

    # daily vmt reduction (TOTAL) = Averge Number of Passengers per Vanpool X Number of Vanpools X Average roundtrip mileage [seperately for Military & Non-Military]
    daily_vmt_reduction_total = np.sum(vanpool_mean_stat_df["Avg_Passenger_ExDriver"] * vanpool_demand_df["Num_Vans_Scen_Total"] * vanpool_mean_stat_df["Avg_RoundTrip_Mileage"])

    # daily vmt reduction (TOTAL) = Averge Number of Passengers per Vanpool X Number of Vanpools X Average roundtrip mileage [seperately for Military & Non-Military]
    daily_vmt_reduction_incounty = np.sum(vanpool_mean_stat_df["Avg_Passenger_ExDriver"] * vanpool_demand_df["Num_Vans_Scen_Total"] * vanpool_mean_stat_df["Avg_InCounty_RoundTrip_Mileage"])

    # GHG reductions associated with cold start (tons)
    coldstart_ghg_reduction_tons = daily_trips_reduction * co2_strex_emission_factor

    # GHG reduction associated with VMT from running vehicles (tons)
    vmt_ghg_reduction_tons = daily_vmt_reduction_incounty * co2_runex_emission_factor

    # total GHG reductions (tons)
    daily_total_ghg_reduction_tons = coldstart_ghg_reduction_tons + vmt_ghg_reduction_tons

    # calculate total population in the scenario year based on mgra input file
    regional_population_scen = mgra_scen_input_df["pop"].sum()

    # daily ghg reduction per capita (lbs)
    daily_ghg_reduction_per_capita_lbs = (daily_total_ghg_reduction_tons * 2000) / regional_population_scen

    # TODO add calculation of daily reduction in percentage
    # GHG reductions (Ib/person) / Daily emissions per capita

    # compile outputs into a df
    regional_results_df = pd.DataFrame(data = {"Regional Population" : [regional_population_scen],
                                    "Total daily vehicle trip reduction" : [daily_trips_reduction], 
                                    "Total daily VMT reduction by vanpooling" : [daily_vmt_reduction_total], 
                                    "VMT reduced in San Diego County by vanpooling" : [daily_vmt_reduction_incounty], 
                                    "GHG reduction due to cold starts (short tons)" : [coldstart_ghg_reduction_tons], 
                                    "GHG reduction due to VMT (short tons)" : [vmt_ghg_reduction_tons],
                                    "Daily Total GHG reduction (short tons)" : [daily_total_ghg_reduction_tons],
                                    "Daily Per capita GHG reduction (lbs/person)" : [daily_ghg_reduction_per_capita_lbs]}
                                    )

    regional_results_df

    # add totals as a new row to the vanpool demand df
    new_row = pd.DataFrame({
        'Industry_Type': ['TOTAL'],
        'Num_Vans_Base': [vanpool_demand_df["Num_Vans_Base"].sum()],
        'Num_Vans_Scen': [vanpool_demand_df["Num_Vans_Scen"].sum()],
        'Num_Vans_Scen_ML': [vanpool_demand_df["Num_Vans_Scen_ML"].sum()],
        'Num_Vans_Scen_Total': [vanpool_demand_df["Num_Vans_Scen_Total"].sum()]
    })

    # Append the new row using concat
    vanpool_demand_df_export = pd.concat([vanpool_demand_df, new_row], ignore_index = True)

    # rename column names and prepare to export in new section
    vanpool_demand_df_export.rename(
        columns={"Num_Vans_Base": "Num_Vans_" + str(base_year),
                "Num_Vans_Scen": "Num_Vans_Scen_" + str(scen_year),
                "Num_Vans_Scen_ML": "Num_Vans_Scen_ML_" + str(scen_year),
                "Num_Vans_Scen_Total": "Num_Vans_Scen_Total_" + str(scen_year)
                },
        inplace = True
    )

    # print outputs
    vanpool_demand_df_export

    #### EXPORT RESULTS

    # export outputs
    # write result tables into a df
    results_dict = {
        "Regional_Results": regional_results_df,
        "Vanpool_Growth": vanpool_demand_df_export,
        "Emission_Factors": emission_df_scen
        }

    # write results to excel file
    with pd.ExcelWriter(os.path.join(output_dir, output_results_filename)) as writer:
        for key, value in results_dict.items():
            value.to_excel(writer, sheet_name = key, index = False)
            
    print_end_time(start_time)
    print("END OF SCRIPT\n")