# %%
# Python 3.11.7

# Import libraries
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

import numpy as np
import openmatrix as omx
import yaml
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

# %%
# Read input/output files
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

    # Read config file manually
    # config_filename = r"C:/Regional Plan ABM3 OMC/carshare/data/input/config_abm3.yml"

    with open(config_filename, "r") as yml_file:
        config = yaml.safe_load(yml_file)

    mgra_scen_input_file = config['inputs']['mgra_scen_input_file']
    mgra_base_input_file = config['inputs']['mgra_base_input_file']
    household_input_file = config['inputs']['household_input_file']
    person_input_file = config['inputs']['person_input_file']
    geography_xwalk_file = config['inputs']['geography_xwalk_file']
    emission_factors_file = config['inputs']['emission_factors_file']
    carshare_file = config['inputs']['carshare_mgra_file']

    # output path
    # read config: outputs
    output_dir = config['outputs']['output_dir']
    output_results_filename = config['outputs']['output_file_name']

    # read config: parameters
    population_density_threshold = config['parameters']['population_density_threshold']
    high_density_carshare_participation = config['parameters']['high_density_mobility_hub_carshare_participation']
    low_density_carshare_participation = config['parameters']['low_density_mobility_hub_carshare_participation']
    college_carshare_participation = config['parameters']['college_carshare_participation']
    military_carshare_participation = config['parameters']['military_carshare_participation']
    base_year = config['parameters']['base_year']
    scen_year = config['parameters']['scen_year']
    daily_vmt_reduction_carshare = config['parameters']['daily_vmt_reduction_carshare']

    ###################################################################
    ########   READ INPUT/OUTPUT FILE NAMES MANUALLY   ################
    ###################################################################

    # mgra_scen_input_file = r"C:\OMC\Carshare OMC\ABM3\inputs\mgra15_based_input2035.csv"
    # mgra_base_input_file = r"C:\OMC\Carshare OMC\ABM3\inputs\mgra15_based_input2022.csv"
    # person_input_file = r"C:\OMC\Carshare OMC\ABM3\inputs\persons.csv"
    # household_input_file = r"C:\OMC\Carshare OMC\ABM3\inputs\households.csv"
    # geography_xwalk_file = r"C:\OMC\Carshare OMC\ABM3\inputs\xref_MGRA_TAZ_MSA.csv"
    # emission_factors_file = r"C:\OMC\Carshare OMC\ABM3\inputs\co2_emissions_rates.xlsx"
    # carshare_file = r"C:\OMC\Carshare OMC\ABM3\inputs\mgra_carshare_inputs_s15.csv"

    # # output file path
    # output_dir = r"C:\OMC\Carshare OMC\ABM3\outputs"
    # output_results_filename = r"carshare_calculator_results.xlsx"

    # # parameters
    # base_year = 2022
    # scen_year = 2035

    # # 2016-2017 San Diego Regional Transportation Study (SANDAG, 2017). The 2016-2017 San Diego Regional Transportation Study reports that approximately 2 percent of the San Diego population are carshare participants. In the San Diego region, coverage areas with a population density greater than 17 persons per acre are assumed to reflect these participation rates.
    # population_density_threshold = 17 

    # # carshare participation rates based on (Petersen et al, 2016). Data for the Puget Sound region indicates that carshare participation in the Seattle-Bellevue-Redmond area is 2 percent in urban neighborhoods and 0.5 percent in suburban neighborhoods. In the San Diego region, coverage areas with a population density less than 17 persons per acre are assumed to reflect the participation rates of lower density neighborhoods in the Puget Sound region.
    # high_density_carshare_participation = 0.02
    # low_density_carshare_participation = 0.005

    # # Local data on the carshare participation at colleges is unavailable. Participation rates are assumed equal to higher density area carshare participation rates
    # college_carshare_participation = 0.02

    # # Local data on the carshare participation at military bases is unavailable. Participation rates are assumed equal to higher density area carshare participation rates.
    # military_carshare_participation = 0.02

    # # Estimated based on data for San Francisco’s City CarShare service (7.0 miles per day)
    # daily_vmt_reduction_carshare = 7

    ##############################################
    ########   OPEN INPUT FILES   ################
    ##############################################

    # read data
    mgra_scen_input_df = pd.read_csv(mgra_scen_input_file)
    mgra_base_input_df = pd.read_csv(mgra_base_input_file)
    household_input_df = pd.read_csv(household_input_file)
    person_input_df = pd.read_csv(person_input_file)
    geo_xwalk_df = pd.read_csv(geography_xwalk_file)
    emission_df = pd.read_excel(emission_factors_file)
    carshare_df = pd.read_csv(carshare_file)
    emission_df = pd.read_excel(emission_factors_file)

    # %%
    # calculate adult population by mgra
    # filter adult population to 18-65 years old
    adults_df = person_input_df[(person_input_df["age"] > 18) & (person_input_df["age"] < 65)]

    # join household information to individuals in the data
    adults_df = pd.merge(adults_df, household_input_df, on = "hhid", how = "left")

    # aggregate population by mgra. The following groups individuals by mgra number, then calculate number of people based on person id in each mgra group
    adults_mgra_df = adults_df.groupby(["mgra"])["perid"].count().reset_index(name = 'adult_pop')

    adults_mgra_df.head()

    # %%
    # extract carshare numbers by mgra for the scenario year
    carshare_mgra_df = carshare_df[carshare_df["year"] == scen_year]
    carshare_mgra_df.head()

    # %%
    # join "adult population data by mgra" and "carshare data by mgra" to the input mgra for scenario year
    data_df = pd.merge(mgra_scen_input_df, adults_mgra_df, on = "mgra", how = "left")
    data_df = pd.merge(data_df, carshare_mgra_df, left_on = "mgra", right_on = "mgra_15", how = "left")

    # replace null values with zero
    data_df[["adult_pop", "MoHub_carshare_flag", "univ_flag", "MLB_flag"]] = data_df[["adult_pop", "MoHub_carshare_flag", "univ_flag", "MLB_flag"]].fillna(0)

    # calculate population density
    data_df["pop_density"] = data_df["pop"] / data_df["acres"]

    # calculate student enrollment
    data_df["student_enrollment"] = data_df["collegeenroll"] + data_df["othercollegeenroll"]

    # filter columns to data necessary for the analysis
    data_df = data_df[["mgra", "pop", "adult_pop", "pop_density", "acres", "student_enrollment", "emp_total", "MoHub_carshare_flag", "univ_flag", "MLB_flag"]]
    data_df.head()

    # %%
    # build df that includes carshare market participations
    # copy data_df into a new df called regional_df
    regional_df = data_df.copy()

    # assign carshare participation for MoHub based on population density threshold
    regional_df = regional_df.assign(MoHub_carshare_participant = [low_density_carshare_participation if ii <= population_density_threshold else high_density_carshare_participation for ii in regional_df['pop_density']])

    # calculate participating population for individual markets
    regional_df["mobility_hub"] = regional_df["adult_pop"] * regional_df["MoHub_carshare_participant"] * regional_df["MoHub_carshare_flag"] 
    regional_df["college_staff"] = regional_df["emp_total"] * college_carshare_participation * regional_df["univ_flag"] 
    regional_df["college_student"] = regional_df["student_enrollment"] * college_carshare_participation * regional_df["univ_flag"] 
    regional_df["military_base"] = regional_df["emp_total"] * military_carshare_participation * regional_df["MLB_flag"]

    # aggregate markets together
    regional_df["total_carshare_market"] = regional_df["mobility_hub"] + regional_df["college_staff"] + regional_df["college_student"] + regional_df["military_base"]

    # filter columns to data necessary for the analysis
    regional_df = regional_df[["mgra", "pop", "mobility_hub", "college_staff", "college_student", "military_base", "total_carshare_market"]]

    regional_df.head()

    # %%
    # read emission factors for the scenario year from the EMFAC outputs
    emission_df_scen = emission_df[(emission_df["Year"] == scen_year) & (emission_df["Vehicle Type"] == "Passenger Car")]
    emission_df_scen.reset_index(inplace = True, drop = True)
    emission_df_scen

    # %%
    # calculate VMT & GHG reduction
    co2_runex_emission_factor = emission_df_scen["CO2 RunEx Emission Factor (tons/mile)"].values[0]

    ###################################################
    ##########   REGIONAL RESUTLS    ##################
    ###################################################

    # regional population
    regional_population = regional_df["pop"].sum()

    # total carshare market
    total_carshare_market = regional_df["total_carshare_market"].sum()

    # total vmt reduction by carsharing
    total_carshare_vmt_reduction = total_carshare_market * daily_vmt_reduction_carshare

    # total ghg reduction in tons
    total_ghg_reduction = total_carshare_vmt_reduction * co2_runex_emission_factor

    # daily ghg reduction per capita in lb/capita (1 ton = 2000 lb)
    daily_ghg_per_capita_reduction = total_ghg_reduction * 2000 / regional_population 

    ###################################################
    ##########   BY SEGMENT RESUTLS    ################
    ###################################################

    # carshare market by segment
    college_staff_carshare_market = regional_df["college_staff"].sum()
    college_student_carshare_market = regional_df["college_student"].sum()
    military_base_carshare_market = regional_df["military_base"].sum()

    # vmt reduction by carsharing by segment
    college_staff_vmt_reduction = college_staff_carshare_market * daily_vmt_reduction_carshare
    college_student_vmt_reduction = college_student_carshare_market * daily_vmt_reduction_carshare
    military_base_vmt_reduction = military_base_carshare_market * daily_vmt_reduction_carshare

    # ghg reduction in tons by segment
    college_staff_ghg_reduction = college_staff_vmt_reduction * co2_runex_emission_factor
    college_student_ghg_reduction = college_student_vmt_reduction * co2_runex_emission_factor
    military_base_ghg_reduction = military_base_vmt_reduction * co2_runex_emission_factor

    # compile outputs into a df
    regional_results_df = pd.DataFrame(data = {"Regional Population" : [regional_population],
                                            "Total Carshare Market" : [total_carshare_market],
                                            "Total VMT Reduction by Carsharing" : [total_carshare_vmt_reduction],
                                            "Total Daily GHG Reduction (short tons)" : [total_ghg_reduction],
                                            "Daily Per Capita GHG Reduction (lbs/person)" : [daily_ghg_per_capita_reduction],
                                            "College Staff Carshare Market" : [college_staff_carshare_market],
                                            "College Student Carshare Market" : [college_student_carshare_market],
                                            "Military Base Carshare Market" : [military_base_carshare_market],
                                            "College Staff VMT Reduction by Carsharing" : [college_staff_vmt_reduction],
                                            "College Student VMT Reduction by Carsharing" : [college_student_vmt_reduction],
                                            "Military Base VMT Reduction" : [military_base_vmt_reduction],
                                            "College Staff GHG Reduction (short tons)" : [college_staff_ghg_reduction],
                                            "College Student GHG Reduction (short tons)" : [college_student_ghg_reduction],
                                            "Military Base GHG Reduction (short tons)" : [military_base_ghg_reduction]}
                                    )

    regional_results_df

    # %%
    # export outputs
    # write result tables into a df
    results_dict = {"Regional_Results": regional_results_df,
                    "Emission_Factors": emission_df}

    # write results to excel file
    with pd.ExcelWriter(os.path.join(output_dir, output_results_filename)) as writer:
        for key, value in results_dict.items():
            value.to_excel(writer, sheet_name = key, index = False)
            
            
    print_end_time(start_time)
    print("END OF SCRIPT")