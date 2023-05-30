import pandas as pd
import numpy as np
import yaml
import math
import os
import sys


def get_adult_population(person_df, houshold_df):
    # compute adult population by mgra

    adults_df = person_df.query("age > 15 and age < 66", engine="python")
    adults_df = pd.merge(adults_df, houshold_df, on="hhid", how="left")
    adults_mgra_df = adults_df.groupby(["mgra"])["hhid"].count().reset_index(name='adult_pop')

    return adults_mgra_df


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


def compute_regional_reductions(regional_careshare_market_df, emission_factors, carshare_vmt_reduction_factor, scenario_year):
    # calculate regional vmt and emission reduction from carshare

    co2_runex_emission_factor = emission_factors.loc[
        emission_factors.Variable == "CO2 RunEx Emission Factor (gr/mile)", "Value"].values[0]

    GRAMS_TO_SHORT_TONS = 0.0000011

    total_carshare_market = regional_careshare_market_df["total_carshare_market"].sum()

    regional_population = regional_careshare_market_df["pop"].sum()

    total_carshare_vmt_reduction = total_carshare_market * carshare_vmt_reduction_factor

    ghg_reduction = total_carshare_vmt_reduction * co2_runex_emission_factor * GRAMS_TO_SHORT_TONS

    daily_ghg_per_capita_reduction = ghg_reduction * 2000 / regional_population

    out_data = {
        "Variable": ["Regional Population",
                     "Total Carshare Market",
                     "Total VMT Reduction by Carsharing",
                     "Daily Total GHG Reduction (short tons)",
                     "Daily Per Capita GHG Reduction (lbs/person)"],
        "Value": [regional_population,
                  total_carshare_market,
                  total_carshare_vmt_reduction,
                  ghg_reduction,
                  daily_ghg_per_capita_reduction]
    }

    out_df = pd.DataFrame(out_data)

    return out_df


def compute_corridor_reductions(corridor_carshare_market_df, emission_factors, carshare_vmt_reduction_factor, scenario_year):
    # calculate vmt and emission reductions from carshare by corridor

    co2_runex_emission_factor = emission_factors.loc[
        emission_factors_df.Variable == "CO2 RunEx Emission Factor (gr/mile)", "Value"].values[0]

    GRAMS_TO_SHORT_TONS = 0.0000011

    work_df = corridor_carshare_market_df.copy()

    work_df["Scenario Year"] = scenario_year

    work_df["Total Daily VMT Reduction"] = work_df["total_carshare_market"] * \
        carshare_vmt_reduction_factor

    work_df["Daily Total GHG Reduction (short tons)"] = work_df["Total Daily VMT Reduction"] * \
        co2_runex_emission_factor * GRAMS_TO_SHORT_TONS

    work_df.rename(columns={"total_carshare_market": "Total Carshare Market"},  inplace=True)

    work_df.set_index('corridor', inplace=True)

    work_df = work_df[["Scenario Year", "Total Carshare Market",
                       "Total Daily VMT Reduction", "Daily Total GHG Reduction (short tons)"]]

    out_df = work_df.transpose()
    out_df.reset_index(inplace=True)
    out_df.rename(columns={"index": "Variable"}, inplace=True)

    return out_df

def compute_region_market_segments_reductions(regional_careshare_market_df, emission_factors, carshare_vmt_reduction_factor, scenario_year):
    # calculate vmt and emission reductions for market segements
    co2_runex_emission_factor = emission_factors.loc[
        emission_factors_df.Variable == "CO2 RunEx Emission Factor (gr/mile)", "Value"].values[0]

    GRAMS_TO_SHORT_TONS = 0.0000011

    work_df = regional_careshare_market_df.copy()
    
    work_df["Region"] = "Region"
    
    work_df["Total daily VMT reduction general population market"] = work_df["mobility_hub"] * carshare_vmt_reduction_factor

    work_df["Total daily VMT reduction college staff market"] = work_df["college_staff"] * carshare_vmt_reduction_factor

    work_df["Total daily VMT reduction college students market"] = work_df["college_student"] * carshare_vmt_reduction_factor

    work_df["Total daily VMT reduction military market"] = work_df["military_base"] * carshare_vmt_reduction_factor


    work_df["Daily Total GHG reduction general population market (short tons)"] = work_df["Total daily VMT reduction general population market"] * co2_runex_emission_factor * GRAMS_TO_SHORT_TONS

    work_df["Daily Total GHG reduction college staff market (short tons)"] = work_df["Total daily VMT reduction college staff market"] * co2_runex_emission_factor * GRAMS_TO_SHORT_TONS

    work_df["Daily Total GHG reduction college students market (short tons)"] = work_df["Total daily VMT reduction college students market"] * co2_runex_emission_factor * GRAMS_TO_SHORT_TONS

    work_df["Daily Total GHG reduction military market (short tons)"] = work_df["Total daily VMT reduction military market"] * co2_runex_emission_factor * GRAMS_TO_SHORT_TONS

    out_df = work_df.groupby(["Region"]).agg({"mobility_hub" : 'sum', "college_staff" : 'sum',
                                                    "college_student" : 'sum', "military_base" : 'sum', 
                                                    "Total daily VMT reduction general population market" : 'sum',
                                                    "Total daily VMT reduction college staff market" : 'sum',
                                                    "Total daily VMT reduction college students market" : 'sum',
                                                    "Total daily VMT reduction military market" : 'sum',
                                                    "Daily Total GHG reduction general population market (short tons)" : 'sum',
                                                    "Daily Total GHG reduction college staff market (short tons)" : 'sum',
                                                    "Daily Total GHG reduction college students market (short tons)" : 'sum',
                                                    "Daily Total GHG reduction military market (short tons)" : 'sum',
                                                     })    
    out_df.reset_index(inplace = True)
    
    out_df.rename(columns = {"mobility_hub" : "general population market", "college_staff" : "college staff market", 
                            "college_student": "college student market", "military_base" : "military market"},  inplace=True)

    out_df.set_index('Region', inplace = True)
    
    out_df = out_df.transpose()
    out_df.reset_index(inplace = True)
    out_df.rename(columns = {"index": "Variable"}, inplace = True)

    return out_df

def compute_market_segments_reductions(corridor_carshare_market_df, emission_factors, carshare_vmt_reduction_factor, scenario_year):
    # calculate vmt and emission reductions for market segements
    co2_runex_emission_factor = emission_factors.loc[
        emission_factors_df.Variable == "CO2 RunEx Emission Factor (gr/mile)", "Value"].values[0]

    GRAMS_TO_SHORT_TONS = 0.0000011

    work_df = corridor_carshare_market_df.copy()

    work_df["Total daily VMT reduction general population market"] = (work_df["mobility_hub"]
                                                                      * carshare_vmt_reduction_factor)

    work_df["Total daily VMT reduction college staff market"] = (work_df["college_staff"]
                                                                 * carshare_vmt_reduction_factor)

    work_df["Total daily VMT reduction college students market"] = (work_df["college_student"]
                                                                    * carshare_vmt_reduction_factor)

    work_df["Total daily VMT reduction military market"] = (work_df["military_base"]
                                                            * carshare_vmt_reduction_factor)

    work_df["Daily Total GHG reduction general population market (short tons)"] = (work_df["Total daily VMT reduction general population market"]
                                                                                   * co2_runex_emission_factor
                                                                                   * GRAMS_TO_SHORT_TONS)

    work_df["Daily Total GHG reduction college staff market (short tons)"] = (work_df["Total daily VMT reduction college staff market"]
                                                                              * co2_runex_emission_factor
                                                                              * GRAMS_TO_SHORT_TONS)

    work_df["Daily Total GHG reduction college students market (short tons)"] = (work_df["Total daily VMT reduction college students market"]
                                                                                 * co2_runex_emission_factor
                                                                                 * GRAMS_TO_SHORT_TONS)

    work_df["Daily Total GHG reduction military market (short tons)"] = (work_df["Total daily VMT reduction military market"]
                                                                         * co2_runex_emission_factor
                                                                         * GRAMS_TO_SHORT_TONS)

    work_df.rename(columns={"mobility_hub": "general population market",
                            "college_staff": "college staff market",
                            "college_student": "college student market",
                            "military_base": "military market"},  inplace=True)

    work_df.set_index('corridor', inplace=True)

    work_df.drop(["total_carshare_market"], axis=1, inplace=True)
    
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
    mgra_scen_input_file = config['inputs']['mgra_scen_input_file']
    mgra_base_input_file = config['inputs']['mgra_base_input_file']
    household_input_file = config['inputs']['household_input_file']
    person_input_file = config['inputs']['person_input_file']
    geography_xwalk_file = config['inputs']['geography_xwalk_file']
    emission_factors_file = config['inputs']['emission_factors_file']
    corridor_mgra_xwalk_file = config['inputs']['corridors_mgra_xwalk_file']
    carshare_mgra_file = config['inputs']['carshare_mgra_file']

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

    # read data
    mgra_scen_input_df = pd.read_csv(mgra_scen_input_file)
    mgra_base_input_df = pd.read_csv(mgra_base_input_file)
    household_input_df = pd.read_csv(household_input_file)
    person_input_df = pd.read_csv(person_input_file)
    geo_xwalk_df = pd.read_csv(geography_xwalk_file)
    emission_df = pd.read_excel(emission_factors_file)
    corridor_mgra_df = pd.read_csv(corridor_mgra_xwalk_file)
    carshare_mgra_df = pd.read_csv(carshare_mgra_file)
    emission_df = pd.read_excel(emission_factors_file)

    # preparing the main input dataset
    adults_mgra_df = get_adult_population(person_input_df, household_input_df)
    carshare_mgra_df = carshare_mgra_df[carshare_mgra_df.year == scen_year]

    data_df = pd.merge(mgra_scen_input_df, adults_mgra_df, on="mgra", how="left")
    data_df = pd.merge(data_df, carshare_mgra_df, left_on="mgra", right_on="mgra_13", how="left")

    data_df[["adult_pop", "MoHub_carshare_flag", "univ_flag", "MLB_flag"]] = data_df[[
        "adult_pop", "MoHub_carshare_flag", "univ_flag", "MLB_flag"]].fillna(0)

    data_df["pop_density"] = data_df["pop"] / data_df["acres"]
    data_df["student_enrollment"] = data_df["collegeenroll"] + data_df["othercollegeenroll"]

    data_df = data_df[["mgra", "pop", "adult_pop", "pop_density", "acres",
                       "student_enrollment", "emp_total", "MoHub_carshare_flag", "univ_flag", "MLB_flag"]]

    # build regional dataset with carshare market participations
    regional_df = data_df.copy()
    regional_df = regional_df.assign(
        MoHub_carshare_participant=[low_density_carshare_participation if a <= population_density_threshold else high_density_carshare_participation for a in regional_df['pop_density']])

    regional_df["mobility_hub"] = (regional_df["adult_pop"]
                                   * regional_df["MoHub_carshare_participant"]
                                   * regional_df["MoHub_carshare_flag"])

    regional_df["college_staff"] = (regional_df["emp_total"]
                                    * college_carshare_participation
                                    * regional_df["univ_flag"])

    regional_df["college_student"] = (regional_df["student_enrollment"]
                                      * college_carshare_participation
                                      * regional_df["univ_flag"])

    regional_df["military_base"] = (regional_df["emp_total"]
                                    * military_carshare_participation
                                    * regional_df["MLB_flag"])

    regional_df["total_carshare_market"] = (regional_df["mobility_hub"]
                                            + regional_df["college_staff"]
                                            + regional_df["college_student"]
                                            + regional_df["military_base"])

    regional_df = regional_df[["mgra", "pop", "mobility_hub", "college_staff",
                               "college_student", "military_base", "total_carshare_market"]]

    # build corridor level dataset with carshare market participations
    corridor_df = pd.merge(regional_df, corridor_mgra_df, on="mgra", how="right")
    corridor_df["mobility_hub"] = corridor_df["mobility_hub"] * corridor_df["weight"]
    corridor_df["college_staff"] = corridor_df["college_staff"] * corridor_df["weight"]
    corridor_df["college_student"] = corridor_df["college_student"] * corridor_df["weight"]
    corridor_df["military_base"] = corridor_df["military_base"] * corridor_df["weight"]

    corridor_df["total_carshare_market"] = (corridor_df["total_carshare_market"]
                                            * corridor_df["weight"])

    corridor_df = corridor_df.groupby(["corridor"]).agg({
        "mobility_hub": 'sum',
        "college_staff": 'sum',
        "college_student": 'sum',
        "military_base": 'sum',
        "total_carshare_market": 'sum'})

    corridor_df.reset_index(inplace=True)

    # get emission factors for scenario year
    emission_factors_df = get_emission_factors(emission_df, scen_year)

    # compute regional emission reduction
    regional_results_df = compute_regional_reductions(
        regional_df,
        emission_factors_df,
        daily_vmt_reduction_carshare,
        scen_year
    )

    # compute corridor level emission reduction
    corridor_results_df = compute_corridor_reductions(
        corridor_df,
        emission_factors_df,
        daily_vmt_reduction_carshare,
        scen_year
    )

    # compute market segments emission reduction
    regional_market_segments_results_df = compute_region_market_segments_reductions(
        regional_df,
        emission_factors_df,
        daily_vmt_reduction_carshare,
        scen_year
    )
    
    market_segments_results_df = compute_market_segments_reductions(
        corridor_df,
        emission_factors_df,
        daily_vmt_reduction_carshare,
        scen_year
    )
    
    market_segments_results_df = pd.merge(market_segments_results_df, regional_market_segments_results_df, on = "Variable", how = "left")





    results_dict = {"Regional_Results": regional_results_df,
                    "Corridor_Results": corridor_results_df,
                    "Market_Segment_Results": market_segments_results_df,
                    "Emission_Factors": emission_factors_df}

    write_results(results_dict, output_results_filename, output_dir)
