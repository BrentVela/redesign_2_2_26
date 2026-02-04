#!/usr/bin/env python
import numpy as np
import pandas as pd
from tc_python import *
from itertools import compress
from tc_python import server
import time
import concurrent.futures
import os.path as path
import time
import os
import traceback
from pathlib import Path



#Pre-Processing Data
def Property(param):
    indices = param["INDICES"]
    comp_df = param["COMP"]
    elements = param["ACT_EL"]

    active_el = elements

    prev_active_el = []

    try:
        with TCPython() as session:
            session.disable_caching()

            eq_calculation = (
                session.
                    select_database_and_elements('TCHEA7', active_el).
                    get_system().
                    with_property_model_calculation("Liquidus and Solidus Temperature").
                    set_argument('upperTemperatureLimit', 4000).
                    set_composition_unit(CompositionUnit.MOLE_PERCENT).set_temperature(300))

            prop_calc = (
                session.
                    select_database_and_elements('TCHEA7', active_el).
                    get_system().
                    with_property_model_calculation('Equilibrium with Freeze-in Temperature').
                    set_temperature(25 + 273.15).
                    set_composition_unit(CompositionUnit.MOLE_PERCENT))

            indices = sorted(indices)
            print('Starting alloys with indices {}-{}'.format(indices[0],indices[-1]))
            for i in indices:
                comp = np.array(comp_df.loc[i][active_el])
                
                if len(active_el) != 1:
                    for j in range(len(active_el) - 1):
                        eq_calculation.set_composition(active_el[j], comp[j])
                        prop_calc.set_composition(active_el[j], comp[j])
                else:
                    eq_calculation.set_dependent_element(active_el[0])
                    prop_calc.set_dependent_element(active_el[0])


                result = eq_calculation.calculate()
                comp_df.at[i,'PROP LT (K)']  = result.get_value_of('Liquidus temperature')
                comp_df.at[i,'PROP ST (K)']  = result.get_value_of('Solidus temperature')

                # #Solidus Props
                # prop_calc = prop_calc.set_argument('Freeze-in-temperature', comp_df.at[i,'PROP ST (K)'])
                # prop_calc = prop_calc.set_argument('Minimization strategy', 'Global minimization only')
                # prop_calc = prop_calc.set_temperature(comp_df.at[i,'PROP ST (K)'])
                # prop_result = prop_calc.calculate()
                # comp_df.at[i,'PROP ST Density (g/cm3)'] = prop_result.get_value_of('Density (g/cm3)')
                # comp_df.at[i,'PROP ST C (J/(mol K))']     = prop_result.get_value_of('Heat capacity (J/(mol K))')
                # comp_df.at[i,'PROP ST THCD (W/(mK))'] = prop_result.get_value_of('Thermal conductivity (W/(mK))')
                # comp_df.at[i,'PROP ST THRS (mK/W)'] = prop_result.get_value_of('Thermal resistivity (mK/W)')
                # comp_df.at[i,'PROP ST TDIV (m2/s)']    = prop_result.get_value_of('Thermal diffusivity (m2/s)')
                # comp_df.at[i,'PROP ST ELRS (Ohm m)'] = prop_result.get_value_of('Electric resistivity (ohm m)')
                # comp_df.at[i,'PROP ST ELCD (S/m)'] = prop_result.get_value_of('Electric conductivity (S/m)')
                #
                # #Liquidus Props
                # prop_calc = prop_calc.set_argument('Freeze-in-temperature',comp_df.at[i,'PROP LT (K)'])
                # prop_calc = prop_calc.set_argument('Minimization strategy', 'Global minimization only')
                # prop_calc = prop_calc.set_temperature(comp_df.at[i,'PROP LT (K)'])
                # prop_result = prop_calc.calculate()
                # comp_df.at[i,'PROP LT Density (g/cm3)'] = prop_result.get_value_of('Density (g/cm3)')
                # comp_df.at[i,'PROP LT C (J/(mol K))']     = prop_result.get_value_of('Heat capacity (J/(mol K))')
                # comp_df.at[i,'PROP LT THCD (W/(mK))'] = prop_result.get_value_of('Thermal conductivity (W/(mK))')
                # comp_df.at[i,'PROP LT THRS (mK/W)'] = prop_result.get_value_of('Thermal resistivity (mK/W)')
                # comp_df.at[i,'PROP LT TDIV (m2/s)']    = prop_result.get_value_of('Thermal diffusivity (m2/s)')
                # comp_df.at[i,'PROP LT ELRS (Ohm m)'] = prop_result.get_value_of('Electric resistivity (ohm m)')
                # comp_df.at[i,'PROP LT ELCD (S/m)'] = prop_result.get_value_of('Electric conductivity (S/m)')
                #
                # #T Props, Specify T in celcius
                # temperatures = [25, 600, 1300, 2000]
                # for temp in temperatures:
                #     prop_calc = prop_calc.set_argument('Freeze-in-temperature',2000+273)
                #     prop_calc = prop_calc.set_argument('Minimization strategy', 'Global minimization only')
                #     prop_calc = prop_calc.set_argument('Reference temperature for technical CTE',25+273)
                #     prop_calc = prop_calc.set_temperature(temp+273)
                #     prop_result = prop_calc.calculate()
                #     comp_df.at[i,'PROP {}C Density (g/cm3)'.format(temp)] = prop_result.get_value_of('Density (g/cm3)')
                #     comp_df.at[i,'PROP {}C C (J/(mol K))'.format(temp)]     = prop_result.get_value_of('Heat capacity (J/(mol K))')
                #     comp_df.at[i,'PROP {}C THCD (W/(mK))'.format(temp)] = prop_result.get_value_of('Thermal conductivity (W/(mK))')
                #     comp_df.at[i,'PROP {}C THRS (mK/W)'.format(temp)] = prop_result.get_value_of('Thermal resistivity (mK/W)')
                #     comp_df.at[i,'PROP {}C TDIV (m2/s)'.format(temp)]    = prop_result.get_value_of('Thermal diffusivity (m2/s)')
                #     comp_df.at[i,'PROP {}C ELRS (Ohm m)'.format(temp)] = prop_result.get_value_of('Electric resistivity (ohm m)')
                #     comp_df.at[i,'PROP {}C ELCD (S/m)'.format(temp)] = prop_result.get_value_of('Electric conductivity (S/m)')

                
                    

    except Exception as e2:
        print('Exception occurred on line {}:'.format(traceback.extract_tb(e2.__traceback__)[0][1]))
        print(e2)

    finally:
        comp_df.to_csv('CalcFiles/PROP_OUT_{}.csv'.format(param["INDICES"][0]),index=False)
        print('Saving alloys with indices {}-{}'.format(indices[0],indices[-1]))

    return 'Complete'


if __name__ == '__main__':

    ##########################################################################
    # Initialize variables and load the data
    results_df = pd.read_excel('BrentRequest.xlsx')  # Read the CSV file into a pandas DataFrame
    elements =  sorted(['Al',	'Cr',	'Hf',	'Mo',	'Nb',	'Re',	'Ta',	'Ti',	'V',	'W',	'Zr',	'B',	'C',	'Co',	'Cu',	'Fe',	'Mn',	'N',	'Ni',	'Si',	'Sn',	'Y'])  # Define the elements of interest

    results_df = results_df[elements]

    elements = sorted(elements)  # Sort the list of elements
    savename = 'prop_out'  # Define the name for saving output
    ##########################################################################

    # Reorganize dataframe to have ordered sequence of element groups 
    # (to minimize Thermo-Calc initializations and improve efficiency)
    tic = time.time()  # Start the timer to measure time taken
    Els = list()  # Initialize a list to hold unique element combinations
    prev_active_el = []  # Store the previously active elements
    for row in range(results_df.shape[0]):  # Iterate over each row in the dataframe
        comp = results_df.iloc[row][elements]  # Get the composition for the current row
        active_el = list(compress(elements, list(comp > 0)))  # Extract elements that have a non-zero composition
        if active_el not in Els:  # If this combination of elements hasn't been added yet
            Els.append(active_el)  # Add it to the list of element combinations
        prev_active_el = active_el  # Update the previous active elements
        toc = time.time()  # Measure time taken for each iteration
        print(f"{round(row / results_df.shape[0] * 100, 3)} % Done Gathering Systems in {round(toc - tic, 3)} secs")

    
    results_df = results_df.reset_index(drop=True)  # Reset the index for the new DataFrame

    # Create a directory for saving calculation files if it doesn't exist
    if not path.exists("CalcFiles"):
        os.mkdir("CalcFiles")

    indices = results_df.index  # Get the index values from the DataFrame

    prev_active_el = []  # Initialize the previous active element tracker
    parameters = []  # List to hold calculation parameters
    count = 0  # Counter for tracking progress

    # Group the calculations into sets for efficient processing
    for i in indices:  # Loop over the rows of the DataFrame
        comp = results_df.loc[i][elements]  # Get the composition for the current row
        active_el = list(compress(elements, list(comp > 0)))  # Extract non-zero elements

        # Check if the current active elements differ from the previous ones or if the count reaches a threshold
        if (active_el != prev_active_el) or (count == 500):
            try:
                new_calc_dict["COMP"] = results_df.loc[new_calc_dict["INDICES"]]  # Assign the composition to the dictionary
                new_calc_dict["ACT_EL"] = prev_active_el  # Assign the previous active elements
                # Check if the result set already exists
                if not os.path.exists(f"CalcFiles/EQUIL-Results-Set-{new_calc_dict['INDICES'][0]}"):
                    parameters.append(new_calc_dict)  # Add the new calculation set to the parameters list
                else:
                    print(f"******Calculation Already Completed: Start Index {new_calc_dict['INDICES'][0]} \n")  # Inform the user
                new_calc_dict = {"INDICES": [], "COMP": [], "ACT_EL": []}  # Reset the calculation dictionary
            except Exception as e:
                new_calc_dict = {"INDICES": [], "COMP": [], "ACT_EL": []}  # Handle any exceptions that occur
            count = 0  # Reset the count for the next group of calculations

        new_calc_dict["INDICES"].append(i)  # Add the index to the current calculation set
        prev_active_el = active_el  # Update the previous active elements
        count += 1  # Increment the counter

    # Add the last calculation set after the loop
    new_calc_dict["COMP"] = results_df.loc[new_calc_dict["INDICES"]]  # Assign the composition to the dictionary
    new_calc_dict["ACT_EL"] = prev_active_el  # Assign the previous active elements
    if not os.path.exists(f"CalcFiles/PROP-Results-Set-{new_calc_dict['INDICES'][0]}"):
        parameters.append(new_calc_dict)  # Add the final calculation set
        print(f"**Calculation Added to list: Start Index {new_calc_dict['INDICES'][0]} \n")
    else:
        print(f"**Calculation Already Completed: Start Index {new_calc_dict['INDICES'][0]} \n")

    print("*****Calculation Sets Generated*****\n")  # Indicate that the calculation sets have been generated

    completed_calculations = []  # List to keep track of completed calculations
    del results_df  # Delete the original DataFrame as it's no longer needed

    # Use a ProcessPoolExecutor to process the calculations concurrently
    with concurrent.futures.ProcessPoolExecutor(2) as executor:
        for result_from_process in zip(parameters, executor.map(Property, parameters)):  # Map the Property function to parameters
            params, results = result_from_process  # Extract the parameters and results from the processed tasks
            if results == "Calculation Completed":
                completed_calculations.append('Completed')  # Mark the calculation as completed
