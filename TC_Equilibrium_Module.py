#!/usr/bin/env python
import os
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


def _ensure_tc_home() -> None:
    if os.environ.get("TC25B_HOME"):
        return
    default_path = "/opt/Thermo-Calc/2025b"
    if os.path.isdir(default_path):
        os.environ["TC25B_HOME"] = default_path
    else:
        raise EnvironmentError("Environment variable 'TC25B_HOME' not found and default path is missing.")


_ensure_tc_home()


#Pre-Processing Data
def EQUIL(param, temps_c):
    indices = param["INDICES"]
    comp_df = param["COMP"]
    elements = param["ACT_EL"]
    active_el = elements

    with TCPython(logging_policy=LoggingPolicy.NONE) as session:
        session.disable_caching()

        system = (
            session.
                select_database_and_elements('TCHEA7', active_el).
                get_system())

        options = SingleEquilibriumOptions().set_required_accuracy(1.0e-2).set_max_no_of_iterations(300).set_smallest_fraction(1.0e-6)
        eq_calculation = system.with_single_equilibrium_calculation(). \
            set_condition(ThermodynamicQuantity.temperature(), 298).with_options(options)
        eq_calculation.disable_global_minimization()

        # Focus on BCC and Laves phases only for faster convergence
        eq_calculation.set_phase_to_suspended(ALL_PHASES)
        for phase in ["BCC_B2", "BCC_B2#2", "C14_LAVES", "C15_LAVES", "C36_LAVES"]:
            eq_calculation.set_phase_to_entered(phase)

        for i in indices:
            # Get the composition list and corresponding active_el list
            solidus = comp_df.loc[i]['PROP ST (K)'] if 'PROP ST (K)' in comp_df.columns else None
            liquidus = comp_df.loc[i]['PROP LT (K)'] if 'PROP LT (K)' in comp_df.columns else None
            comp = np.array(comp_df.loc[i][active_el])
            #Create equilibrium calculation object and set conditions
            try: # with TCAL7, if that fails then we will try another database

                # Check Point 1
                if len(active_el) == 1:
                    # Do not do anything for unaries
                    eq_calculation.set_dependent_element(active_el[0])
                else:
                    for j in range(len(active_el)-1):
                        eq_calculation.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(active_el[j]),
                                                     comp[j])

                for temp in temps_c:
                    eq_calculation =  eq_calculation.set_condition(ThermodynamicQuantity.temperature(), temp+273)
                    eq_result = eq_calculation.calculate()

                    #Get all possible phases
                    pnames = eq_result.get_phases()

                    for phase in pnames:
                        #Get mol fraction list for each phase
                        if eq_result.get_value_of('NPM(' + phase + ')') > 0: #Only report phases that are present
                            comp_df.at[i,'EQ {}C {} MOL'.format(temp,phase)] = eq_result.get_value_of('NPM(' + phase + ')')


                if solidus is not None:
                    # Properties at Solidus
                    eq_calculation = eq_calculation.set_condition(ThermodynamicQuantity.temperature(), solidus-1)
                    eq_result = eq_calculation.calculate()
                    comp_df.at[i, 'EQ ST H (J/mol)']     = eq_result.get_value_of('HM')
                    comp_df.at[i, 'EQ ST H (J)']         = eq_result.get_value_of('H')
                    comp_df.at[i, 'EQ ST THCD (W/mK)']   = eq_result.get_value_of('THCD')
                    comp_df.at[i, 'EQ ST Density (g/cc)'] =eq_result.get_value_of('BM') / eq_result.get_value_of('VM') / 10 ** 6
                    comp_df.at[i, 'EQ ST MASS (g/mol)'] = eq_result.get_value_of('BM')
                    comp_df.at[i, 'EQ ST VOL (m3/mol)'] = eq_result.get_value_of('VM')

                if liquidus is not None:
                    # Properties at Liquidus
                    eq_calculation = eq_calculation.set_condition(ThermodynamicQuantity.temperature(), liquidus+1)
                    eq_result = eq_calculation.calculate()
                    comp_df.at[i, 'EQ LT H (J/mol)']    = eq_result.get_value_of('HM')  # J/mol
                    comp_df.at[i, 'EQ LT H (J)']        = eq_result.get_value_of('H')
                    comp_df.at[i, 'EQ LT THCD (W/mK)']  = eq_result.get_value_of('THCD')  # W/mK
                    comp_df.at[i, 'EQ LT DVIS (Pa-s)']  = eq_result.get_value_of('DVIS (liquid)')  # Pa-s
                    comp_df.at[i, 'EQ LT KVIS (m2/s)']  = eq_result.get_value_of('KVIS (liquid)')
                    comp_df.at[i, 'EQ LT Density (g/cc)'] = eq_result.get_value_of('BM') / eq_result.get_value_of('VM') / 10 ** 6
                    comp_df.at[i, 'EQ LT VOL (m3/mol)'] = eq_result.get_value_of('VM')
                    comp_df.at[i, 'EQ LT Surface Tension (N/m)'] = eq_result.get_value_of('SURF(LIQUID)')

            except Exception as e2:
                print('Exception occurred on line {}:'.format(traceback.extract_tb(e2.__traceback__)[0][1]))
                print(e2)
                
            finally:
                comp_df.to_csv('CalcFiles/EQUIL_AP_OUT_{}.csv'.format(param["INDICES"][0]))
                continue

    return 'Complete'


def _equil_with_temps(args):
    param, temps_c = args
    return EQUIL(param, temps_c)


if __name__ == '__main__':

    ##########################################################################
    # Initialize variables and load the data
    results_df = pd.read_csv('trajectory.csv')
    elements = ['Cr','Mo','Nb','Ta','Ti','V','W']
    elements = sorted(elements)
    temps_c = [500, 600, 1300, 2000]
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
        if 'PROP LT (K)' in results_df.columns and 'PROP ST (K)' in results_df.columns:
            comp = results_df.loc[i][elements + ['PROP LT (K)', 'PROP ST (K)']]
        else:
            comp = results_df.loc[i][elements]
        active_el = list(compress(elements, list(comp > 0)))  # Extract non-zero elements

        # Check if the current active elements differ from the previous ones or if the count reaches a threshold
        if (active_el != prev_active_el) or (count == 5):
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
        for result_from_process in zip(parameters, executor.map(_equil_with_temps, [(p, temps_c) for p in parameters])):
            params, results = result_from_process  # Extract the parameters and results from the processed tasks
            if results == "Calculation Completed":
                completed_calculations.append('Completed')  # Mark the calculation as completed
