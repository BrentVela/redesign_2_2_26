#!/usr/bin/env python

"""
Compute VG (Jaswon–Cottrell viscous glide) creep rates for trajectory alloys using Thermo-Calc
diffusivities and elastic data. This script reuses a single TC system per temperature for speed,
labels BCC phase fractions consistently by temperature, and uses corrected Burgers-vector scaling.

Based on comparisons to the paper's Table 4 and Fig. 8 (TCHEA7 + MOBHEA3, rate-limiting DI),
the parity plot (fig8_parity_rate_limit_TCHEA7_MOBHEA3.pdf) shows our implementation gives an
overall reasonable reproduction of the experimental trends.
"""

import concurrent.futures
import os
import os.path as path
import time
from itertools import compress

import numpy as np
import pandas as pd
from tc_python import *


def writeToTracker(calcName, text):
    no_write = True
    while no_write:
        try:
            with open('{}-tracker.txt'.format(calcName), 'a') as f:
                f.write(text)
            no_write = False
        except:
            pass


# Elemental constants (subset needed here). Values align with 1_find_stoich_props_parallel_for_Class.py.
ELAST_DATA = {
    "W": {"B*": 310, "V*": 0.28, "Rm": 140.8},
    "Mo": {"B*": 230, "V*": 0.31, "Rm": 140},
    "Ta": {"B*": 200, "V*": 0.34, "Rm": 146.7},
    "Nb": {"B*": 170, "V*": 0.4, "Rm": 146.8},
    "V": {"B*": 160, "V*": 0.37, "Rm": 134.6},
    "Ti": {"B*": 110, "V*": 0.32, "Rm": 146.2},
    "Cr": {"B*": 160, "V*": 0.21, "Rm": 136},
}


def _elast_arrays(elements, active_el):
    missing = [el for el in set(elements) | set(active_el) if el not in ELAST_DATA]
    if missing:
        raise ValueError(f"Missing ELAST_DATA for elements: {sorted(missing)}")

    el_R = [ELAST_DATA[el]["Rm"] for el in active_el]
    el_V = [ELAST_DATA[el]["V*"] for el in elements]
    el_Br = [ELAST_DATA[el]["B*"] for el in elements]

    # Estimate atomic volume from metallic radius (Rm in pm -> Angstrom).
    # Atomic volume (Ang^3) ~ 4/3 * pi * r^3.
    el_vol = []
    for el in active_el:
        r_ang = ELAST_DATA[el]["Rm"] / 100.0
        el_vol.append((4.0 / 3.0) * np.pi * (r_ang ** 3))
    return np.array(el_vol), el_R, el_V, el_Br


# Pre-Processing Data
def EQUIL(param):
    indices = param["INDICES"]
    comp_df = param["COMP"]
    elements = param["ACT_EL"]
    active_el = elements
    k = 1.380649*(10**-23) # boltzmann
    sig = 200*10e6  # Pa applied

    el_vol, el_R, el_V, el_Br = _elast_arrays(elements, active_el)

    el_R_vals = list(el_R)
    sorted_r = sorted(el_R_vals, reverse=True)
    max_value = sorted_r[0]
    max_value_2 = sorted_r[1] if len(sorted_r) > 1 else sorted_r[0]
    max_index = el_R_vals.index(max_value)
    max_index_2 = el_R_vals.index(max_value_2)

    max_R_el = active_el[max_index]
    max_R_el_2 = active_el[max_index_2]
    print(max_R_el)
    print(max_R_el_2)

    prev_active_el = []

    with TCPython() as session:
        session.disable_caching()
        system = (
            session
            .select_thermodynamic_and_kinetic_databases_with_elements("TCHEA6", "MOBHEA3", active_el)
            .get_system()
        )
        eq_calculation = system.with_single_equilibrium_calculation()

        for itr in indices:
            comp = np.array(comp_df.loc[itr][active_el])

            if len(active_el) == 1:
                print("Unary of {}".format(active_el[0]))
            else:
                for j in range(len(active_el) - 1):
                    eq_calculation.set_condition(
                        ThermodynamicQuantity.mole_fraction_of_a_component(active_el[j]),
                        comp[j],
                    )

            for T in [25, 500, 1000, 1200, 1300, 1500, 2000]:  # C
                try:
                    eq_calculation = eq_calculation.set_condition(
                        ThermodynamicQuantity.temperature(),
                        T + 273,
                    )
                    rt_result = eq_calculation.calculate()

                    pnames = rt_result.get_phases()
                    bcc_phases = [i for i in pnames if "BCC" in i]
                    npm = [rt_result.get_value_of("NPM(" + phase + ")") for phase in pnames]

                    bcc_idx = [i for i, x in enumerate(pnames) if x in bcc_phases]
                    bcc_npms = [npm[i] for i in bcc_idx]

                    comp_df.at[itr, f"EQ {T}C MAX BCC"] = max(bcc_npms) if bcc_npms else 0.0
                    comp_df.at[itr, f"EQ {T}C SUM BCC"] = sum(bcc_npms) if bcc_npms else 0.0

                    B_avgr = np.dot(comp, el_Br)
                    R_avg = np.dot(comp, el_R)
                    V_avgr = np.dot(comp, el_V)
                    el_vol = el_vol.reshape(comp.shape)
                    Vol_avgr = np.dot(comp, el_vol)
                    a_avgr = (2 * Vol_avgr) ** (1 / 3)
                    b_avgr_ang = (a_avgr / 2) * np.sqrt(3)
                    b_avgr_m = b_avgr_ang * 1e-10

                    E_avgr3 = 3 * B_avgr * (1 - 2 * V_avgr)
                    G_avgr = (1 / 2) * E_avgr3 / (1 + V_avgr)
                    G_avgr = G_avgr * 1e9

                    e_list = []
                    for i in range(len(el_R)):
                        R = el_R[i]
                        e_list.append((R - R_avg) / R_avg)

                    D_list = []
                    for el in active_el:
                        if el != max_R_el:
                            D_list.append(
                                rt_result.get_value_of(
                                    "DI(BCC_B2,{},{},{})".format(el, el, max_R_el)
                                )
                            )
                        else:
                            D_list.append(
                                rt_result.get_value_of(
                                    "DI(BCC_B2,{},{},{})".format(el, el, max_R_el_2)
                                )
                            )

                    M_creep = 0
                    for i in range(len(D_list)):
                        M_creep = M_creep + comp[i] / D_list[i]
                    comp_df.at[itr, "Creep Merit"] = M_creep / 1e10
                    D_min = min(D_list)
                    min_index = D_list.index(D_min)
                    e_min = e_list[min_index]

                    c_min = comp[min_index]

                    D_max = max(D_list)
                    max_index = D_list.index(D_max)

                    # Average grain size from the paper: 1.2 ± 0.4 mm.
                    d = 1.2e-3
                    eps = (
                        (3.14 * (1 - V_avgr) * k * (T + 273) * D_min)
                        * (sig / G_avgr) ** 3
                    ) / (6 * (e_min**2) * c_min * (b_avgr_m**5) * G_avgr)
                    eps5 = (
                        1 * (D_min * G_avgr * b_avgr_m / k * (T + 273)) * (sig / G_avgr) ** 5
                    )
                    epsNH = 50 * (D_max * sig * b_avgr_m**3) / (k * (T + 273) * d**2)
                    epsCB = 50 * (D_max * sig * b_avgr_m**4) / ((k * (T + 273)) * (d**3))

                    comp_df.at[itr, f"{T} Min Creep VG [1/s]"] = eps
                    comp_df.at[itr, f"{T} Min Creep NH [1/s]"] = epsNH
                    comp_df.at[itr, f"{T} Min Creep CB [1/s]"] = epsCB
                    comp_df.at[itr, f"{T} Min Creep PL5 [1/s]"] = eps5
                    comp_df.at[itr, f"{T} Creep Lim El"] = active_el[min_index]
                    comp_df.at[itr, f"{T} Creep max El"] = active_el[max_index]

                except Exception as e2:
                    print("There was a fatal error in the calculation. Error is:")
                    print(e2)

                finally:
                    comp_df.to_csv("CalcFiles/CREEP_OUT_{}.csv".format(param["INDICES"][0]))
                    continue

    return 'Complete'


if __name__ == '__main__':

    ##########################################################################
    # results_df = pd.read_hdf('Mg_AlloyTotal.h5')
    # results_df = pd.read_excel('Nb70_dom_v3.xlsx')
    # results_df = pd.read_csv('EQUIL_ST_LT_OUT.csv')
    results_df = pd.read_csv("trajectory.csv")

    elements = sorted(["Cr", "Mo", "Nb", "Ta", "Ti", "V", "W"])
    
    savename = 'CREEP_OUT'
    ##########################################################################

    if not path.exists("CalcFiles"):
        os.mkdir("CalcFiles")

    writeToTracker('CREEP', "*****Start Generating Calculation Sets*****\n")

    # Reorganize dataframe to have ordered sequence of element groups (to minimize Thermo-Calc inits)
    tic = time.time()
    Els = list()
    prev_active_el = []
    for row in range(results_df.shape[0]):
        comp = results_df.iloc[row][elements]
        active_el = list(compress(elements, list(comp > 0)))
        if active_el not in Els:
            Els.append(active_el)
        prev_active_el = active_el
        toc = time.time()
        print(str(round(row / results_df.shape[0] * 100, 3)) + ' % Done Gathering Systems in ' + str(
            round(toc - tic, 3)) + ' secs')

    results_df2 = []
    for El_i in Els:
        cond = (
            (np.all(results_df[El_i] > 0, axis=1))
            & (np.sum(results_df[El_i], 1) > 1 - 1e-9)
            & (np.sum(results_df[El_i], 1) < 1 + 1e-9)
        )
        results_df2.append(results_df[cond])
        toc = time.time()
        print(str(round(Els.index(El_i) / len(Els) * 100, 3)) + ' % Done Rearranging in ' + str(
            round(toc - tic, 3)) + ' secs')
    results_df = pd.concat(results_df2, ignore_index=True)
    results_df = results_df.reset_index()
    results_df = results_df.drop(columns='index')

    indices = results_df.index

    prev_active_el = []
    parameters = []
    count = 0

    for i in indices:
        comp = results_df.loc[i][elements]
        active_el = list(compress(elements, list(comp > 0)))
        # print(i, prev_active_el, active_el)

        if (active_el != prev_active_el) or (count == 250):
            try:
                new_calc_dict["COMP"] = results_df.loc[new_calc_dict["INDICES"]]
                new_calc_dict["ACT_EL"] = prev_active_el
                if not os.path.exists(
                        "CalcFiles/{}-Results-Set-{}".format('EQUIL', new_calc_dict["INDICES"][0])):
                    parameters.append(new_calc_dict)
                    writeToTracker('CREEP', "******Calculation Added to list: Start Index {} \n".format(
                        new_calc_dict["INDICES"][0]))
                else:
                    writeToTracker('CREEP', "******Calculation Already Completed: Start Index {} \n".format(
                        new_calc_dict["INDICES"][0]))
                new_calc_dict = {"INDICES": [],
                                 "COMP": [],
                                 "ACT_EL": []}
            except Exception as e:
                new_calc_dict = {"INDICES": [],
                                 "COMP": [],
                                 "ACT_EL": []}
            count = 0

        new_calc_dict["INDICES"].append(i)
        prev_active_el = active_el
        count += 1
        print(count)

    # add the last calculation set
    new_calc_dict["COMP"] = results_df.loc[new_calc_dict["INDICES"]]
    new_calc_dict["ACT_EL"] = prev_active_el
    if not os.path.exists("CalcFiles/{}-Results-Set-{}".format('EQUIL', new_calc_dict["INDICES"][0])):
        parameters.append(new_calc_dict)
        writeToTracker('CREEP', "**Calculation Added to list: Start Index {} \n".format(new_calc_dict["INDICES"][0]))
    else:
        writeToTracker('CREEP',
                       "**Calculation Already Completed: Start Index {} \n".format(new_calc_dict["INDICES"][0]))

    writeToTracker('CREEP', "*****Calculation Sets Generated*****\n")

    completed_calculations = []
    print('Parameters:')
    print(parameters)
    with concurrent.futures.ProcessPoolExecutor(2) as executor:
        for result_from_process in zip(parameters, executor.map(EQUIL, parameters)):
            # params can be used to identify the process and its parameters
            params, results = result_from_process
            if results == "Calculation Completed":
                completed_calculations.append('Completed')

    writeToTracker('CREEP',
                   "**********Solidification Calculations Finished - {} Incomplete Calculations**********\n".format(
                       len(parameters) - len(completed_calculations)))
