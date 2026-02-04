import os
import sys
import pandas as pd
import numpy as np
from tc_python import TCPython, CompositionUnit, ThermodynamicQuantity, LoggingPolicy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT_LANDING = os.path.join(ROOT, "root_landing")
if ROOT_LANDING not in sys.path:
    sys.path.insert(0, ROOT_LANDING)
import strength_model

PROP_DATA = {
    'W': {'Density [g/cm^3]': 19.25, 'Melting Temperature [K]': 3695, 'Valence Electrons': 6},
    'Mo': {'Density [g/cm^3]': 10.28, 'Melting Temperature [K]': 2896, 'Valence Electrons': 6},
    'Ta': {'Density [g/cm^3]': 16.69, 'Melting Temperature [K]': 3290, 'Valence Electrons': 5},
    'Nb': {'Density [g/cm^3]': 8.57, 'Melting Temperature [K]': 2750, 'Valence Electrons': 5},
    'V': {'Density [g/cm^3]': 6.11, 'Melting Temperature [K]': 2183, 'Valence Electrons': 5},
    'Ti': {'Density [g/cm^3]': 4.506, 'Melting Temperature [K]': 1941, 'Valence Electrons': 4},
    'Zr': {'Density [g/cm^3]': 6.52, 'Melting Temperature [K]': 2128, 'Valence Electrons': 4},
    'Cr': {'Density [g/cm^3]': 7.15, 'Melting Temperature [K]': 2180, 'Valence Electrons': 6},
}

ELAST_DATA = {
    'W': {'C11': 517.8, 'C12': 201.7, 'C44': 139.4, 'B*': 310, 'G*': 161, 'V*': 0.28},
    'Mo': {'C11': 466, 'C12': 165.2, 'C44': 99.5, 'B*': 230, 'G*': 20, 'V*': 0.31},
    'Ta': {'C11': 260.9, 'C12': 165.2, 'C44': 70.4, 'B*': 200, 'G*': 67, 'V*': 0.34},
    'Nb': {'C11': 247.2, 'C12': 140, 'C44': 14.2, 'B*': 170, 'G*': 38, 'V*': 0.4},
    'V': {'C11': 272, 'C12': 144.8, 'C44': 17.6, 'B*': 160, 'G*': 47, 'V*': 0.37},
    'Ti': {'C11': 95.9, 'C12': 115.9, 'C44': 40.3, 'B*': 110, 'G*': 44, 'V*': 0.32},
    'Zr': {'C11': 81.8, 'C12': 94.3, 'C44': 30.2, 'B*': 91.1, 'G*': 33, 'V*': 0.34},
    'Cr': {'C11': 247.6, 'C12': 73.4, 'C44': 48.3, 'B*': 160, 'G*': 115, 'V*': 0.21},
}


def ensure_tc_home() -> None:
    if os.environ.get('TC25B_HOME'):
        return
    default = '/opt/Thermo-Calc/2025b'
    if os.path.isdir(default):
        os.environ['TC25B_HOME'] = default
    else:
        raise EnvironmentError("Environment variable 'TC25B_HOME' not found and default path is missing.")


def compute_avg_props(df: pd.DataFrame, elements: list[str], test_temp_c: int = 1300) -> pd.DataFrame:
    el_density = np.array([PROP_DATA[e]['Density [g/cm^3]'] for e in elements])
    el_tm = np.array([PROP_DATA[e]['Melting Temperature [K]'] for e in elements])
    el_ve = np.array([PROP_DATA[e]['Valence Electrons'] for e in elements])
    el_c11 = np.array([ELAST_DATA[e]['C11'] for e in elements])
    el_c12 = np.array([ELAST_DATA[e]['C12'] for e in elements])
    el_c44 = np.array([ELAST_DATA[e]['C44'] for e in elements])
    el_br = np.array([ELAST_DATA[e]['B*'] for e in elements])
    el_vr = np.array([ELAST_DATA[e]['V*'] for e in elements])

    comp = df[elements].values

    df['Density Avg'] = comp.dot(el_density)
    df['Tm Avg'] = comp.dot(el_tm)
    df['VEC Avg'] = comp.dot(el_ve)

    c11_avg = comp.dot(el_c11)
    c12_avg = comp.dot(el_c12)
    c44_avg = comp.dot(el_c44)
    df['Cauchy Pres Avg'] = c12_avg - c44_avg

    b_avgr = comp.dot(el_br)
    v_avgr = comp.dot(el_vr)
    e_avgr3 = 3 * b_avgr * (1 - 2 * v_avgr)
    g_avgr3 = 0.5 * e_avgr3 / (1 + v_avgr)
    df['Pugh_Ratio_PRIOR'] = b_avgr / g_avgr3

    ys_t = []
    for i in range(len(df)):
        comp_i = df.iloc[i][elements].values
        _, _, tau = strength_model.model_Control(elements, comp_i, T=test_temp_c + 273)
        ys_t.append(3000 * tau)
    df['YS T C PRIOR'] = ys_t

    return df


def compute_tc_props(
    df: pd.DataFrame,
    elements: list[str],
    temps_c: tuple[int, int] = (600, 1300),
    freeze_in_c: int = 2000,
) -> pd.DataFrame:
    ensure_tc_home()
    active = elements

    bcc_cols: dict[str, list[float]] = {}
    for t in temps_c:
        bcc_cols[f'EQ {t}C BCC_B2 MOL'] = []
        bcc_cols[f'EQ {t}C BCC_B2#2 MOL'] = []

    solidus = []
    liquidus = []
    thcd_25 = []

    with TCPython(logging_policy=LoggingPolicy.NONE) as session:
        session.disable_caching()
        system = session.select_database_and_elements('TCHEA7', active).get_system()

        ls_calc = (
            system
            .with_property_model_calculation('Liquidus and Solidus Temperature')
            .set_argument('upperTemperatureLimit', 4000)
            .set_composition_unit(CompositionUnit.MOLE_PERCENT)
            .set_temperature(300)
        )

        eq = system.with_single_equilibrium_calculation().set_condition(
            ThermodynamicQuantity.temperature(), 298
        )
        eq.disable_global_minimization()
        eq.set_phase_to_suspended('*')
        for phase in ['BCC_B2', 'BCC_B2#2', 'C14_LAVES', 'C15_LAVES', 'C36_LAVES']:
            eq.set_phase_to_entered(phase)

        thcd_calc = (
            system
            .with_property_model_calculation('Equilibrium with Freeze-in Temperature')
            .set_composition_unit(CompositionUnit.MOLE_PERCENT)
            .set_temperature(25 + 273.15)
        )

        for _, row in df.iterrows():
            comp = np.array([row[e] for e in active]) * 100.0
            for j, el in enumerate(active[:-1]):
                ls_calc.set_composition(el, comp[j])
                eq.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(el), comp[j] / 100.0)
                thcd_calc.set_composition(el, comp[j])

            res = ls_calc.calculate()
            liquidus.append(res.get_value_of('Liquidus temperature'))
            solidus.append(res.get_value_of('Solidus temperature'))

            for t in temps_c:
                res_eq = eq.set_condition(ThermodynamicQuantity.temperature(), t + 273).calculate()
                b2 = res_eq.get_value_of('NPM(BCC_B2)') if 'BCC_B2' in res_eq.get_phases() else 0.0
                b22 = res_eq.get_value_of('NPM(BCC_B2#2)') if 'BCC_B2#2' in res_eq.get_phases() else 0.0
                bcc_cols[f'EQ {t}C BCC_B2 MOL'].append(b2)
                bcc_cols[f'EQ {t}C BCC_B2#2 MOL'].append(b22)

            thcd_calc = thcd_calc.set_argument('Freeze-in-temperature', freeze_in_c + 273.15)
            thcd_calc = thcd_calc.set_argument('Minimization strategy', 'Global minimization only')
            thcd_calc = thcd_calc.set_argument('Reference temperature for technical CTE', 25 + 273.15)
            res_th = thcd_calc.calculate()
            thcd_25.append(res_th.get_value_of('Thermal conductivity (W/(mK))'))

    df['PROP LT (K)'] = liquidus
    df['PROP ST (K)'] = solidus
    df['Solidification Range (K)'] = df['PROP LT (K)'] - df['PROP ST (K)']
    for k, v in bcc_cols.items():
        df[k] = v
    df['BCC_600C_TOTAL_MOL'] = df['EQ 600C BCC_B2 MOL'].fillna(0) + df['EQ 600C BCC_B2#2 MOL'].fillna(0)
    df['BCC_1300C_TOTAL_MOL'] = df['EQ 1300C BCC_B2 MOL'].fillna(0) + df['EQ 1300C BCC_B2#2 MOL'].fillna(0)
    df['Thermal Conductivity 25C (W/mK)'] = thcd_25
    return df


def process(input_csv: str, out_csv: str, alloy_prefix: str) -> None:
    df = pd.read_csv(input_csv)
    for col in ['Cr', 'Mo', 'Nb', 'Ta', 'Ti', 'V', 'W', 'Zr']:
        if col not in df.columns:
            df[col] = 0.0
    df['Alloy_ID'] = [f"{alloy_prefix}_{i:03d}" for i in range(len(df))]

    elements = ['Cr', 'Mo', 'Nb', 'Ta', 'Ti', 'V', 'W', 'Zr']
    df = compute_avg_props(df, elements)
    df = compute_tc_props(df, elements)

    props = [
        'Density Avg', 'Tm Avg', 'VEC Avg', 'Cauchy Pres Avg', 'Pugh_Ratio_PRIOR', 'YS T C PRIOR',
        'BCC_600C_TOTAL_MOL', 'BCC_1300C_TOTAL_MOL', 'Solidification Range (K)',
        'Thermal Conductivity 25C (W/mK)'
    ]
    out_cols = ['Alloy_ID'] + elements + props
    df[out_cols].to_csv(out_csv, index=False)
    print(f'Wrote {out_csv} with {len(df)} rows')


def main() -> None:
    process('trajectory_search/search_results.csv', 'trajectory_search/properties_list.csv', 'TS')
    process('trajectory_search_grid2/search_results.csv', 'trajectory_search_grid2/properties_list.csv', 'TSG2')


if __name__ == '__main__':
    main()
