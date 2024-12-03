#!/usr/bin/env python3
import os
import argparse
from constants import SOLVER_DICT
from utils import loadDict
from decimal import Decimal
from ipdb import set_trace as stop

def getNum(num):
    return float((f"{Decimal(num):.1E}").split("E")[0])


def getPwr(num):
    return int((f"{Decimal(num):.1E}").split("E")[1])


def main(args):
    solverList = ["8wave", "bouchut-split", "Roe", "HLLD", "HLLC", "bk-usm"]
    table1Data = ''''''
    table2Data = ''''''

    c_nu = 0.025
    c_nu_pos_er = 0.005
    c_nu_neg_er = -0.006
    c_nu_pos = c_nu + c_nu_pos_er
    c_nu_neg = c_nu + c_nu_neg_er

    coef = 2.3
    coef_pos_er = 0.8
    coef_neg_er = -0.5
    coef_pos = coef + coef_pos_er
    coef_neg = coef + coef_neg_er

    solverKinFit = loadDict(f"{args.dict}/kinFitDict.pkl")
    solverCurFit = loadDict(f"{args.dict}/curFitDict.pkl")
    growthDict   = loadDict(f"{args.dict}/growthDict.pkl")

    for solver in solverList:
        if solver not in solverKinFit or solver not in solverCurFit:
            continue
        kinFit = solverKinFit[solver]
        curFit = solverCurFit[solver]
        
        p_kin = -1.7
        p_nu = kinFit["p_nu"][0]
        p_nu_pos_er = kinFit["p_nu"][2]
        p_nu_neg_er = kinFit["p_nu"][1]
        p_nu_pos = p_nu + p_nu_pos_er
        p_nu_neg = p_nu + p_nu_neg_er
        p_bn = kinFit["p_bn"][0]
        p_bn_pos_er = kinFit["p_bn"][2]
        p_bn_neg_er = kinFit["p_bn"][1]
        k_bn = kinFit["k_bn"][0]
        k_bn_pos_er = kinFit["k_bn"][2]
        k_bn_neg_er = kinFit["k_bn"][1]
        k_nu_tilde = kinFit["k_nu_tilde"][0]
        k_nu_tilde_pos_er = kinFit["k_nu_tilde"][2]
        k_nu_tilde_neg_er = kinFit["k_nu_tilde"][1]
        k_nu_tilde_pos = k_nu_tilde + k_nu_tilde_pos_er
        k_nu_tilde_neg = k_nu_tilde + k_nu_tilde_neg_er
        k_nu = k_nu_tilde**(1/p_nu)
        k_nu_pos = max(k_nu_tilde_pos**(1/p_nu_neg), k_nu_tilde_pos**(1/p_nu_pos), k_nu_tilde_neg**(1/p_nu_neg), k_nu_tilde_neg**(1/p_nu_pos))
        k_nu_neg = min(k_nu_tilde_pos**(1/p_nu_neg), k_nu_tilde_pos**(1/p_nu_pos), k_nu_tilde_neg**(1/p_nu_neg), k_nu_tilde_neg**(1/p_nu_pos))
        k_nu_pos_er = k_nu_pos - k_nu
        k_nu_neg_er = k_nu_neg - k_nu

        k_eta = curFit["k_eta_23"][0]
        k_eta_pos_er = curFit["k_eta_23"][2]
        k_eta_neg_er = curFit["k_eta_23"][1]
        k_eta_pos = k_eta + k_eta_pos_er
        k_eta_neg = k_eta + k_eta_neg_er

        Re = (k_nu / (c_nu * 2))**(4/3)
        Re_pos = ((k_nu_pos / (c_nu_neg * 2))**(4/3))
        Re_neg = ((k_nu_neg / (c_nu_pos * 2))**(4/3))
        Re_pos_er = Re_pos - Re
        Re_neg_er = Re_neg - Re
        Re_pos_er_disp = round(Re_pos_er / 10**getPwr(Re), 1)
        Re_neg_er_disp = round(Re_neg_er / 10**getPwr(Re), 1)

        Pm = (k_eta / (coef * k_nu))**2
        Pm_pos = ((k_eta_pos / (coef_neg * k_nu_neg))**2)
        Pm_neg = ((k_eta_neg / (coef_pos * k_nu_pos))**2)
        Pm_pos_er = Pm_pos - Pm
        Pm_neg_er = Pm_neg - Pm
        Pm_pos_er_disp = f"{Pm_pos_er:.1f}"
        Pm_neg_er_disp = f"{Pm_neg_er:.1f}"

        Rm = Pm * Re
        Rm_pos = Pm_pos * Re_pos
        Rm_neg = Pm_neg * Re_neg
        Rm_pos_er = Rm_pos - Rm
        Rm_neg_er = Rm_neg - Rm
        Rm_pos_er_disp = round(Rm_pos_er / 10**getPwr(Rm), 1)
        Rm_neg_er_disp = round(Rm_neg_er / 10**getPwr(Rm), 1)

        tau = growthDict[solver]["alpha"][0]
        tau_neg_er = growthDict[solver]["alpha"][1]
        tau_pos_er = growthDict[solver]["alpha"][2]

        sat = growthDict[solver]["sat"][0]
        sat_neg_er = growthDict[solver]["sat"][1]
        sat_pos_er = growthDict[solver]["sat"][2]

        table1Data += f"{SOLVER_DICT[solver]} & ${tau:.3f}^{{+{tau_pos_er:.3f}}}_{{{tau_neg_er:.3f}}}$ & ${sat:.2f}^{{+{sat_pos_er:.2f}}}_{{{sat_neg_er:.2f}}}$ & ${p_bn:.2f}^{{+{p_bn_pos_er:.2f}}}_{{{p_bn_neg_er:.2f}}}$ & ${k_bn:.2f}^{{+{k_bn_pos_er:.2f}}}_{{{k_bn_neg_er:.2f}}}$ & ${k_nu:.2f}^{{+{k_nu_pos_er:.2f}}}_{{{k_nu_neg_er:.2f}}}$ &  ${k_eta:.1f}^{{+{k_eta_pos_er:.1f}}}_{{{k_eta_neg_er:.1f}}}$ \\\\ \n"
        table2Data += f"{SOLVER_DICT[solver]} & ${getNum(Re)}\\times10^{{{getPwr(Re)}}}$ & ${getNum(Rm)}\\times10^{{{getPwr(Rm)}}}$ & ${getNum(Pm)}$\\\\ \n"

    table1Str = f'''
    \\begin{{table*}}
    \\centering
    \\setlength{{\\tabcolsep}}{{1.8pt}}
    \\renewcommand{{\\arraystretch}}{{1.5}}
    \\begin{{tabular}}{{lcccccccc}}
    \\hline
    Name & $\\Gamma(t_\\mathrm{{turb}}^{{-1}})$ & $(E_\\mathrm{{M  }}/E_\\mathrm{{K}})_\\mathrm{{sat}}$ & $p_\\mathrm{{bn}}$ & $k_\\mathrm{{bn}}$ & ${{k}}_\\mathrm{{\\nu}}$ & $k_\\mathrm{{\\eta}}$\\\\
    (1) & (2) & (3) & (4) & (5) & (6) & (7)\\\\
    \\hline
    {table1Data}
    \\hline
    \\end{{tabular}}
    \\caption{{All parameters were measured/derived by averaging over the kinematic phase of the dynamo when $t>3t_\\mathrm{{turb}}$ and $10^{{-5}}\\le E_\\mathrm{{mag}}/E_\\mathrm{{kin}} \\le 10^{{-3}}$. Columns: \\textbf{{(1)}} Name of the numerical scheme as described in Table~\\ref{{tab:solvers}}. \\textbf{{(2)}} Exponent of the power law part of the kinetic spectra. \\textbf{{3}} Average value of the ratio of the magnetic energy to the kinetic energy in the saturation phase of the dynamo $(t>60t_\\mathrm{{turb}})$. \\textbf{{(4)}} Exponent of the bottleneck effect. \\textbf{{(5)}} Scaling wave number of the bottleneck effect. \\textbf{{(6)}} Viscous dissipation wave number. \\textbf{{(7)}} Resistive dissipation wave number.}}
    \\label{{tab:Turbulent dynamo fit parameters}}
    \\end{{table*}}
    '''
    table2Str = f'''
    \\begin{{table*}}
    \\centering
    \\setlength{{\\tabcolsep}}{{1.8pt}}
    \\renewcommand{{\\arraystretch}}{{1.5}}
    \\begin{{tabular}}{{lcccc}}
    \\hline
    Name & $\\mathrm{{Re}}$ & $\\mathrm{{Rm}}$ & $\\mathrm{{Pm}}$\\\\
    (1) & (2) & (3) & (4) \\\\
    \\hline
    {table2Data}
    \\hline
    \\end{{tabular}}
    \\caption{{Columns: \\textbf{{(1)}} Name of the numerical scheme as described in Table~\\ref{{tab:solvers}}. \\textbf{{(2)}} Effective Hydrodynamic Reynolds number. \\textbf{{(3)}} Effective Magnetic Reynolds number. \\textbf{{(4)}} Effective Prandtl number.}}
    \\label{{tab:Turbulent dynamo effective flow numbers}}
    \\end{{table*}}
    '''
    table1File = open(f"{args.o}/Table1.txt", "w")
    table1File.write(table1Str)
    table1File.close()
    table2File = open(f"{args.o}/Table2.txt", "w")
    table2File.write(table2Str)
    table2File.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate table creation")
    parser.add_argument("-i", type=str, help="Input Directory")
    parser.add_argument("-dict", type=str, help="Dictionary file location")
    parser.add_argument("-o", type=str, help="Output Directory")


    args = parser.parse_args()
    if args.dict is None:
        args.dict = args.i

    main(args)