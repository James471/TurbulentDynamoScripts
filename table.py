#!/usr/bin/env python3
import os
import argparse
from constants import SOLVER_DICT
from utils import loadDict
from decimal import Decimal
import numpy as np
import pickle
from ipdb import set_trace as stop

def getNum(num):
    return float((f"{Decimal(num):.1E}").split("E")[0])


def getPwr(num):
    return int((f"{Decimal(num):.1E}").split("E")[1])


base = "/scratch/ek9/jw5893/TurbulentDynamo/Mach0.1IdealSolHighRes"
dir_dct = {"8wave": base+"/SplitRoe256", "bouchut-split": base+"/SplitBouchut256", "Roe": base+"/Roe256", "HLLD": base+"/HLLD256", "HLLC": base+"/HLLC256", "bk-usm": base+"/BK256_MCut0.4_M0.1"}


def generate_random_gaussian_numbers(n=100, mu=0.0, sigma=1.0, seed=None):
    from random import seed as random_seed, gauss as random_gauss
    random_seed(seed) # set the random seed; if None, random uses the system time
    random_numbers = [random_gauss(mu, sigma) for _ in range(n)]
    return np.array(random_numbers)
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

    coef = 2.1
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
        
        p_nu = kinFit["p_nu"][0]
        p_nu_pos_er = kinFit["p_nu"][2]
        p_nu_neg_er = kinFit["p_nu"][1]
        
        p_bn = kinFit["p_bn"][0]
        p_bn_pos_er = kinFit["p_bn"][2]
        p_bn_neg_er = kinFit["p_bn"][1]
        k_bn = kinFit["k_bn"][0]
        k_bn_pos_er = kinFit["k_bn"][2]
        k_bn_neg_er = kinFit["k_bn"][1]
        k_nu_tilde = kinFit["k_nu_tilde"][0]
        k_nu_tilde_pos_er = kinFit["k_nu_tilde"][2]
        k_nu_tilde_neg_er = kinFit["k_nu_tilde"][1]
        
        with open(f"{dir_dct[solver]}/k_nu_sample.pkl", "rb") as f:
            k_nu_sample = pickle.load(f)
        with open(f"{dir_dct[solver]}/res_scl.pkl", "rb") as f:
            res_scl = pickle.load(f)

        n = len(k_nu_sample)
        k_eta_sample = np.random.choice(res_scl, n)

        k_nu = np.nanpercentile(k_nu_sample, 50)
        k_nu_pos_er = np.nanpercentile(k_nu_sample, 84) - k_nu
        k_nu_neg_er = k_nu - np.nanpercentile(k_nu_sample, 16)

        k_eta = curFit["k_eta_23"][0]
        k_eta_pos_er = curFit["k_eta_23"][2]
        k_eta_neg_er = curFit["k_eta_23"][1]
        print("Using samples:", n)
        c_nu_sample = generate_random_gaussian_numbers(n, 0.025, 0.006)
        with open("/home/100/jw5893/LowMach/LMDsync/coef_dist.pkl", "rb") as f:
            coef_sample = pickle.load(f)
        coef_sample = np.random.choice(coef_sample, n)
        Re_sample = (k_nu_sample / (c_nu_sample * 2))**(4/3)
        Re = np.nanpercentile(Re_sample, 50)
        Re_pos_er = np.nanpercentile(Re_sample, 84) - Re
        Re_neg_er = np.nanpercentile(Re_sample, 16) - Re
        # stop()
        Re_pwr = getPwr(Re)
        Re_pos_er_disp = f"{Re_pos_er/10**Re_pwr:.1f}"
        Re_neg_er_disp = f"{Re_neg_er/10**Re_pwr:.1f}"

        Pm_sample = (k_eta_sample / (coef_sample * k_nu_sample))**2
        Pm = np.nanpercentile(Pm_sample, 50)
        Pm_pos_er = np.nanpercentile(Pm_sample, 84) - Pm
        Pm_neg_er = np.nanpercentile(Pm_sample, 16) - Pm
        
        Rm_sample = Re_sample * Pm_sample
        Rm = np.nanpercentile(Rm_sample, 50)
        Rm_pos_er = np.nanpercentile(Rm_sample, 84) - Rm
        Rm_neg_er = np.nanpercentile(Rm_sample, 16) - Rm
        Rm_pwr = getPwr(Rm)
        Rm_pos_er_disp = f"{Rm_pos_er/10**Rm_pwr:.1f}"
        Rm_neg_er_disp = f"{Rm_neg_er/10**Rm_pwr:.1f}"

        tau = growthDict[solver]["alpha"][0]
        tau_neg_er = growthDict[solver]["alpha"][1]
        tau_pos_er = growthDict[solver]["alpha"][2]

        sat = growthDict[solver]["sat"][0]
        sat_neg_er = growthDict[solver]["sat"][1]
        sat_pos_er = growthDict[solver]["sat"][2]

        table1Data += f"{SOLVER_DICT[solver]} & ${tau:.2f}^{{+{tau_pos_er:.2f}}}_{{{tau_neg_er:.2f}}}$ & ${sat:.2f}^{{+{sat_pos_er:.2f}}}_{{{sat_neg_er:.2f}}}$ & ${p_bn:.2f}^{{+{p_bn_pos_er:.2f}}}_{{{p_bn_neg_er:.2f}}}$ & ${k_bn:.1f}^{{+{k_bn_pos_er:.1f}}}_{{{k_bn_neg_er:.1f}}}$ & ${k_nu_tilde:.1f}^{{+{k_nu_tilde_pos_er:.1f}}}_{{{k_nu_tilde_neg_er:.1f}}}$ & ${p_nu:.1f}^{{+{p_nu_pos_er:.1f}}}_{{{p_nu_neg_er:.1f}}}$ & ${k_nu:.1f}^{{+{k_nu_pos_er:.1f}}}_{{{-k_nu_neg_er:.1f}}}$ & ${k_eta:.0f}^{{+{k_eta_pos_er:.0f}}}_{{{k_eta_neg_er:.0f}}}$ \\\\ \n"
        table2Data += f"{SOLVER_DICT[solver]} & ${getNum(Re)}^{{+{Re_pos_er_disp}}}_{{{Re_neg_er_disp}}}\\times10^{getPwr(Re)}$ & ${getNum(Rm)}^{{+{Rm_pos_er_disp}}}_{{{Rm_neg_er_disp}}}\\times10^{getPwr(Rm)}$ & ${Pm:.1f}^{{+{Pm_pos_er:.1f}}}_{{{Pm_neg_er:.1f}}}$\\\\ \n"

    table1Str = f'''
    {table1Data}
    '''
    table2Str = f'''
    {table2Data}
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