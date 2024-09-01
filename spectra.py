#!/usr/bin/env python3
import sys
import os
import glob
import socket

import matplotlib
matplotlib.use('Agg')

import numpy as np  
import argparse
import matplotlib.pyplot as pl
import pandas as pd
from scipy.special import kn
from scipy.optimize import fminbound
from scipy.interpolate import CubicSpline

from constants import *
from datautils import *
from utils import *
from myconfig import *
import designParams
sys.path.append(PYTHON_PATH)
import cfpack as cfp
import flashlib as fl
import turblib as tl

from ipdb import set_trace as stop

def Re_fun(k_nu, c_nu):
    return (k_nu/(c_nu * 2))**(4/3)


def k_eta_fun(k_tilde_eta, p_eta):
    return (k_tilde_eta)**(1/p_eta)


def Rm_fun(k_nu, c_nu, k_tilde_eta, p_eta, c_eta):
    k_eta = k_eta_fun(k_tilde_eta, p_eta)
    Pm = (k_eta/(c_eta*k_nu))**2
    return Pm * Re_fun(k_nu, c_nu)


def Log10_P_kin(k, A_kin, p_bn, k_bn, k_tilde_nu, p_nu):
    p_kin = -1.7
    return np.log10(A_kin * ((k / k_bn)**p_kin + (k / k_bn)**p_bn) * np.exp(-(k / k_tilde_nu)**p_nu))


def Log10_P_mag(k, A_mag, p_mag, p_eta, k_tilde_eta):
    y = A_mag * (k**p_mag) * kn(0, (k / k_tilde_eta)**p_eta)
    return np.log10(y)


def generateSpectra(simDir, verbose, spectType, infoDict, lf, uf, stf, nProcs=1):
    '''
    verbsoe  : Could be 0 or 1 or 2
    spectType: Could be 'vels' or 'curr' or 'mags'
    '''

    if verbose == 2:
        datFileVerbose = 1
        devNull = ""
    else:
        datFileVerbose = 0
        devNull = " > /dev/null 2>&1"

    spectraDir = simDir + "/spectra"
    if not os.path.isdir(spectraDir):
        print("creating spectra dir: '"+spectraDir+"'")
        os.mkdir(spectraDir)

    files = sorted(glob.glob(simDir+"/Turb_hdf5_plt_cnt_????"))

    t_turb = 1 / (2 * infoDict['v'])
    if not os.path.exists(simDir+"/Turb.dat_cleaned"):
        if verbose > 0: print("Loading turb.dat for cleaning")
        datfile_obj = fl.datfile(simDir+"/Turb.dat", verbose=datFileVerbose)
        datfile_obj.write_cleaned()
        if verbose > 0: print("Turb.dat cleaned")
    else:
        print("Turb.dat_cleaned already exists")

    n = getNForTurbDat(simDir, res=5e-2, filename="Turb.dat_cleaned")
    if verbose > 0: print("Loading cleaned turb.dat")
    dat = loadFile(simDir + "/Turb.dat_cleaned", shift=10, n=n)
    if verbose > 0: print("Loaded Turb.dat_cleaned")

    emag_over_ekin = dat[E_MAG_COLUMN_INDEX].flatten() / dat[E_KIN_COLUMN_INDEX].flatten()
    t = dat[TIME_COLUMN_INDEX].flatten()
    ind = (emag_over_ekin >= lf) & (emag_over_ekin <= uf) & (t/t_turb >= stf)
    t_valid = t[ind]

    dumpStart = fl.get_dump(files, time=t_valid[0])
    dumpNStart = int(dumpStart[0][-4:])
    dumpEnd = fl.get_dump(files, time=t_valid[-1])
    dumpNEnd = int(dumpEnd[0][-4:])

    if verbose > 0: 
        print("Start dump:", dumpNStart)
        print("End dump:", dumpNEnd)

    if spectType == "cur":
        types = "0 -dsets curx cury curz"
        for x in range(dumpNStart, dumpNEnd+1):
            if verbose > 0: print("Adding current variable to dump:", x)
            f = simDir+"/Turb_hdf5_plt_cnt_{:04d}".format(x)
            print("In current spectra")
            if "nid" in socket.gethostname():
                derivCmd = f"srun -n {nProcs} /software/projects/pawsey0810/jwatt/flash-tools/tools/derivative_var/derivative_var {f} -current {devNull}"
            else:
                derivCmd = f"mpirun -np {nProcs} /home/100/jw5893/flash-tools/tools/derivative_var/derivative_var {f} -current {devNull}"
            os.system(derivCmd)
    elif spectType == "vels":
        types = "1"
    elif spectType == "mags":
        types = "2"

    if spectType == "cur":
        spectFileString = "dset_curx_cury_curz"
    elif spectType == "mags":
        spectFileString = "mags"
    elif spectType == "vels":
        spectFileString = "vels"

    for x in range(dumpNStart, dumpNEnd+1):
        if verbose > 0: print("Processing dump to generate spectra:", x)
        f = simDir+"/Turb_hdf5_plt_cnt_{:04d}".format(x)
        if "nid" in socket.gethostname():
            spectrCmd = f"srun -n {nProcs} /software/projects/pawsey0810/jwatt/flash-tools/tools/spectra_mpi/spectra_mpi {f} -types {types} {devNull}"
        else:
            spectrCmd = f"mpirun -np {nProcs} /home/100/jw5893/flash-tools/tools/spectra_mpi/spectra_mpi {f} -types {types} {devNull}"
        print(f"Running: {spectrCmd}")
        os.system(spectrCmd)
    mvCmd = "mv "+simDir+"/Turb_hdf5_plt_cnt_????_spect_*.dat "+spectraDir
    os.system(mvCmd)

    spectFiles = [simDir+"/spectra/Turb_hdf5_plt_cnt_{:04d}".format(x)+"_spect_"+spectFileString+".dat" for x in range(dumpNStart, dumpNEnd+1)]
    if spectType == "vels": normalise = False
    if spectType == "mags" or spectType == "cur": normalise = True # normalise magnetic spectra and current spectra before averaging, because of field growth with time
    averDat, headerAver = tl.aver_spect(spectFiles, normalise=normalise, verbose=0)
    outfile = simDir+"/spectra/aver_spect_"+spectFileString+".dat"
    tl.write_spect(outfile, averDat, headerAver, verbose=verbose)


def plotSpectra(simDir, verbose, spectType, fact, infoDict, params, outdir, compensate=False, fit=True):
    '''
    verbsoe  : Could be 0 or 1
    spectType: Could be 'vels' or 'mags' or 'curr'
    '''

    if spectType == "mags":
        spectVar = "mags"
    elif spectType == "vels":
        spectVar = "vels"
    elif spectType == "cur":
        spectVar = "dset_curx_cury_curz"

    averFile = simDir + f"/spectra/aver_spect_{spectVar}.dat"
    averDf   = pd.read_csv(averFile, sep="\s+", header=0, skiprows=5)
    
    log10PTot = averDf['log10_Ptot']
    deltaLog10PTot = averDf['sigma_log10_Ptot']
    pTot = 10**log10PTot
    k = averDf['k']
    err = np.array([pTot * np.log(10) * deltaLog10PTot, pTot * np.log(10) * deltaLog10PTot])
    
    kFitMin = 3
    kFitMax = max(k) / 2
    kFit = np.array(k[(k >= kFitMin) & (k <= kFitMax)].tolist())
    log10PTotFit = np.array(log10PTot[(k >= kFitMin) & (k <= kFitMax)].tolist())
    deltaLog10PTotFit = np.array(deltaLog10PTot[(k >= kFitMin) & (k <= kFitMax)].tolist())

    if compensate and spectType == "vels":
        compensateFact = k**1.7
        compensateFitFact = kFit**1.7

    else:
        compensateFact = 1
        compensateFitFact = 1

    if compensate and spectType == "vels":
        plObj = cfp.plot(compensateFact * fact * pTot, k, yerr=fact * err * np.array([compensateFact, compensateFact]), shaded_err=True, label=SOLVER_DICT[infoDict['solver']], color=COLOR_DICT[infoDict['solver']])
    else:
        plObj = cfp.plot(compensateFact * fact * pTot, k, yerr=fact * err, shaded_err=True, label=SOLVER_DICT[infoDict['solver']], color=COLOR_DICT[infoDict['solver']])

    if fit:
        fitDict = {}
        fitDict = fit_func(spectType, simDir, kFit, log10PTotFit, deltaLog10PTotFit, params, compensateFitFact, fact, fitDict, verbose, infoDict)
    else:
        if spectType == "mags":
            fitDict = loadDict(f"{outdir}/magFitDict.pkl")[infoDict['solver']]
            # cfp.plot(compensateFitFact * fact * 10**Log10_P_mag(kFit, *(fitDict["A_mag"][0], fitDict["p_mag"][0], fitDict["p_eta"][0], fitDict["k_tilde_eta"][0])), kFit, color="black")
        elif spectType == "vels":
            fitDict = loadDict(f"{outdir}/kinFitDict.pkl")[infoDict['solver']]
            cfp.plot(compensateFitFact * fact * 10**Log10_P_kin(kFit, *(fitDict["A_kin"][0], fitDict["p_bn"][0], fitDict["k_bn"][0], fitDict["k_tilde_nu"][0]), fitDict["p_nu"][0]), kFit, color="black")
        elif spectType == "cur":
            fitDict = loadDict(f"{outdir}/curFitDict.pkl")[infoDict['solver']]


    return plObj, fitDict


def fit_func(spectType, simDir, kFit, log10PTotFit, deltaLog10PTotFit, params, compensateFitFact, fact, fitDict, verbose, infoDict):

    if spectType == "mags":
        # if verbose: print("Mag Spectra-> Fitting for:", infoDict['solver'])
        # magFit = cfp.fit(Log10_P_mag, kFit, log10PTotFit, xerr=None, yerr=deltaLog10PTotFit,
                            # params=params, max_nfev=100000, n_random_draws=1000)
        # cfp.plot(compensateFitFact * fact * 10**Log10_P_mag(kFit, *(magFit.popt)), kFit, color="black")
        fitDict["A_mag"] = (0,0,0)
        fitDict["p_mag"] = (0,0,0)
        fitDict["p_eta"] = (0,0,0)
        fitDict["k_tilde_eta"] = (0,0,0)

    elif spectType == "cur":
        if verbose: print("Current Spectra-> Working with:", infoDict['solver'])
        infoDict = getInfoDict(simDir)
        spectraDir = simDir+"/spectra"
        kList = []
        for file in os.listdir(spectraDir):
            if "curx" in file and "aver" not in file:
                spectraFile = spectraDir + "/" + file
                plotFileNum = file.split("_")[4].split("_")[0]
                spectraDf   = pd.read_csv(spectraFile, sep="\s+", header=0, skiprows=5)
                pTot = spectraDf['#15_SpectFunctTot']
                k = spectraDf['#01_KStag']
                spl = CubicSpline(k, pTot)
                k1, k2 = k[0], k[len(k)-1]
                kList.append(fminbound(lambda x: -spl(x), k1, k2))
                # kList.append(k[pTot==max(pTot)])

        # print(np.percentile(kList, 84), np.percentile(kList, 50), np.percentile(kList, 16))
        # print(kList)
        k_eta_23 = np.percentile(kList, 50)
        k_eta_23_pos_er = np.percentile(kList, 84) - k_eta_23
        k_eta_23_neg_er = np.percentile(kList, 16) - k_eta_23
        fitDict["k_eta_23"] = (k_eta_23, k_eta_23_neg_er, k_eta_23_pos_er)
        print(f"{SOLVER_DICT[infoDict['solver']]}: {k_eta_23:.2f}^{{{k_eta_23_pos_er:.2f}}}_{{{k_eta_23_neg_er:.2f}}}")
        # ylim = plCurObj.gca().get_ylim()
        # plCurObj.plot([k_nu_23, k_nu_23], [ylim[0], 2e-10], color=colorDict[infoDict['solver']], scaley=False)

    elif spectType == "vels":
        if verbose: print("Kin Spectra-> Fitting for:", infoDict['solver'])
        kinFit = cfp.fit(Log10_P_kin, kFit, log10PTotFit, xerr=None, yerr=deltaLog10PTotFit, params=params, n_random_draws=1000)
        cfp.plot(fact * 10**Log10_P_kin(kFit, *(kinFit.popt)) * compensateFitFact, kFit, color="black")
        fitDict["A_kin"] = (kinFit.popt[0], kinFit.perr[0][0], kinFit.perr[0][1])
        fitDict["p_bn"] = (kinFit.popt[1], kinFit.perr[1][0], kinFit.perr[1][1])
        fitDict["k_bn"] = (kinFit.popt[2], kinFit.perr[2][0], kinFit.perr[2][1])
        fitDict["k_tilde_nu"] = (kinFit.popt[3], kinFit.perr[3][0], kinFit.perr[3][1])
        fitDict["p_nu"] = (kinFit.popt[4], kinFit.perr[4][0], kinFit.perr[4][1])
    
    return fitDict

def postPlot(plObj, spectType, compensated=False):
    ax = plObj.ax()
    ax.set_xscale('log')
    ax.set_yscale('log')
    if spectType == "mags":
        ax.set_ylabel(r'$P_\mathrm{mag}$')
        ax.get_xaxis().set_ticks([])
        # ax.set_xlabel(r'$k$')
        # ax.legend(loc='best')
    elif spectType == "vels":
        if compensated:
            ax.set_ylabel(r'$k^{1.7}P_\mathrm{kin}$')
        else:
            ax.set_ylabel(r'$P_\mathrm{kin}$')
        ax.get_xaxis().set_ticks([])
        ax.legend(loc='best')
    elif spectType == "cur":
        ax.set_ylabel(r'$P_\mathrm{cur}$')
        ax.set_xlabel(r'$k$')
        # ax.legend(loc='best')

def plotScaleLoc(plObj, solverFit, type):
    ax = plObj.ax()
    ylim = ax.get_ylim()
    if type == "mags":
        ax.plot([10, 10], [ylim[0], 2*ylim[0]], color="white", scaley=False)
        ax.text(10, 2.5*ylim[0], r"$k_\nu$", color="white")
        # maxKEta = 0
        # for solver in solverFit:
        #     val = solverFit[solver]
        #     ax.plot([val["k_tilde_eta"][0]**(1/val["p_eta"][0]), val["k_tilde_eta"][0]**(1/val["p_eta"][0])], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
        #     if val["k_tilde_eta"][0]**(1/val["p_eta"][0]) > maxKEta:
        #         maxKEta = val["k_tilde_eta"][0]**(1/val["p_eta"][0])
        # ax.text(maxKEta, 2.5*ylim[0], r"$k_\eta$", color="black")
    elif type == "vels":
        maxKNu = 0
        for solver in solverFit:
            val = solverFit[solver]
            k_nu = val["k_tilde_nu"][0]**(1/val["p_nu"][0])
            print(f"Value of k_nu for {solver}: {k_nu}")
            ax.plot([k_nu, k_nu], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
            if k_nu > maxKNu:
                maxKNu = k_nu
        ax.text(maxKNu, 2.5*ylim[0], r"$k_\nu$", color="black")
    elif type == "cur":
        maxKEta = 0
        for solver in solverFit:
            val = solverFit[solver]
            ax.plot([val["k_eta_23"][0], val["k_eta_23"][0]], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
            if val["k_eta_23"][0] > maxKEta:
                maxKEta = val["k_eta_23"][0]
        ax.text(maxKEta, 2.5*ylim[0], r"$k_\eta$", color="black")

def main(args):

    simList = getSolverSortedList(args.i)

    if args.kin_spect:
        for sim in simList:
            generateSpectra(sim, args.v, "vels", getInfoDict(sim), args.lf, args.uf, args.stf, args.n)
    if args.cur_spect:
        for sim in simList:
            generateSpectra(sim, args.v, "cur", getInfoDict(sim), args.lf, args.uf, args.stf, args.n)
    if args.mag_spect:
        for sim in simList:
            generateSpectra(sim, args.v, "mags", getInfoDict(sim), args.lf, args.uf, args.stf, args.n)

    if args.kin_plot:
        solverKinFit = {}
        for simDir in simList:
            infoDict = getInfoDict(simDir)
            if args.no_shift:
                fact = 1
            else:
                fact = FACT_DICT[infoDict['solver']]
            fit = True
            if not args.refit and os.path.exists(f"{args.o}/kinFitDict.pkl"):
                fit = False
            kinParams = {"A_kin": [0, 0.0015, np.inf], "p_bn": [0, 1, np.inf], "k_bn": [0.1, 4.0, 128], "k_tilde_nu": [0.1, 4.0, 128], "p_nu": [1, 1, 1+1e-6]}
            if os.path.exists(simDir+"/spectra/kinFitInit.txt"):
                kinParams = txtToCfpDict(simDir+"/spectra/kinFitInit.txt")
                print("Using kin dict:", kinParams)
            plKinObj, fitDict = plotSpectra(simDir, 1, "vels", fact, infoDict, kinParams, args.o, compensate=args.compensate, fit=fit)
            solverKinFit[infoDict['solver']] = fitDict
            postPlot(plKinObj, "vels", args.compensate)
        dumpDict(solverKinFit, f"{args.o}/kinFitDict.pkl")
        plotScaleLoc(plKinObj, solverKinFit, "vels")
        if args.compensate:
            plKinObj.ax().figure.savefig(f"{args.o}/Compensated Kinetic Spectra.pdf")
        else:
            plKinObj.ax().figure.savefig(f"{args.o}/Kinetic Spectra.pdf")
        plKinObj.ax().figure.clf(); plKinObj.ax().cla(); pl.close(); plKinObj = None

    
    if args.mag_plot:
        solverMagFit = {}
        for simDir in simList:
            infoDict = getInfoDict(simDir)
            if args.no_shift:
                fact = 1
            else:
                fact = FACT_DICT[infoDict['solver']]
            fit = True
            if not args.refit and os.path.exists(f"{args.o}/magFitDict.pkl"):
                fit = False
            # Mag fit has been disabled for now
            magParams = {"A_mag": [0, 0.0001, np.inf], "p_mag": [0, 1, np.inf], "p_eta": [0, 1, np.inf], "k_tilde_eta": [0, 4.0, 128]}
            plMagObj, fitDict = plotSpectra(simDir, 1, "mags", fact, infoDict, magParams, args.o, compensate=False, fit=fit)
            solverMagFit[infoDict['solver']] = fitDict
            postPlot(plMagObj, "mags")
        dumpDict(solverMagFit, f"{args.o}/magFitDict.pkl")
        plotScaleLoc(plMagObj, solverMagFit, "mags")
        plMagObj.ax().figure.savefig(f"{args.o}/Magnetic Spectra.pdf")
        plMagObj.ax().figure.clf(); plMagObj.ax().cla(); pl.close(); plMagObj = None

    if args.cur_plot:
        solverCurFit = {}
        for simDir in simList:
            infoDict = getInfoDict(simDir)
            if args.no_shift:
                fact = 1
            else:
                fact = FACT_DICT[infoDict['solver']]
            fit = True
            if not args.refit and os.path.exists(f"{args.o}/curFitDict.pkl"):
                fit = False
            plCurObj, fitDict = plotSpectra(simDir, 1, "cur", fact, infoDict, None, args.o, compensate=False, fit=fit)
            solverCurFit[infoDict['solver']] = fitDict
            postPlot(plCurObj, "cur")
        dumpDict(solverCurFit, f"{args.o}/curFitDict.pkl")
        plotScaleLoc(plCurObj, solverCurFit, "cur")
        plCurObj.ax().figure.savefig(f"{args.o}/Current Spectra.pdf")
        plCurObj.ax().figure.clf(); plCurObj.ax().cla(); pl.close(); plCurObj = None



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="./", help="Output Directory")
    parser.add_argument("-kin_spect", action="store_true", help="Generate kinetic spectra")
    parser.add_argument("-cur_spect", action="store_true", help="Generate current spectra")
    parser.add_argument("-mag_spect", action="store_true", help="Generate magnetic spectra")
    parser.add_argument("-kin_plot", action="store_true", help="Plot kinetic spectra")
    parser.add_argument("-compensate", action="store_true", help="Compensate for k^1.7. Works only for the kinetic spectra.")
    parser.add_argument("-cur_plot", action="store_true", help="Plot current spectra")
    parser.add_argument("-mag_plot", action="store_true", help="Plot magnetic spectra")
    parser.add_argument("-v", type=int, default=2, help="Verbose")
    parser.add_argument("-lf", type=float, default=1e-7, help="Lower bound for kinematic phase")
    parser.add_argument("-uf", type=float, default=1e-3, help="Upper bound for kinematic phase")
    parser.add_argument("-stf", type=float, default=5.0, help="Start time for kinematic phase")
    parser.add_argument("-n", type=int, default=1, help="Number of processors for spectra generation")
    parser.add_argument("-refit", action="store_true", help="Refit the spectra even if a fit file is present")
    parser.add_argument("-no_shift", action="store_true", help="Do not shift the spectra for different solvers")

    commonKeys = ["n", "o", "v", "mag_spect", "kin_spect", "cur_spect", "mag_plot", "kin_plot", "cur_plot", "legend", "lf", "uf", "stf", "refit", "compensate", "no_shift"]

    args = parseArgs(parser.parse_args(), commonKeys)

    main(args)