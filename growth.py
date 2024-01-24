#!/usr/bin/env python3

import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.optimize import curve_fit
import argparse
import pickle
import os

E_MAG_COLUMN_INDEX = 11
E_KIN_COLUMN_INDEX = 9
V_RMS_COLUMN_INDEX = 13
Cs_RMS_COLUMN_INDEX = 14
TIME_COLUMN_INDEX = 0

def getTurnOverTime(v):
    return 1/(2*v)


def getNonDimensionalTime(f, v):
    turnOverTime = getTurnOverTime(v)
    return (f[TIME_COLUMN_INDEX] / turnOverTime).flatten()


def getEMag(f):
    return f[E_MAG_COLUMN_INDEX].flatten()


def getEKin(f):
    return f[E_KIN_COLUMN_INDEX].flatten()


def getVRMS(f):
    return f[V_RMS_COLUMN_INDEX].flatten()


def getCsRMS(f):
    return f[Cs_RMS_COLUMN_INDEX].flatten()


def getEMagOverEKin(f):
    return (getEMag(f) / getEKin(f)).flatten()


def loadFile(path, shift):
    return np.loadtxt(path, unpack=True, skiprows=shift)


def getInfoDict(dirPath):
    # if os.path.exists(args.i + "/info.pkl"):
    with open(dirPath + "/info.pkl", "rb") as f:
        return pickle.load(f)
    # return None


def model(t, alpha, lnA, t0):
    return lnA + alpha * (t - t0)


def addPlot(fig, axMach, axRatio, f, label, color, low, high, velocity, stf):
    t = getNonDimensionalTime(f, velocity)
    mach = getVRMS(f) / getCsRMS(f)
    ratio = getEMagOverEKin(f)
    axMach.plot(t, mach, color=color)
    axRatio.plot(t, ratio, label=label, color=color)
    
    fitMask = (t > stf) & (ratio > low) & (ratio < high)
    t_fit = t[fitMask]
    ratio_fit = ratio[fitMask]
    alphaGuess = np.log(ratio_fit[-1]/ratio_fit[0])/(t_fit[-1]-t_fit[0])
    lnAGuess = np.log(ratio_fit[0])

    popt, pcov = curve_fit(lambda t, alpha, lnA: model(t, alpha, lnA, t_fit[0]), t_fit, np.log(ratio_fit), p0=[alphaGuess, lnAGuess])
    print(f"Label:{label}, popt:{popt[0], np.exp(popt[1])}, err:{np.sqrt(pcov[0][0])}, {np.exp(popt[1]) * np.sqrt(pcov[1][1])}")
    axRatio.plot(t_fit, np.exp(model(t_fit, *(list(popt)+[t_fit[0]]))), color="black")


def parseArgs(args):
    numSim = -1
    for key, value in vars(args).items():
        if key not in commonKeys:
            if value is not None and len(value) > 1:
                numSim = len(value)
                break
            
    for key, value in vars(args).items():
        if key not in commonKeys:
            if value is None:
                setattr(args, key, [None for i in range(numSim)])
            elif len(value) == 1:
                setattr(args, key, [value[0] for i in range(numSim)])
    return args


def main(args):

    colorDict = {"8wave": "#377eb8", "HLLD": "#984ea3", "HLLC": "#4daf4a", 
                 "Roe": "#a65628", "bouchut-split": "#f781bf", "bk-usm": "#ff7f00"}
    labelDict = {"8wave": "Split-Roe", "bouchut-split": "Split-Bouchut", "Roe": "USM-Roe",
                 "bk-usm": "USM-BK", "HLLC": "USM-HLLC", "HLLD": "USM-HLLD"}
    orderDict = {"8wave": 0, "bouchut-split": 1, "Roe": 2, "HLLD": 3, "HLLC": 4, "bk-usm": 5}
    
    fig, axes = pl.subplots(nrows=2, sharex=True, gridspec_kw={'hspace': 0.05}, figsize=(8, 8))
    axMach, axRatio = axes
    axRatio.set_xlabel(r"$t / t_{\mathrm{turb}}$")
    axRatio.set_ylabel(r"$E_\mathrm{mag}/E_\mathrm{kin}$")
    axMach.set_ylabel(r"$\mathcal{M}$")
    axRatio.tick_params(axis="x",direction="in")
    axRatio.tick_params(axis="y",direction="in")
    axMach.tick_params(axis="x",direction="in")
    axMach.tick_params(axis="y",direction="in")
    axRatio.set_aspect('auto')
    axMach.set_aspect('auto')

    axRatio.set_yscale("log")

    order = [0 for i in range(len(args.i))]

    fileList = args.i
    for index, file in enumerate(fileList):
        infoDict = getInfoDict(file)
        velocity = infoDict["v"]
        label = labelDict[infoDict["solver"]]
        color = colorDict[infoDict["solver"]]
        addPlot(fig, axMach, axRatio, loadFile(file + "/Turb.dat", args.sr[index]), label, color, args.lf[index], args.uf[index], velocity, args.stf[index])
        order[index] = orderDict[infoDict["solver"]]
    order = np.argsort(order)
    order[1], order[2], order[3], order[4] = order[3], order[1], order[4], order[2]
    handles, labels = axRatio.get_legend_handles_labels()
    axRatio.legend([handles[i] for i in order], [labels[i] for i in order], ncol=3)
    if args.ld is not None:
        axRatio.set_ylim(bottom=args.ld)
    fig.savefig(args.o)


commonKeys = ["ld"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="Turb.pdf", help="Output Directory")
    parser.add_argument("-lf", type=float, nargs="*", help="Lower bound for fit")
    parser.add_argument("-uf", type=float, nargs="*", help="Upper bound for fit")
    parser.add_argument("-sr", type=int, nargs="*", help="Skip rows")
    parser.add_argument("-stf", type=float, nargs="*", default=[1.0], help="Skip turnover time before fit")
    parser.add_argument("-ld", type=float, help="Low bound to display")

    args = parseArgs(parser.parse_args())
    print(args)

    main(args)