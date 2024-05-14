#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

from constants import *
from datautils import *
from utils import *
from myconfig import *
import designParams

sys.path.append(PYTHON_PATH)
import cfp


def model(t, alpha, lnA, t0):
    return lnA + alpha * (t - t0)


def addPlot(axMach, axRatio, f, label, color, low, high, velocity, stf, fitMethod):
    t = getNonDimensionalTime(f, velocity)
    mach = getVRMS(f) / getCsRMS(f)
    ratio = getEMagOverEKin(f)
    axMach.plot(t, mach, color=color, label=label)
    axRatio.plot(t, ratio, color=color)
    
    fitMask = (t > stf) & (ratio > low) & (ratio < high)
    t_fit = t[fitMask]
    ratio_fit = ratio[fitMask]
    alphaGuess = np.log(ratio_fit[-1]/ratio_fit[0])/(t_fit[-1]-t_fit[0])
    lnAGuess = np.log(ratio_fit[0])
    guessParam = {"alpha": alphaGuess, "lnA": lnAGuess}

    fit = {}

    if fitMethod == "scp":
        popt, pcov = curve_fit(lambda t, alpha, lnA: model(t, alpha, lnA, t_fit[0]), t_fit, np.log(ratio_fit), p0=[alphaGuess, lnAGuess])
        fit["alpha"] = (popt[0], -np.sqrt(pcov[0][0]), np.sqrt(pcov[0][0]))
        fit["lnA"] = (popt[1], -np.sqrt(pcov[1][1]), np.sqrt(pcov[1][1]))
    elif fitMethod == "asym":
        temp = cfp.fit(lambda t, alpha, lnA: model(t, alpha, lnA, t_fit[0]), t_fit, np.log(ratio_fit), params=guessParam)
        fit["alpha"] = (temp.popt[0], temp.perr[0][0], temp.perr[0][1])
        fit["lnA"] = (temp.popt[1], temp.perr[1][0], temp.perr[1][1])
    elif fitMethod == "syst":
        temp = cfp.fit(lambda t, alpha, lnA: model(t, alpha, lnA, t_fit[0]), t_fit, np.log(ratio_fit), params=guessParam, perr_method="systematic", dat_frac_for_systematic_perr=0.2)
        fit["alpha"] = (temp.popt[0], temp.perr[0][0], temp.perr[0][1])
        fit["lnA"] = (temp.popt[1], temp.perr[1][0], temp.perr[1][1])
    
    print(f"{label}: {fit["alpha"]}")
    
    axRatio.plot(t_fit, np.exp(model(t_fit, fit["alpha"][0], fit["lnA"][0], t_fit[0])), color="black")

    return fit


def main(args):
    
    figMach = pl.figure()
    figRatio = pl.figure()
    axMach = figMach.add_subplot(111)
    axRatio = figRatio.add_subplot(111)
    axRatio.set_xlabel(r"$t / t_{\mathrm{turb}}$")
    axRatio.set_ylabel(r"$E_\mathrm{mag}/E_\mathrm{kin}$")
    axMach.set_ylabel(r"$\mathcal{M}$")

    axRatio.set_yscale("log")

    simList = getSolverSortedList(args.i)

    fitDict = {}

    for index, file in enumerate(simList):
        infoDict = getInfoDict(file)
        velocity = infoDict["v"]
        label = SOLVER_DICT[infoDict["solver"]]
        color = COLOR_DICT[infoDict["solver"]]
        fit = addPlot(axMach, axRatio, loadFile(file + "/Turb.dat", args.sr[index]), label, color, args.lf[index], args.uf[index], velocity, args.stf[index], args.fit)
        fitDict[label] = fit

    dumpDict(fitDict, args.o+"/fitDict.pkl")

    order = [i for i in range(0, 6)]
    order[1], order[2], order[3], order[4] = order[3], order[1], order[4], order[2]
    handles, labels = axMach.get_legend_handles_labels()
    axMach.legend([handles[i] for i in order], [labels[i] for i in order], ncol=3)
    if args.ld is not None:
        axRatio.set_ylim(bottom=args.ld)
    if args.ud is not None:
        axRatio.set_ylim(top=args.ud)
    figMach.savefig(args.o+"/MachGrowth.pdf")
    figRatio.savefig(args.o+"/RatioGrowth.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="./", help="Output Directory")
    parser.add_argument("-lf", type=float, nargs="*", default=[1e-7], help="Lower bound for fit")
    parser.add_argument("-uf", type=float, nargs="*", default=[1e-3], help="Upper bound for fit")
    parser.add_argument("-sr", type=int, nargs="*", default=[20000], help="Skip rows")
    parser.add_argument("-stf", type=float, nargs="*", default=[5.0], help="Skip turnover time before fit")
    parser.add_argument("-ld", type=float, help="Low bound to display")
    parser.add_argument("-ud", type=float, help="Upper bound to display")
    parser.add_argument("-fit", type=str, default="syst", help="Fit method to use")

    commonKeys = ["ld", "ud", "fit"]

    args = parseArgs(parser.parse_args(), commonKeys)
    print(args)

    main(args)