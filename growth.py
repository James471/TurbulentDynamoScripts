#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

from constants import *
from datautils import *
from utils import *
from myconfig import *
import designParams

sys.path.append(PYTHON_PATH)
import cfpack as cfp


def model(t, alpha, lnA, t0):
    return lnA + alpha * (t - t0)


def addPlot(axMach, axRatio, f, label, color, low, high, velocity, stf, sbs, fitMethod, fitParam=None, nofit=False):
    '''
    fitParam is valid only if fitMethod="none"
    '''
    t = getNonDimensionalTime(f, velocity)
    mach = getVRMS(f) / getCsRMS(f)
    ratio = getEMagOverEKin(f)
    axMach.plot(t, mach, color=color)
    axRatio.plot(t, ratio, color=color, label=label)

    fit = None
    
    if not nofit:
        fitMask = (t > stf) & (ratio > low) & (ratio < high)
        t_fit = t[fitMask]
        ratio_fit = ratio[fitMask]
        alphaGuess = np.log(ratio_fit[-1]/ratio_fit[0])/(t_fit[-1]-t_fit[0])
        lnAGuess = np.log(ratio_fit[0])
        guessParam = {"alpha": [0, alphaGuess, np.inf], "lnA": [-np.inf, lnAGuess, np.inf]}

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
        elif fitMethod == "none":
            fit = fitParam

        if fitMethod != "none":
            sat_level = np.mean(ratio[t>sbs])
            sat_er = np.std(ratio[t>sbs])
            fit["sat"] = (sat_level, -sat_er, sat_er)

        print(f"{label}_tau: {fit['alpha']}")
        print(f"{label}_sat: {fit['sat']}")
    
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
    axMach.axes.xaxis.set_ticklabels([])

    axRatio.set_yscale("log")

    simList = getSolverSortedList(args.i)

    fitDict = {}

    if not args.refit and os.path.exists(args.o+"/growthDict.pkl"):
        fitDict = loadDict(args.o+"/growthDict.pkl")
        fitMethod = "none"
    else:
        fitMethod = args.fit
        fitDict = {solver: {} for solver in SOLVER_DICT.keys()}

    for index, file in enumerate(simList):
        infoDict = getInfoDict(file)
        velocity = infoDict["v"]
        label = SOLVER_DICT[infoDict["solver"]]
        color = COLOR_DICT[infoDict["solver"]]
        print(f"Processing {label}")
        n = getNForTurbDat(file, res=5e-3)
        name = "/Turb.dat_cleaned" if os.path.exists(file + "/Turb.dat_cleaned") else "/Turb.dat"
        fit = addPlot(axMach, axRatio, loadFile(file + name, args.sr[index], n), label, color, args.lf[index], 
                      args.uf[index], velocity, args.stf[index], args.sbs[index], fitMethod, 
                      fitParam=fitDict[infoDict["solver"]], nofit=args.nofit)
        if not args.nofit:
            fitDict[infoDict["solver"]] = fit

    if not args.nofit:
        dumpDict(fitDict, args.o+"/growthDict.pkl")

    # order = [i for i in range(0, 6)]
    # order[1], order[2], order[3], order[4] = order[3], order[1], order[4], order[2]
    # handles, labels = axMach.get_legend_handles_labels()
    # axMach.legend([handles[i] for i in order], [labels[i] for i in order], ncol=3, loc="best")
    axRatio.legend(loc="best")
    # axMach.legend(ncol=2, loc="best")
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
    parser.add_argument("-sr", type=int, nargs="*", default=[12000], help="Skip rows")
    parser.add_argument("-stf", type=float, nargs="*", default=[5.0], help="Skip turnover time before fit")
    parser.add_argument("-sbs", type=float, nargs="*", default=[60.0], help="Skip turnover time before measuring saturation")
    parser.add_argument("-ld", type=float, nargs="?", default=1e-8, help="Low bound to display")
    parser.add_argument("-ud", type=float, default=5e0, help="Upper bound to display")
    parser.add_argument("-fit", type=str, default="syst", help="Fit method to use")
    parser.add_argument("-refit", action="store_true", help="Refit the data")
    parser.add_argument("-nofit", action="store_true", help="Do not fit the data")

    commonKeys = ["ld", "ud", "fit", "refit", "nofit"]

    args = parseArgs(parser.parse_args(), commonKeys)
    print(args)

    main(args)