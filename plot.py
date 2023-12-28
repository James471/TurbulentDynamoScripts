#!/usr/bin/env python3
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.optimize import curve_fit
import argparse
import pickle
import os

E_MAG_AXIS_INDEX = 0
E_RATIO_AXIS_INDEX = 1
E_VRMS_AXIS_INDEX = 2

E_MAG_COLUMN_INDEX = 11
E_KIN_COLUMN_INDEX = 9
V_RMS_COLUMN_INDEX = 13
TIME_COLUMN_INDEX = 0

LOW_CUT = 10**10
HIGH_AMPLIFY = 10

FIT_COLOR_CYCLE = ["b", "g", "r", "c", "m", "y", "k"]


def getInfoDict(args):
    if os.path.exists(args.i + "/info.pkl"):
        with open(args.i + "/info.pkl", "rb") as f:
            return pickle.load(f)
    return None

def getFitLabel(args, alpha):
    infoDict = getInfoDict(args)
    if args.e is not None:
        return args.e+f", alpha={round(alpha, 2)}"
    if infoDict is None and args.e is None:
        return ("Unknown" if args.e is None else args.e)+f", alpha={round(alpha, 2)}"
    if infoDict is not None and args.e is None:
        if infoDict["solver"] == "bk-usm":
            return f"BK-USM, Cut={infoDict['mcut']}, alpha={round(alpha, 2)}"
        else:
            return f"{infoDict['solver']}, alpha={round(alpha, 2)}"
    

def getPlotLabel(args):
    infoDict = getInfoDict(args)
    if args.e is not None:
        return args.e
    if args.e is None and infoDict is None:
        return "Unknown"
    if args.e is None and infoDict is not None:
        if infoDict["solver"] == "bk-usm":
            return f"BK-USM, Cut={infoDict['mcut']}"
        else:
            return infoDict["solver"]


def addPlotInfo(args, fig, axes, isNewFig):
    # Title for old plots will not be changed
    if isNewFig:
        infoDict = getInfoDict(args)
        if args.title is not None:
            fig.suptitle(args.title)
        elif args.title is None and infoDict is None:
            fig.suptitle("Unknown")
        elif args.title is None:
            fig.suptitle(f"v={infoDict['v']}")


def getNewFigYLim(data):
    dataMin = np.min(np.ma.masked_invalid(data))
    dataMax = np.max(np.ma.masked_invalid(data))

    low = max(dataMin, dataMax / LOW_CUT)
    high = dataMax * HIGH_AMPLIFY

    return (low, high)


def getOldFigYLim(data, fig, ax):
    oldYLow = ax.get_ylim()[0]
    oldYHigh = ax.get_ylim()[1]

    dataMin = np.min(np.ma.masked_invalid(data))
    dataMax = np.max(np.ma.masked_invalid(data))

    high = max(oldYHigh, dataMax * HIGH_AMPLIFY)
    low = max(oldYLow, dataMin)

    return (low, high)


def adjustPlotAxis(args, fig, axes, isNewFig, f):
    if not args.no_adj_mag:
        ylim_mag = (
            args.ylim_mag
            if args.ylim_mag is not None
            else getNewFigYLim(getEMag(f))
            if isNewFig
            else getOldFigYLim(getEMag(f), fig, axes[E_MAG_AXIS_INDEX])
        )
        axes[E_MAG_AXIS_INDEX].set_ylim(ylim_mag[0], ylim_mag[1])
        print("Adjusted E_mag y-limits.")
    print(f"Set ylim_mag to {axes[E_MAG_AXIS_INDEX].get_ylim()}")

    if not args.no_adj_ratio:
        ylim_ratio = (
            args.ylim_ratio
            if args.ylim_ratio is not None
            else getNewFigYLim(getEMagOverEKin(f))
            if isNewFig
            else getOldFigYLim(getEMagOverEKin(f), fig, axes[E_RATIO_AXIS_INDEX])
        )
        axes[E_RATIO_AXIS_INDEX].set_ylim(ylim_ratio[0], ylim_ratio[1])
        print("Adjusted E_ratio y-limits.")
    print(f"Set ylim_ratio to {axes[E_RATIO_AXIS_INDEX].get_ylim()}")


def getTurnOverTime(args):
    infoDict = getInfoDict(args)
    if infoDict == None:
        return args.t
    return infoDict["turnover_time"]


def getNonDimensionalTime(f, args):
    turnOverTime = getTurnOverTime(args)
    return (f[TIME_COLUMN_INDEX] / turnOverTime).flatten()


def getEMag(f):
    return f[E_MAG_COLUMN_INDEX].flatten()


def getEKin(f):
    return f[E_KIN_COLUMN_INDEX].flatten()


def getVRMS(f):
    return f[V_RMS_COLUMN_INDEX].flatten()


def getEMagOverEKin(f):
    return (f[E_MAG_COLUMN_INDEX] / f[E_KIN_COLUMN_INDEX]).flatten()


def addPlot(args, fig, axes, isNewFig):
    f = np.loadtxt(args.i + "/Turb.dat", unpack=True)

    ax_emag = axes[E_MAG_AXIS_INDEX]
    ax_ratio = axes[E_RATIO_AXIS_INDEX]
    ax_vrms = axes[E_VRMS_AXIS_INDEX]

    if isNewFig:
        print("Addig new x and y labels to new figure")
        ax_emag.set_yscale("log")
        ax_ratio.set_yscale("log")

        ax_emag.set_xlabel("Turnover Time", fontsize=14)
        ax_emag.set_ylabel(r"$E_{mag}$", fontsize=14)

        ax_ratio.set_xlabel("Turnover Time", fontsize=14)
        ax_ratio.set_ylabel(r"$\frac{E_{mag}}{E_{kin}}$", fontsize=18)

        ax_vrms.set_xlabel("Turnover Time", fontsize=14)
        ax_vrms.set_ylabel(r"$v_{rms}$", fontsize=14)

        addPlotInfo(args, fig, axes, isNewFig)
    else:
        print("Reusing old labels")


    ax_emag.plot(
        getNonDimensionalTime(f, args),
        getEMag(f),
        label=getPlotLabel(args),
    )
    color = ax_emag.get_lines()[-1].get_color()
    ax_emag.legend()
    ax_ratio.plot(
        getNonDimensionalTime(f, args),
        getEMagOverEKin(f),
        label=getPlotLabel(args),
        color=color
    )
    ax_ratio.legend()
    ax_vrms.plot(
        getNonDimensionalTime(f, args),
        getVRMS(f),
        label=getPlotLabel(args),
        color=color
    )
    ax_vrms.legend()

    adjustPlotAxis(args, fig, axes, isNewFig, f)

def model(t, alpha, A, t0):
    return A*np.exp(alpha*(t-t0))

def getTransient(args):
    path = args.i + "/Turb.dat"
    f = np.loadtxt(path, unpack=True, skiprows=args.skiprows)
    ratio = getEMagOverEKin(f)
    return f[:, np.where(np.logical_and(ratio > args.fit_range[0], ratio < args.fit_range[1]))]

def getFitParams(args, f):
    x = getNonDimensionalTime(f, args)
    y = getEMagOverEKin(f)
    alphaGuess = np.log(y[-1]/y[0])/(x[-1]-x[0])
    popt, pcov = curve_fit(lambda t, alpha: model(t, alpha, y[0], x[0]), x, y, p0=[alphaGuess])
    return popt[0], pcov[0][0]

def addFit(args, fig, axes):
    f = getTransient(args)
    ax_ratio = axes[E_RATIO_AXIS_INDEX]
    alpha, alphaErr = getFitParams(args, f)
    ax_ratio.plot(getNonDimensionalTime(f, args), 
                  model(getNonDimensionalTime(f, args), alpha, getEMagOverEKin(f)[0], 
                        getNonDimensionalTime(f, args)[0]), label=getFitLabel(args, alpha))
    ax_ratio.legend()

def savePlot(args, fig):
    if args.save:
        if args.outdir is None:
            print("Saving the figure at", args.i + "/Turb.png")
            fig.savefig(args.i + "/Turb.png", dpi=250)
        else:
            print("Saving the figure at", args.outdir + "/Turb.png")
            fig.savefig(args.outdir + "/Turb.png", dpi=250)
    if args.show:
        pl.show()

def parseArgs(args):
    return args


def main(args, fig=None, axes=None):
    print("Input directory=", args.i)
    isNewFig = False

    if fig == None:
        fig = pl.figure(figsize=(12, 12))
        gs = gridspec.GridSpec(2, 4, figure=fig)
        gs.update(wspace=0.8)
        axes = [0, 0, 0]
        axes[E_MAG_AXIS_INDEX] = pl.subplot(gs[0, :2])
        axes[E_RATIO_AXIS_INDEX] = pl.subplot(gs[0, 2:])
        axes[E_VRMS_AXIS_INDEX] = pl.subplot(gs[1, 1:3])
        isNewFig = True

    addPlot(args, fig, axes, isNewFig)
    if args.fit:
        addFit(args, fig, axes)
    savePlot(args, fig)

    return fig, axes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to automate visualization of FLASH simulation"
    )
    parser.add_argument("-i", required=True, type=str, help="Input directory")
    parser.add_argument("-outdir", type=str, help="Output directory")
    parser.add_argument("-t", type=float, default=1, help="1 turn over time")
    parser.add_argument("-save", action="store_true", help="Save the figure")
    parser.add_argument(
        "-ylim_mag", type=float, nargs=2, help="Y-axis limits for E_mag plot"
    )
    parser.add_argument(
        "-ylim_ratio", type=float, nargs=2, help="Y-axis limits for E_mag/ERot plot"
    )
    parser.add_argument(
        "-no_adj_ratio", action="store_true", help="Don't adjust Ekin/Emag axis"
    )
    parser.add_argument(
        "-no_adj_mag", action="store_true", help="Don't adjust Emag axis"
    )
    parser.add_argument("-show", action="store_true", help="Show the figure")
    parser.add_argument("-e", type=str, help="Info to put in the legend")
    parser.add_argument("-title", type=str, help="Title of the plot")
    parser.add_argument("-fit", action="store_true", help="Fit the data to a power law")
    parser.add_argument("-fit_range", type=float, nargs=2, help="Range of data to fit")
    parser.add_argument("-skiprows", type=int, default=0, help="Number of rows to skip in the data file")

    args = parseArgs(parser.parse_args())
    main(args)
