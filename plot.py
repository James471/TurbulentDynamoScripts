#!/usr/bin/env python3
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import numpy as np
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


def getInfoDict(args):
    if os.path.exists(args.i + "/info.pkl"):
        with open(args.i + "/info.pkl", "rb") as f:
            return pickle.load(f)
    return None


def getPlotLabel(args):
    infoDict = getInfoDict(args)
    if infoDict == None:
        return "Unknown" if args.e is None else args.e
    elif infoDict["solver"] == "bk-usm":
        return f"BK-USM, Cut={infoDict['mcut']}"
    else:
        return f"{infoDict['solver']}"


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
    return f[TIME_COLUMN_INDEX] / turnOverTime


def getEMag(f):
    return f[E_MAG_COLUMN_INDEX]


def getEKin(f):
    return f[E_KIN_COLUMN_INDEX]


def getVRMS(f):
    return f[V_RMS_COLUMN_INDEX]


def getEMagOverEKin(f):
    return f[E_MAG_COLUMN_INDEX] / f[E_KIN_COLUMN_INDEX]


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
    ax_emag.legend()
    ax_ratio.plot(
        getNonDimensionalTime(f, args),
        getEMagOverEKin(f),
        label=getPlotLabel(args),
    )
    ax_ratio.legend()
    ax_vrms.plot(
        getNonDimensionalTime(f, args),
        getVRMS(f),
        label=getPlotLabel(args),
    )
    ax_vrms.legend()

    adjustPlotAxis(args, fig, axes, isNewFig, f)


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
    savePlot(args, fig)

    return fig, axes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to automate visualization of FLASH simulation"
    )
    parser.add_argument("-i", required=True, type=str, help="Input directory")
    parser.add_argument("-outdir", type=str, help="Output directory")
    parser.add_argument("-t", type=float, deafult=1, help="1 turn over time")
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
    parser.add_argument("-e", type=str, help="Extra info to put in the legend")
    parser.add_argument("-title", type=str, help="Title of the plot")

    args = parseArgs(parser.parse_args())
    main(args)
