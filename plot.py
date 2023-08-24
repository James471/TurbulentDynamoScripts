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


def getInfoDict(args):
    if os.path.exists(args.i + "/info.pkl"):
        with open(args.i + "/info.pkl", "rb") as f:
            return pickle.load(f)
    return None


def getSolverInfo(args):
    infoDict = getInfoDict(args)
    if infoDict == None:
        return "Unknown" if args.e is None else args.e
    if infoDict["solver"] == "bk-usm":
        return f"BK-USM, Cut={infoDict['mcut']}"
    else:
        return f"{infoDict['solver']}"


def addPlotInfo(args, fig, axes, isNewFig):
    infoDict = getInfoDict(args)
    if infoDict == None:
        fig.suptitle("Unknown" if args.e is None else args.e)
        return
    if isNewFig:
        fig.suptitle(f"v={infoDict['v']}")
    if args.e is not None:
        fig.suptitle(fig.suptitle().get_text() + "\n" + args.e)


def getNewFigYLim(data):
    dataMin = np.min(data)
    dataMax = np.max(data)

    low = max(dataMin, dataMax / (10**10))
    high = dataMax * 10

    print(dataMax)

    return (low, high)


def getOldFigYLim(data, fig, ax):
    oldYLow = ax.get_ylim()[0]
    oldYHigh = ax.get_ylim()[1]

    dataMin = np.min(data)
    dataMax = np.max(data)

    high = max(oldYHigh, dataMax * 10)
    low = max(oldYLow, dataMin)

    return (low, high)


def adjustPlotAxis(args, fig, axes, isNewFig, f):
    if not args.no_adj_mag:
        ylim_mag = (
            args.ylim_mag
            if args.ylim_mag is not None
            else getNewFigYLim(f[E_MAG_COLUMN_INDEX])
            if isNewFig
            else getOldFigYLim(f[E_MAG_COLUMN_INDEX], fig, axes[E_MAG_AXIS_INDEX])
        )
        print(f"Setting ylim_mag to {ylim_mag}")
        axes[E_MAG_AXIS_INDEX].set_ylim(ylim_mag[0], ylim_mag[1])

    if not args.no_adj_ratio:
        ylim_ratio = (
            args.ylim_ratio
            if args.ylim_ratio is not None
            else getNewFigYLim(f[E_KIN_COLUMN_INDEX] / [E_MAG_COLUMN_INDEX])
            if isNewFig
            else getOldFigYLim(
                f[E_KIN_COLUMN_INDEX] / [E_MAG_COLUMN_INDEX], fig, axes[E_RATIO_AXIS_INDEX]
            )
        )
        print(f"Setting ylim_ratio to {ylim_ratio}")
        axes[E_RATIO_AXIS_INDEX].set_ylim(ylim_ratio[0], ylim_ratio[1])


def getTurnOverTime(args):
    infoDict = getInfoDict(args)
    if infoDict == None:
        return args.t
    return infoDict["turn_overtime"]


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

    turnOverTime = getTurnOverTime(args)

    ax_emag.plot(
        f[TIME_COLUMN_INDEX] / turnOverTime,
        f[E_MAG_COLUMN_INDEX],
        label=getSolverInfo(args),
    )
    ax_emag.legend()
    ax_ratio.plot(
        f[TIME_COLUMN_INDEX] / turnOverTime,
        f[E_KIN_COLUMN_INDEX] / f[E_MAG_COLUMN_INDEX],
        label=getSolverInfo(args),
    )
    ax_ratio.legend()
    ax_vrms.plot(
        f[TIME_COLUMN_INDEX] / turnOverTime,
        f[V_RMS_COLUMN_INDEX],
        label=getSolverInfo(args),
    )
    ax_vrms.legend()

    adjustPlotAxis(args, fig, axes, isNewFig, f)
    pl.show()


def parseArgs(args):
    return args


def main(args, fig=None, axes=None):
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

    return fig, axes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to automate visualization of FLASH simulation"
    )
    parser.add_argument("-i", required=True, type=str, help="Input directory")
    parser.add_argument("-t", type=float, help="1 turn over time")
    parser.add_argument(
        "-ylim_mag", type=float, nargs=2, help="Y-axis limits for E_mag plot"
    )
    parser.add_argument(
        "-ylim_ratio", type=float, nargs=2, help="Y-axis limits for E_mag/ERot plot"
    )
    parser.add_argument("-no_adj_ratio", action="store_true", help="Don't adjust Ekin/Emag axis")
    parser.add_argument("-no_adj_mag", action="store_true", help="Don't adjust Emag axis")
    parser.add_argument("-e", type=str, help="Extra info to put in the title")

    args = parseArgs(parser.parse_args())
    main(args)
