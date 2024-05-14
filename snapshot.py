#!/usr/bin/env python
import argparse
import numpy as np
import textwrap
import pickle
import subprocess
import os
import re
import sys
sys.path.append("/home/100/jw5893/flash-tools/python")
import flashlib


E_MAG_COLUMN_INDEX = 11
E_KIN_COLUMN_INDEX = 9
V_RMS_COLUMN_INDEX = 13
CS_RMS_COLUMN_INDEX = 14
TIME_COLUMN_INDEX = 0

SOLVER_DICT = {"8wave": "Split-Roe", "bouchut-split": "Split-Bouchut", "Roe": "USM-Roe", 
                "HLLD": "USM-HLLD", "HLLC": "USM-HLLC", "bk-usm": "USM-BK"}

X_LABEL_DICT = {"Split-Roe": "", "Split-Bouchut": "", "USM-Roe": "",
                "USM-HLLD": "\$x/L\$", "USM-HLLC": "\$x/L\$", "USM-BK": "\$x/L\$"}

Y_LABEL_DICT = {"Split-Roe": "\$y/L\$", "Split-Bouchut": "", "USM-Roe": "",
                "USM-HLLD": "\$y/L\$", "USM-HLLC": "", "USM-BK": ""}

X_TICK_DICT = {"Split-Roe": "\"\"", "Split-Bouchut": "\"\"", "USM-Roe": "\"\"",
                "USM-HLLD": None, "USM-HLLC": None, "USM-BK": None}

Y_TICK_DICT = {"Split-Roe": None, "Split-Bouchut": "\"\"", "USM-Roe": "\"\"",
                "USM-HLLD": None, "USM-HLLC": "\"\"", "USM-BK": "\"\""}


def loadFile(path, shift=0):
    return np.loadtxt(path, unpack=True, skiprows=shift)


def getInfoDict(dirPath):
    filePath = dirPath + "/info.pkl"
    if os.path.exists(filePath):
        with open(filePath, "rb") as f:
            return pickle.load(f)
    else:
        raise FileNotFoundError(f"No info.pkl file found in {dirPath}.")


def getTurnOverTime(v):
    return 1/(2*v)


def getNonDimensionalTime(f, v):
    turnOverTime = getTurnOverTime(v)
    return (f[TIME_COLUMN_INDEX] / turnOverTime).flatten()


def getEMag(f):
    return f[E_MAG_COLUMN_INDEX].flatten()


def getEKin(f):
    return f[E_KIN_COLUMN_INDEX].flatten()


def getEMagOverEKin(f):
    return (getEMag(f) / getEKin(f)).flatten()


def getPlotFileNumber(turbFile, velocity, r, stf, lf, dt):
    t = getNonDimensionalTime(turbFile, velocity)
    ratio = getEMagOverEKin(turbFile)
    transientMask = (t > stf) & (ratio > lf)
    location = np.argmin(np.abs(ratio[transientMask] - r))
    snapshotTime = t[transientMask][location]
    plotFileNumber = int(round(snapshotTime / dt))

    print(f"Snapshot Time: {snapshotTime}t_turb")

    return plotFileNumber


def makePlots(simDirList, type, stream_var, r, stf, lf, outdir, cbar_type, flashplotlib_path):
    
    plotFileDict = {}
    boundDict = {}
    turbFileDict = {}

    if type == 'mag':
        datasetName = 'emag'
    elif type == 'kin':
        datasetName = 'ekin'

    for index, simDir in enumerate(simDirList):
        
        infoDict = getInfoDict(simDir)
        velocity = infoDict["v"]
        solver = SOLVER_DICT[infoDict["solver"]]
        dt = infoDict["dt"]

        print(f"Processing for {solver}")

        print("Loading Turb.dat")
        turbFile = loadFile(simDir + "/Turb.dat", 10)
        print("Loaded Turb.dat")

        plotFileNumber = getPlotFileNumber(turbFile, velocity, 1e-6, 5, 1e-7, dt)
        plotFile = simDir + "/Turb_hdf5_plt_cnt_" + str(plotFileNumber).zfill(4)
        plotFileDict[solver] = plotFile

        plotFileObj = flashlib.FlashGG(plotFile)
        snapshotTime = plotFileObj.scalars['time']

        time = turbFile[TIME_COLUMN_INDEX]
        e = turbFile[E_MAG_COLUMN_INDEX]
        eTot = e[time==snapshotTime][0]

        cmd = f'python3 {flashplotlib_path} -i {plotFile} -d {datasetName} -verbose 2 -outtype pdf -outdir {outdir} -outname {solver}-{type} -direction z -ncpu 1 -colorbar 0 -data_transform q/{eTot}'
        output = subprocess.check_output(cmd, shell=True).decode('utf-8')
        searchString = "flashplotlib.py: prep_map: min, max of data.*"
        iterms = re.findall(searchString, output)
        boundString = iterms[0]
        boundString = boundString.replace("flashplotlib.py: prep_map: min, max of data = ", "")
        minBound, maxBound = boundString.split(",")
        minBound, maxBound = float(minBound), float(maxBound)
        boundDict[solver] = (minBound, maxBound)
        turbFileDict[solver] = turbFile

    print(f"{datasetName} bounds:", boundDict)

    vMin = min([boundDict[solver][0] for solver in boundDict.keys()])
    vMax = max([boundDict[solver][1] for solver in boundDict.keys()])
    print(f"vMin: {vMin}, vMax: {vMax}")

    streamStr = ""
    if stream_var is not None:
        streamStr = f" -stream -stream_var {stream_var} -stream_color white"

    for index, simDir in enumerate(simDirList):
        infoDict = getInfoDict(simDir)
        velocity = infoDict["v"]
        solver = SOLVER_DICT[infoDict["solver"]]
        dt = infoDict["dt"]

        print(f"Plotting for {solver}")

        plotFile = plotFileDict[solver]
        plotFileObj = flashlib.FlashGG(plotFile)
        snapshotTime = plotFileObj.scalars['time']
        turbFile = turbFileDict[solver]
        
        time = turbFile[TIME_COLUMN_INDEX]
        e = turbFile[E_MAG_COLUMN_INDEX]
        eTot = e[time==snapshotTime][0]
        
        plotFile = plotFileDict[solver]
        cmd = f'python3 {flashplotlib_path} -i {plotFile} -d {datasetName} -verbose 2 -outtype pdf -outdir {outdir} -outname {solver}-{type} -direction z -ncpu 1 -colorbar 0 -vmin {vMin} -vmax {vMax} -axes_label "{X_LABEL_DICT[solver]}" "{Y_LABEL_DICT[solver]}" "" -axes_unit "" "" "" -plotlabel {solver} -time_unit "\$t_\mathrm{{turb}}\$" -time_scale {getTurnOverTime(velocity)} -axes_format {X_TICK_DICT[solver]} {Y_TICK_DICT[solver]} -data_transform q/{eTot} -cmap {cbar_type} {streamStr}'
        os.system(cmd)

    print("Making colorbars")
    cbarCmd = f'python3 {flashplotlib_path} -i {simDirList[0]}/Turb_hdf5_plt_cnt_0000 -d dens -verbose 2 -outtype pdf -outdir {outdir} -outname {type}Bar -ncpu 1 -colorbar only -vmin {vMin} -vmax {vMax} -cmap_label "\$E_\mathrm{{{type}}}/E_\mathrm{{{type};total}}\$ (Projection length \$z=L\$)" -doubleHeight 1 -cmap {cbar_type}'
    os.system(cbarCmd)


def makeTable(simList, outdir, type):
    texScript = f"""\\documentclass{{article}}
        \\usepackage{{graphicx}}
        \\usepackage{{multirow}}
        \\usepackage{{geometry}}

        \\begin{{document}}

        \\newgeometry{{left=1cm, right=1cm}}
        \\begin{{table}}[h!]
            \\centering
            \\setlength\\tabcolsep{{1.5pt}}
            \\begin{{tabular}}{{cccc}}
                \\includegraphics[width=0.317\\linewidth]{{{outdir}/Split-Roe-{type}.pdf}} & \\includegraphics[width=0.28\\linewidth]{{{outdir}/Split-Bouchut-{type}.pdf}} &
                \\includegraphics[width=0.28\\linewidth]{{{outdir}/USM-Roe-{type}.pdf}} & \\multirow{{2}}{{*}}[5.22cm]{{\\includegraphics[height=11.1cm]{{{outdir}/{type}Bar.pdf}}}}\\\\
                \\includegraphics[width=0.317\\linewidth]{{{outdir}/USM-HLLD-{type}.pdf}} &
                \\includegraphics[width=0.28\\linewidth]{{{outdir}/USM-HLLC-{type}.pdf}} &
                \\includegraphics[width=0.28\\linewidth]{{{outdir}/USM-BK-{type}.pdf}} &
            \\end{{tabular}}
        \\end{{table}}

        \\end{{document}}"""
    
    with open(f"{type}.tex", "w") as f:
        f.write(textwrap.dedent(texScript))
    os.system(f"pdflatex -interaction=nonstopmode --halt-on-error -output-directory {args.o} {args.o}/{type}.tex")
    os.system(f"rm {outdir}/{type}.aux {outdir}/{type}.log {outdir}/{type}.tex")



def main(args):

    simList = [sim for sim in args.i if os.path.isdir(sim)]
    outdir = args.o

    if args.kin:
        makePlots(simList, "kin", args.kin_stream, args.r, args.stf, args.lf, outdir, args.cbar_type, args.flashplotlib_path)
        if args.table:
            makeTable(simList, outdir, "kin")
    if args.mag:
        makePlots(simList, "mag", args.mag_stream, args.r, args.stf, args.lf, outdir, args.cbar_type, args.flashplotlib_path)
        if args.table:
            makeTable(simList, outdir, "mag")
    

def parseArgs(args):
    numSim = 1
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


commonKeys = ["r", "flashplotlib_path", "o", "kin", "mag", "cbar_type", "kin_stream", "mag_stream", "table"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="./", help="Output Directory")
    parser.add_argument("-r", type=float, default=1e-5, help="Energy ratio for projections")
    parser.add_argument("-stf", type=float, nargs="*", help="Start time for looking for ratio")
    parser.add_argument("-lf", type=float, nargs="*", help="Lower bound for looking for ratio")
    parser.add_argument("-kin", action='store_true', help="Plot kinetic energy")
    parser.add_argument("-mag", action='store_true', help="Plot magnetic energy")
    parser.add_argument("-table", action='store_true', help="Want the plots to be arranged in a table")
    parser.add_argument("-cbar_type", type=str, default="magma", help="Colorbar type")
    parser.add_argument("-kin_stream", type=str, default=None, choices=['vel', 'mag'], help="Streamline on kinetic energy")
    parser.add_argument("-mag_stream", type=str, default=None, choices=['vel', 'mag'], help="Streamline on magnetic energy")
    parser.add_argument("-flashplotlib_path", type=str, default="/home/100/jw5893/flash-tools/python/flashplotlib.py", help="Path to flashplotlib")

    args = parseArgs(parser.parse_args())
    print(args)

    main(args)




######################################################################
                                #OLD CODE
# import matplotlib.pyplot as pl
# import matplotlib.image as mpimg
# from matplotlib.gridspec import GridSpec


# def getGrid():
#     fig = pl.figure(figsize=(21.0, 10.0))
#     gs1 = GridSpec(2, 3, figure=fig, width_ratios=[6.5, 6.5, 6.5], height_ratios=[1, 1])
#     gs2 = GridSpec(1, 1, figure=fig)

#     # Add images to the grid
#     axes = [
#         fig.add_subplot(gs1[0, 0]),
#         fig.add_subplot(gs1[0, 1]),
#         fig.add_subplot(gs1[0, 2]),
#         fig.add_subplot(gs1[1, 0]),
#         fig.add_subplot(gs1[1, 1]),
#         fig.add_subplot(gs1[1, 2]),
#         fig.add_subplot(gs2[0, 0])
#     ]

#     gs1.update(right=0.8, wspace=-0.4, hspace=0.03)
#     gs2.update(left=0.62)

#     return fig, axes

# def addImage(ax, dirPath, solver, xlabel, ylabel, xtick, ytick, type):
#     ax.imshow(mpimg.imread(f"{dirPath}/{solver}-{type}.pdf"), extent=[-0.52, 0.52, -0.52, 0.52])
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_xticks(xtick)
#     ax.set_yticks(ytick)
#     ax.set_xticklabels(xtick)
#     ax.set_yticklabels(ytick)
#     ax.set_frame_on(False)
#     ax.tick_params(axis="x", which='both', bottom=False, top=False)
#     ax.tick_params(axis="y", which='both', left=False, right=False)

# def makeCBar(ax, dirPath, type):
#     ax.imshow(mpimg.imread(f"{dirPath}/{type}Bar.pdf"))
#     ax.set_frame_on(False)
#     ax.tick_params(axis="x", which='both', bottom=False, top=False)
#     ax.tick_params(axis="y", which='both', left=False, right=False)
#     ax.set_xticks([])
#     ax.set_yticks([])
#     ax.set_xticklabels([])
#     ax.set_yticklabels([])


# def main():



    # Do data stuff here



    # kinFig, kinAxes = getGrid()
    # magFig, magAxes = getGrid()


    # xtick = [-0.4, 0, 0.4]
    # ytick = [-0.4, 0, 0.4]
    # xlabel = r"x/L"
    # ylabel = r"y/L"
    # kinAxesDict = {"Split-Roe": kinAxes[0], "Split-Bouchut": kinAxes[1], "USM-Roe": kinAxes[2],
    #             "USM-HLLD": kinAxes[3], "USM-HLLC": kinAxes[4], "USM-BK": kinAxes[5], "cBar": kinAxes[6]}
    # magAxesDict = {"Split-Roe": magAxes[0], "Split-Bouchut": magAxes[1], "USM-Roe": magAxes[2],
    #             "USM-HLLD": magAxes[3], "USM-HLLC": magAxes[4], "USM-BK": magAxes[5], "cBar": magAxes[6]}
    # xtickDict = {"Split-Roe": [], "Split-Bouchut": [], "USM-Roe": [],
    #             "USM-HLLD": xtick, "USM-HLLC": xtick, "USM-BK": xtick}
    # ytickDict = {"Split-Roe": ytick, "Split-Bouchut": [], "USM-Roe": [],
    #             "USM-HLLD": ytick, "USM-HLLC": [], "USM-BK": []}
    # xlabelDict = {"Split-Roe": "", "Split-Bouchut": "", "USM-Roe": "",
    #                 "USM-HLLD": xlabel, "USM-HLLC": xlabel, "USM-BK": xlabel}
    # ylabelDict = {"Split-Roe": ylabel, "Split-Bouchut": "", "USM-Roe": "",
    #                 "USM-HLLD": ylabel, "USM-HLLC": "", "USM-BK": ""}

    # for index, sim in enumerate(simList):
    #     infoDict = getInfoDict(sim)
    #     solver = solverDict[infoDict["solver"]]
    #     addImage(kinAxesDict[solver], args.o, solver, xlabelDict[solver], ylabelDict[solver], xtickDict[solver], ytickDict[solver], "kinetic")
    #     addImage(magAxesDict[solver], args.o, solver, xlabelDict[solver], ylabelDict[solver], xtickDict[solver], ytickDict[solver], "magnetic")

    # makeCBar(kinAxesDict["cBar"], args.o, "kin")
    # makeCBar(magAxesDict["cBar"], args.o, "mag")

    # kinFig.savefig(f"{args.o}/kinetics.pdf", dpi=300, bbox_inches="tight")
    # magFig.savefig(f"{args.o}/magnetics.pdf", dpi=300, bbox_inches="tight")


######################################################################
