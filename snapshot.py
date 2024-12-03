#!/usr/bin/env python
import argparse
import numpy as np
import textwrap
import subprocess
import os
import re
import sys
import h5py
from decimal import Decimal
from ipdb import set_trace as stop

from datautils import *
from constants import *
from utils import *
from myconfig import *

import flashlib


X_LABEL_DICT = {"Split-Roe": "", "Split-Bouchut": "", "USM-Roe": "",
                "USM-HLLD": "\$x/L\$", "USM-HLLC": "\$x/L\$", "USM-BK": "\$x/L\$"}

Y_LABEL_DICT = {"Split-Roe": "\$y/L\$", "Split-Bouchut": "", "USM-Roe": "",
                "USM-HLLD": "\$y/L\$", "USM-HLLC": "", "USM-BK": ""}

X_TICK_DICT = {"Split-Roe": "\"\"", "Split-Bouchut": "\"\"", "USM-Roe": "\"\"",
                "USM-HLLD": None, "USM-HLLC": None, "USM-BK": None}

Y_TICK_DICT = {"Split-Roe": None, "Split-Bouchut": "\"\"", "USM-Roe": "\"\"",
                "USM-HLLD": None, "USM-HLLC": "\"\"", "USM-BK": "\"\""}



def getSnapshotInfo(turbFile, velocity, val, usetime, stf, lf, simDir, type):
    t = getNonDimensionalTime(turbFile, velocity)
    ratio = getEMagOverEKin(turbFile)
    transientMask = (t > stf) & (ratio > lf)
    if usetime:
        location = np.argmin(np.abs(t[transientMask] - val))
    else:
        location = np.argmin(np.abs(ratio[transientMask] - val))
    snapshotTime = t[transientMask][location]

    plotFileList = [f for f in os.listdir(simDir) if f.startswith("Turb_hdf5_plt_cnt_")]
    plotFileList.sort()
    timeList = []
    for plotFile in plotFileList:
        dataset = hdf5ToDict(simDir + "/" + plotFile, "real scalars")
        timeList.append(dataset["time"]/getTurnOverTime(velocity))
    
    timeList = np.array(timeList)
    plotFile = plotFileList[np.argmin(np.abs(timeList - snapshotTime))]
    plotFileObj = flashlib.FlashGG(simDir + "/" + plotFile)
    plotFileTime = plotFileObj.scalars["time"]

    if type == "mag":
        eTot = turbFile[E_MAG_COLUMN_INDEX][np.argmin(np.abs(turbFile[TIME_COLUMN_INDEX]-plotFileTime))]
    elif type == "kin":
        eTot = turbFile[E_KIN_COLUMN_INDEX][np.argmin(np.abs(turbFile[TIME_COLUMN_INDEX]-plotFileTime))]

    print(f"Snapshot Time: {snapshotTime}t_turb")
    print(f"Plotfile: {plotFile}")
    print(f"Minimum discrepancy: {(snapshotTime - timeList)[np.argmin(np.abs(timeList - snapshotTime))]}t_turb")

    return plotFile, eTot


def makePlots(simDirList, type, stream_var, val, stf, lf, outdir, cbar_type, redo, usetime=False, fontsize=1, clow=None, chigh=None):
    
    plotFileDict = {}
    boundDict = {}

    if type == 'mag':
        datasetName = 'emag'
    elif type == 'kin':
        datasetName = 'ekin'

    if not redo and os.path.exists(f"{outdir}/{type}Bounds.pkl") and os.path.exists(f"{outdir}/plotFiles_{type}.pkl"):
        boundDict = loadDict(f"{outdir}/{type}Bounds.pkl")
        plotFileDict = loadDict(f"{outdir}/plotFiles_{type}.pkl")
        print(f"Loaded {type} bounds:", boundDict)
        print(f"Loaded {type} dict:", plotFileDict)
    else:
        for index, simDir in enumerate(simDirList):

            infoDict = getInfoDict(simDir)
            velocity = infoDict["v"]
            solver = SOLVER_DICT[infoDict["solver"]]

            print(f"Processing for {solver}")

            print("Loading Turb.dat")
            n, s = getNForTurbDat(simDir, res=5e-1, stop=10)
            turbFile = loadFile(simDir + "/Turb.dat", 10, n, s)
            print("Loaded Turb.dat")

            plotFile, eTot = getSnapshotInfo(turbFile, velocity, val, usetime, stf[index], lf[index], simDir, type)
            plotFile = simDir + "/" + plotFile
            if type == 'mag':
                plotFileDict[solver] = {"plotFile": plotFile, "eMagTot": eTot}
            elif type == 'kin':
                plotFileDict[solver] = {"plotFile": plotFile, "eKinTot": eTot}

            cmd = f'python3 {PYTHON_PATH}/flashplotlib.py -i {plotFile} -d {datasetName} -fontsize {fontsize} -verbose 2 -outtype pdf -outdir {outdir} -outname {solver}-{type} -direction z -ncpu 1 -colorbar 0 -data_transform q/{eTot}'
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            searchString = "flashplotlib.py: prep_map: min, max of data.*"
            iterms = re.findall(searchString, output)
            boundString = iterms[0]
            boundString = boundString.replace("flashplotlib.py: prep_map: min, max of data = ", "")
            minBound, maxBound = boundString.split(",")
            minBound, maxBound = float(minBound), float(maxBound)
            boundDict[solver] = (minBound, maxBound)
            print(f"Bounds acquired for {solver}")

        dumpDict(boundDict, f"{outdir}/{type}Bounds.pkl")
        dumpDict(plotFileDict, f"{outdir}/plotFiles_{type}.pkl")

    print(f"{datasetName} bounds:", boundDict)
    print(f"plot files:", plotFileDict)

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
        plotFile = plotFileDict[solver]["plotFile"]
        eTot = plotFileDict[solver][f"e{type.capitalize()}Tot"]

        cmd = f'python3 {PYTHON_PATH}/flashplotlib.py -fontsize {fontsize} -i {plotFile} -d {datasetName} -verbose 2 -outtype pdf -outdir {outdir} -outname {solver}-{type} -direction z -ncpu 1 -colorbar 0 -vmin {vMin} -vmax {vMax} -axes_label "{X_LABEL_DICT[solver]}" "{Y_LABEL_DICT[solver]}" "" -axes_unit "" "" "" -plotlabel {solver} -time_unit t_{{turb}} -labels_inside -time_scale {getTurnOverTime(velocity)} -axes_format {X_TICK_DICT[solver]} {Y_TICK_DICT[solver]} -data_transform q/{eTot} -cmap {cbar_type} {streamStr} &>/dev/null'
        os.system(cmd)

    print("Making colorbar")
    if clow is None:
        pwr_low = int((f"{Decimal(vMin):.1E}").split("E")[1])
        vmin = 10**pwr_low
    else:
        vmin = clow
    if chigh is None:
        pwr_high = int((f"{Decimal(vMax):.1E}").split("E")[1])+1
        vmax = 10**pwr_high
    else:
        vmax = chigh
    cbarCmd = f'python3 {PYTHON_PATH}/flashplotlib.py -fontsize {fontsize} -i {simDirList[0]}/Turb_hdf5_plt_cnt_0000 -d dens -verbose 2 -outtype pdf -outdir {outdir} -outname {type}Bar -ncpu 1 -colorbar "only" -vmin {vmin} -vmax {vmax} -cmap_label "\$E_\\mathrm{{{type}}}/E_\\mathrm{{{type},\,total}}\$ (Projection length \$z=L\$)" -cmap {cbar_type} &>/dev/null'
    if len(simDirList) > 3:
        cbarCmd = f'python3 {PYTHON_PATH}/flashplotlib.py -fontsize {fontsize} -i {simDirList[0]}/Turb_hdf5_plt_cnt_0000 -d dens -verbose 2 -outtype pdf -outdir {outdir} -outname {type}Bar -ncpu 1 -colorbar "panels2 only" -vmin {vmin} -vmax {vmax} -cmap_label "\$E_\\mathrm{{{type}}}/E_\\mathrm{{{type},\,total}}\$ (Projection length \$z=L\$)" -cmap {cbar_type} &>/dev/null'
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
    # os.system(f"pdflatex -interaction=nonstopmode --halt-on-error -output-directory {args.o} {args.o}/{type}.tex")
    # os.system(f"rm {outdir}/{type}.aux {outdir}/{type}.log {outdir}/{type}.tex")



def main(args):

    simList = [sim for sim in args.i if os.path.isdir(sim)]
    outdir = args.o

    if args.t is None and args.r is None:
        print("Invalid arguments. Please provide either -t or -r")
        sys.exit(1)
    if args.t is not None and args.r is not None:
        print("Invalid arguments. Please provide either -t or -r")
        sys.exit(1)

    if args.kin:
        if args.t is not None:
            val = args.t
            usetime = True
        else:
            val = args.r
            usetime = False
        makePlots(simList, "kin", args.kin_stream, val, args.stf, args.lf, outdir, args.cbar_type, args.redo, usetime, args.fontsize, args.clow, args.chigh)
        if args.table:
            makeTable(simList, outdir, "kin")
    if args.mag:
        if args.t is not None:
            val = args.t
            usetime = True
        else:
            val = args.r
            usetime = False
        makePlots(simList, "mag", args.mag_stream, val, args.stf, args.lf, outdir, args.cbar_type, args.redo, usetime, args.fontsize, args.clow, args.chigh)
        if args.table:
            makeTable(simList, outdir, "mag")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="./", help="Output Directory")
    parser.add_argument("-r", type=float, default=None, help="Energy ratio for projections")
    parser.add_argument("-t", type=float, default=None, help="Time for projections")
    parser.add_argument("-stf", type=float, nargs="*", default=[3.0], help="Start time for looking for ratio")
    parser.add_argument("-lf", type=float, nargs="*", default=[1e-7], help="Lower bound for looking for ratio")
    parser.add_argument("-kin", action='store_true', help="Plot kinetic energy")
    parser.add_argument("-mag", action='store_true', help="Plot magnetic energy")
    parser.add_argument("-table", action='store_true', help="Want the plots to be arranged in a table")
    parser.add_argument("-cbar_type", type=str, default="magma", help="Colorbar type")
    parser.add_argument("-kin_stream", type=str, default=None, choices=['vel', 'mag'], help="Streamline on kinetic energy")
    parser.add_argument("-mag_stream", type=str, default=None, choices=['vel', 'mag'], help="Streamline on magnetic energy")
    parser.add_argument("-clow", type=float, default=None, help="Lower bound for colorbar")
    parser.add_argument("-chigh", type=float, default=None, help="Upper bound for colorbar")
    parser.add_argument("-redo", action='store_true', help="Redo the bounds calculation")
    parser.add_argument("-fontsize", type=float, default=1.5, help="Fontsize to pass to flashplotlib")

    commonKeys = ["r", "o", "kin", "mag", "cbar_type", "kin_stream", "mag_stream", "table", "redo", "fontsize", "clow", "chigh", "t"]

    args = parseArgs(parser.parse_args(), commonKeys)

    print(args)

    main(args)