import argparse
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import numpy as np
import pickle


E_MAG_COLUMN_INDEX = 11
E_KIN_COLUMN_INDEX = 9
V_RMS_COLUMN_INDEX = 13
Cs_RMS_COLUMN_INDEX = 14
TIME_COLUMN_INDEX = 0


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


def loadFile(path, shift=0):
    return np.loadtxt(path, unpack=True, skiprows=shift)


def getInfoDict(dirPath):
    # if os.path.exists(args.i + "/info.pkl"):
    with open(dirPath + "/info.pkl", "rb") as f:
        return pickle.load(f)
    # return None


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


def getFigGrid():
    fig = pl.figure()
    gs = gridspec.GridSpec(2, 3)
    axSplitRoe = fig.add_subplot(gs[0, 0])
    axSplitBouchut = fig.add_subplot(gs[0, 1])
    axUSMRoe = fig.add_subplot(gs[0, 2])
    axUSMHLLD = fig.add_subplot(gs[1, 0])
    axUSMHLLC = fig.add_subplot(gs[1, 1])
    axUSMBK = fig.add_subplot(gs[1, 2])
    return fig, [axSplitRoe, axSplitBouchut, axUSMRoe, axUSMHLLD, axUSMHLLC, axUSMBK]


def getPlotFile(sim, velocity, r, stf, lf, dt):
    
    turbFile = loadFile(sim + "/Turb.dat")
    t = getNonDimensionalTime(turbFile, velocity)
    ratio = getEMagOverEKin(turbFile)
    
    transientMask = (t > stf) & (t < lf)
    location = np.argmin(np.abs(ratio[transientMask] - r))
    snapshotTime = t[transientMask][location]
    plotFileNumber = int(round(snapshotTime / dt))

    return "Turb_hdf5_plt_cnt_" + str(plotFileNumber).zfill(4)


def main(args):

    simList = args.i
    
    figKin, axesKin = getFigGrid()
    figMag, axesMag = getFigGrid()
    solverDict = {"8wave": "Split-Roe", "bouchut-split": "Split-Bouchut", "Roe": "USM-Roes", 
                  "HLLD": "USM-HLLD", "HLLC": "USM-HLLC", "bk-usm": "USM-BK"}
    kinDict = {"Split-Roe": axesKin[0], "Split-Bouchut": axesKin[1], "USM-Roe": axesKin[2], 
               "USM-HLLD": axesKin[3], "USM-HLLC": axesKin[4], "USM-BK": axesKin[5]}
    magDict = {"Split-Roe": axesMag[0], "Split-Bouchut": axesMag[1], "USM-Roe": axesMag[2], 
               "USM-HLLD": axesMag[3], "USM-HLLC": axesMag[4], "USM-BK": axesMag[5]}

    for index, sim in enumerate(simList):
        infoDict = getInfoDict(sim)
        velocity = infoDict["v"]
        solver = solverDict[infoDict["solver"]]
        dt = infoDict["dt"]
        plotFile = getPlotFile(sim, velocity, args.r, args.stf, args.lf, dt)
        print("Plot file:", sim + plotFile)
        kinAx = kinDict[solver]
        magAx = magDict[solver]

        


commonKeys = [""]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="Turb.pdf", help="Output Directory")
    parser.add_argument("-r", type=float, default=1e-6, help="Energy ratio for projections")
    parser.add_argument("-stf", type=float, nargs="*", help="Start time for looking for ratio")
    parser.add_argument("-lf", type=float, nargs="*", help="Lower bound for looking for ratio")
    parser.add_argument("-projection_mpi_path", type=str, default="/home/james471/Academics/Projects/MHD/Code/tools/projection_mpi/projection_mpi", help="Path to mpiProjection")

    args = parseArgs(parser.parse_args())
    print(args)

    main(args)

