import os
import pickle
from myconfig import *


E_MAG_COLUMN_INDEX = 11
E_KIN_COLUMN_INDEX = 9
V_RMS_COLUMN_INDEX = 13
TIME_COLUMN_INDEX = 0


SOLVER_DICT = {"8wave": "Split-Roe", "bouchut-split": "Split-Bouchut", "Roe": "USM-Roe", 
              "HLLD": "USM-HLLD", "HLLC": "USM-HLLC", "bk-usm": "USM-BK"}
COLOR_DICT = {"8wave": "#377eb8", "HLLD": "#984ea3", "HLLC": "#4daf4a", 
              "Roe": "#a65628", "bouchut-split": "#f781bf", "bk-usm": "#ff7f00"}
ORDER_DICT = {"8wave": 0, "bouchut-split": 1, "Roe": 2, "HLLD": 3, "HLLC": 4, "bk-usm": 5}
FACT_DICT = {"bk-usm": 1.0, "Roe": 0.1, "HLLC": 0.01, "HLLD": 0.001, "bouchut-split": 0.0001, "8wave": 0.00001}

def getInfoDict(sim):
    '''
    sim: Path to simulation directory
    Raises exception if pkl file not found
    '''
    if os.path.exists(sim + "/info.pkl"):
        with open(sim + "/info.pkl", "rb") as f:
            return pickle.load(f)
    else:
        raise Exception("No info.pkl file found in " + sim)
    

def argsToOutdirName(args):
    outputDirName = (
        args.outdir
        + "/Turb_v"
        + str(args.v)
        + "_auto-adj"
        + str(args.auto_adjust)
        + "_visc-"
        + str(args.useVisc)
        + "_Re-"
        + str(args.Re)
        + "_mgRes-"
        + str(args.useMgRes)
        + "_Prm-"
        + str(args.Prm)
        + "_sol-wt"
        + str(args.sol_weight)
        + "_solver-"
        + str(args.solver)
        + "_emd-"
        + str(args.E_method)
        + "_mcut-"
        + str(args.mcut)
        + "_cfl-"
        + str(args.cfl)
        + "_nt-"
        + str(args.nt)
        + "_dt-"
        + str(args.dt)
        + "_iprocs-"
        + str(args.iprocs)
        + "_jprocs-"
        + str(args.jprocs)
        + "_kprocs-"
        + str(args.kprocs)
        + "_nxb-"
        + str(args.nxb)
        + "_nyb-"
        + str(args.nyb)
        + "_nzb-"
        + str(args.nzb)
    )
    if args.extra != None:
        outputDirName += "_" + args.extra
    return outputDirName


def argsToSimulationObjectDirectory(args):
    return argsToOutdirName(args) + "/objStirFromFile"


Object = lambda **kwargs: type("Object", (), kwargs)


def parseArgs(args, commonKeys):
    numSim = -1
    for key, value in vars(args).items():
        if key not in commonKeys:
            if value is not None and len(value) >= 1:
                numSim = len(value)
                break

    for key, value in vars(args).items():
        if key not in commonKeys:
            if value is None:
                setattr(args, key, [None for _ in range(numSim)])
            elif len(value) == 1:
                setattr(args, key, [value[0] for _ in range(numSim)])

    return args