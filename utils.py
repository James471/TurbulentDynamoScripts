import os
import pickle
import numpy as np
from constants import *


def getSolverSortedList(pathList):
    order = []
    cleanPathList = [x for x in pathList if os.path.isdir(x)]
    for sim in cleanPathList:
        if not os.path.exists(sim + "/info.pkl"):
            raise Exception("No info.pkl file found in " + sim)
        infoDict = getInfoDict(sim)
        order.append(ORDER_DICT[infoDict["solver"]])
    srt = np.argsort(np.array(order))
    return [cleanPathList[i] for i in srt]

def txtToCfpDict(path):
    with open(path, "r") as f:
        lines = f.readlines()
    cfpDict = {}
    for line in lines:
        key = str(line.split(":")[0])
        temp = line.split(":")[1].split(",")
        if line == "":
            continue
        values = []
        for i in temp:
            i = i.strip()
            if i == "inf":
                value = np.inf
            elif i == "-inf":
                value = -np.inf
            else:
                value = float(i)
            values.append(value)
        cfpDict[key] = values
    return cfpDict


def dumpDict(dictionary, filename):
    with open(filename, "wb") as f:
        pickle.dump(dictionary, f)


def loadDict(filePath):
    with open(filePath, "rb") as f:
        return pickle.load(f)


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