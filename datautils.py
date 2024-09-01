from constants import *
import numpy as np
import pickle
import os
import h5py
import itertools

def getTurnOverTime(v):
    return 1/(2*v)


def getNonDimensionalTime(f, v):
    turnOverTime = getTurnOverTime(v)
    return (f[TIME_COLUMN_INDEX] / turnOverTime).flatten()


def getEMag(f):
    return f[E_MAG_COLUMN_INDEX].flatten()


def getEKin(f):
    return f[E_KIN_COLUMN_INDEX].flatten()


def getVRMS(f):
    return f[V_RMS_COLUMN_INDEX].flatten()


def getCsRMS(f):
    return f[Cs_RMS_COLUMN_INDEX].flatten()


def getEMagOverEKin(f):
    return (getEMag(f) / getEKin(f)).flatten()


def loadFile(path, shift=0, n=1):
    #Loads every nth row
    if n == 1:
        return np.loadtxt(path, unpack=True, skiprows=shift)
    else:
        with open(path, "r") as f:
            iterator = itertools.islice(f, shift, None, n)
            return np.loadtxt(iterator, unpack=True)
        
def hdf5ToDict(path, datasetname):
    f = h5py.File(path, "r")
    d = f[datasetname]
    dct = {}
    for i, key in enumerate(d["name"]):
        dct[key.decode("utf-8").strip()] = d["value"][i]
    return dct

    
def getAverageDt(simPath):
    chkFileList = [f for f in os.listdir(simPath) if f.startswith("Turb_hdf5_chk_")]
    dtList = []
    for plotFile in chkFileList:
        dt = hdf5ToDict(f"{simPath}/{plotFile}", "real scalars")["dt"]
        dtList.append(dt)
    return np.mean(np.array(dtList))


def getNForTurbDat(simPath, res=1e-2, filename="Turb.dat"):
    
    infoDict = getInfoDict(simPath)
    v = infoDict["v"]
    nt = infoDict["nt"]
    tTurb = getTurnOverTime(v)
    os.system(f"wc {simPath}/{filename} -l > {simPath}/temp.txt")
    lineCount = int(open(f"{simPath}/temp.txt", "r").read().split()[0])
    dt = getAverageDt(simPath)
    print(f"Line count is {lineCount}")
    if lineCount < 50_000 or dt > res * tTurb:
        return 1
    else:
        return round(lineCount * res / nt)
        

def getInfoDict(dirPath):
    # if os.path.exists(args.i + "/info.pkl"):
    with open(dirPath + "/info.pkl", "rb") as f:
        return pickle.load(f)
    # return None
