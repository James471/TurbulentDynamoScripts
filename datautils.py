from constants import *
import numpy as np
import pickle

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


def loadFile(path, shift):
    return np.loadtxt(path, unpack=True, skiprows=shift)


def getInfoDict(dirPath):
    # if os.path.exists(args.i + "/info.pkl"):
    with open(dirPath + "/info.pkl", "rb") as f:
        return pickle.load(f)
    # return None
