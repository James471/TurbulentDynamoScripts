from constants import *
import numpy as np
import pickle
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
    if n == 1:
        return np.loadtxt(path, unpack=True, skiprows=shift)
    else:
        with open(path, "r") as f:
            iterator = itertools.islice(f, shift, None, n)
            return np.loadtxt(iterator, unpack=True)


def getInfoDict(dirPath):
    # if os.path.exists(args.i + "/info.pkl"):
    with open(dirPath + "/info.pkl", "rb") as f:
        return pickle.load(f)
    # return None
