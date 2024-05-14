import numpy as np  
import argparse
import matplotlib.pyplot as pl
import sys
import os
import glob
import pandas as pd
from common import *
import designParams
sys.path.append(PYTHON_PATH)
import cfpack as cfp
import flashlib as fl
from decimal import Decimal
import turblib as tl
from scipy.special import kn
from scipy.optimize import curve_fit
from scipy.optimize import fminbound
from scipy.interpolate import CubicSpline


def Re_fun(k_tilde_nu, c_nu):
    k_nu = k_tilde_nu
    return (k_nu/(c_nu * 2))**(4/3)


def k_eta_fun(k_tilde_eta, p_eta):
    return (k_tilde_eta)**(1/p_eta)


def Rm_fun(k_tilde_nu, c_nu, k_tilde_eta, p_eta, c_eta):
    k_eta = k_eta_fun(k_tilde_eta, p_eta)
    k_nu = k_tilde_nu
    Pm = (k_eta/(c_eta*k_nu))**2
    return Pm * Re_fun(k_tilde_nu, c_nu)

def getNum(num):
    return float((f"{Decimal(num):.1E}").split("E")[0])


def getPwr(num):
    return int((f"{Decimal(num):.1E}").split("E")[1])


def Log10_P_kin(k, A_kin, p_bn, k_bn, k_tilde_nu):
    p_kin = -1.7
    p_nu  = 1.0
    return np.log10(A_kin * ((k / k_bn)**p_kin + (k / k_bn)**p_bn) * np.exp(-(k / k_tilde_nu)**p_nu))


def Log10_P_mag(k, A_mag, p_mag, p_eta, k_tilde_eta):
    y = A_mag * (k**p_mag) * kn(0, (k / k_tilde_eta)**p_eta)
    return np.log10(y)


def generateSpectra(simDir, verbose, spectType, infoDict, lf, uf, stf, n=1):
    '''
    verbsoe  : Could be 0 or 1 or 2
    spectType: Could be 'vels' or 'curr' or 'mags'
    '''

    if verbose == 2:
        datFileVerbose = 1
        devNull = " > /dev/null 2>&1"
    else:
        datFileVerbose = 0
        devNull = ""

    spectraDir = simDir + "/spectra"
    if not os.path.isdir(spectraDir):
        print("creating spectra dir: '"+spectraDir+"'")
        os.mkdir(spectraDir)

    # NOTE: You made some changes here and it's different from how things were done earlier. You just removed [:101] from the end. If things break, look here first.
    files = sorted(glob.glob(simDir+"/Turb_hdf5_plt_cnt_????"))

    dt = infoDict['dt']
    t_turb = 1 / (2 * infoDict['v'])
    if verbose > 0: print("Loading turb.dat for cleaning")
    datfile_obj = fl.datfile(simDir+"/Turb.dat", verbose=datFileVerbose)
    datfile_obj.write_cleaned()
    if verbose > 0: print("Turb.dat cleaned")

    if verbose > 0: print("Loading cleaned turb.dat")
    dat = np.loadtxt(simDir + "/Turb.dat_cleaned", unpack=True, skiprows=10)
    if verbose > 0: print("Loaded Turb.dat_cleaned")

    emag_over_ekin = dat[E_MAG_COLUMN_INDEX].flatten() / dat[E_KIN_COLUMN_INDEX].flatten()
    t = dat[TIME_COLUMN_INDEX].flatten()
    ind = (emag_over_ekin >= lf) & (emag_over_ekin <= uf) & (t/t_turb >= stf)
    t_valid = t[ind]

    dumpStart = fl.get_dump(files, time=t_valid[0])
    dumpNStart = int(dumpStart[0][-4:])
    dumpEnd = fl.get_dump(files, time=t_valid[-1])
    dumpNEnd = int(dumpEnd[0][-4:])

    if verbose > 0: 
        print("Start dump:", dumpNStart)
        print("End dump:", dumpNEnd)

    if spectType == "cur":
        types = "0 -dsets curx cury curz"
        for x in range(dumpNStart, dumpNEnd+1):
            if verbose > 0: print("Adding current variable to dump:", x)
            f = simDir+"/Turb_hdf5_plt_cnt_{:04d}".format(x)
            derivCmd = f"mpirun -np {n} /home/100/jw5893/flash-tools/tools/derivative_var/derivative_var {f} -current {devNull}"
            os.system(derivCmd)
    elif spectType == "vels":
        types = "1"
    elif spectType == "mags":
        types = "2"

    if spectType == "cur":
        spectFileString = "dset_curx_cury_curz"
    elif spectType == "mags":
        spectFileString = "mags"
    elif spectType == "vels":
        spectFileString = "vels"

    for x in range(dumpNStart, dumpNEnd+1):
        if verbose > 0: print("Processing dump to generate spectra:", x)
        f = simDir+"/Turb_hdf5_plt_cnt_{:04d}".format(x)
        spectrCmd = f"mpirun -np {n} /home/100/jw5893/flash-tools/tools/spectra_mpi/spectra_mpi {f} -types {types} {devNull}"
        os.system(spectrCmd)
    mvCmd = "mv "+simDir+"/Turb_hdf5_plt_cnt_????_spect_*.dat "+spectraDir
    os.system(mvCmd)

    spectFiles = [simDir+"/spectra/Turb_hdf5_plt_cnt_{:04d}".format(x)+"_spect_"+spectFileString+".dat" for x in range(dumpNStart, dumpNEnd+1)]
    if spectType == "vels": normalise = False
    if spectType == "mags" or spectType == "cur": normalise = True # normalise magnetic spectra and current spectra before averaging, because of field growth with time
    averDat, headerAver = tl.aver_spect(spectFiles, normalise=normalise, verbose=0)
    outfile = simDir+"/spectra/aver_spect_"+spectFileString+".dat"
    tl.write_spect(outfile, averDat, headerAver, verbose=verbose)


def plotSpectra(simDir, verbose, spectType, fact, infoDict, params, scp=False, compensate=False):
    '''
    verbsoe  : Could be 0 or 1
    spectType: Could be 'vels' or 'mags' or 'curr'
    '''

    if spectType == "mags":
        spectVar = "mags"
    elif spectType == "vels":
        spectVar = "vels"
    elif spectType == "cur":
        spectVar = "dset_curx_cury_curz"

    averFile = simDir + f"/spectra/aver_spect_{spectVar}.dat"
    averDf   = pd.read_csv(averFile, sep="\s+", header=0, skiprows=5)
    
    log10PTot = averDf['log10_Ptot']
    deltaLog10PTot = averDf['sigma_log10_Ptot']
    pTot = 10**log10PTot
    k = averDf['k']
    err = np.array([pTot * np.log(10) * deltaLog10PTot, pTot * np.log(10) * deltaLog10PTot])
    
    kFitMin = 3
    kFitMax = max(k) / 2
    kFit = np.array(k[(k >= kFitMin) & (k <= kFitMax)].tolist())
    log10PTotFit = np.array(log10PTot[(k >= kFitMin) & (k <= kFitMax)].tolist())
    deltaLog10PTotFit = np.array(deltaLog10PTot[(k >= kFitMin) & (k <= kFitMax)].tolist())

    if compensate and spectType == "vels":
        compensateFact = k**2
        compensateFitFact = kFit**2

    else:
        compensateFact = 1
        compensateFitFact = 1

    if compensate and spectType == "vels":
        plObj = cfp.plot(compensateFact * fact * pTot, k, yerr=fact * err * np.array([compensateFact, compensateFact]), shaded_err=True, label=SOLVER_DICT[infoDict['solver']], color=COLOR_DICT[infoDict['solver']])
    else:
        plObj = cfp.plot(compensateFact * fact * pTot, k, yerr=fact * err, shaded_err=True, label=SOLVER_DICT[infoDict['solver']], color=COLOR_DICT[infoDict['solver']])

    fitDict = {}

    if spectType == "mags":
        if verbose: print("Fitting for:", infoDict['solver'])
        if scp is False:
            magFit = cfp.fit(Log10_P_mag, kFit, log10PTotFit, xerr=None, yerr=deltaLog10PTotFit,
                              params=params, max_nfev=100000, n_random_draws=10000)
            cfp.plot(compensateFitFact * fact * 10**Log10_P_mag(kFit, *(magFit.popt)), kFit, color="black")
            fitDict["A_mag"] = (magFit.popt[0], magFit.perr[0][1], magFit.perr[0][0])
            fitDict["p_mag"] = (magFit.popt[1], magFit.perr[1][1], magFit.perr[1][0])
            fitDict["p_eta"] = (magFit.popt[2], magFit.perr[2][1], magFit.perr[2][0])
            fitDict["k_tilde_eta"] = (magFit.popt[3], magFit.perr[3][1], magFit.perr[3][0])
        else:
            A_mag = params['A_mag']
            p_mag = params['p_mag']
            p_eta = params['p_eta']
            k_tilde_eta = params['k_tilde_eta']
            magPopt, magCor = curve_fit(Log10_P_mag, kFit, log10PTotFit, sigma=deltaLog10PTotFit, absolute_sigma=True,
                                         p0=[A_mag[1], p_mag[1], p_eta[1], k_tilde_eta[1]], 
                                         bounds=([A_mag[0], p_mag[0], p_eta[0], k_tilde_eta[0]], [A_mag[2], p_mag[2], p_eta[2], k_tilde_eta[2]]))
            cfp.plot(fact * 10**Log10_P_mag(kFit, *magPopt) * compensateFitFact, kFit, color="black")
            print(f"A_mag: {magPopt[0]}±{magCor[0, 0]**0.5}")
            print(f"p_mag: {magPopt[1]}±{magCor[1, 1]**0.5}")
            print(f"p_eta: {magPopt[2]}±{magCor[2, 2]**0.5}")
            print(f"k_tilde_eta: {magPopt[3]}±{magCor[3, 3]**0.5}")
            fitDict["A_mag"] = (magPopt[0], magCor[0, 0]**0.5, -magCor[0, 0]**0.5)
            fitDict["p_mag"] = (magPopt[1], magCor[1, 1]**0.5, -magCor[1, 1]**0.5)
            fitDict["p_eta"] = (magPopt[2], magCor[2, 2]**0.5, -magCor[2, 2]**0.5)
            fitDict["k_tilde_eta"] = (magPopt[3], magCor[3, 3]**0.5, -magCor[3, 3]**0.5)

    elif spectType == "cur":
        infoDict = getInfoDict(simDir)
        spectraDir = simDir+"/spectra"
        kList = []
        for file in os.listdir(spectraDir):
            if "curx" in file and "aver" not in file:
                spectraFile = spectraDir + "/" + file
                plotFileNum = file.split("_")[4].split("_")[0]
                spectraDf   = pd.read_csv(spectraFile, sep="\s+", header=0, skiprows=5)
                pTot = spectraDf['#15_SpectFunctTot']
                k = spectraDf['#01_KStag']
                spl = CubicSpline(k, pTot)
                k1, k2 = k[0], k[len(k)-1]
                kList.append(fminbound(lambda x: -spl(x), k1, k2))
                # kList.append(k[pTot==max(pTot)])

        # print(np.percentile(kList, 84), np.percentile(kList, 50), np.percentile(kList, 16))
        # print(kList)
        k_eta_23 = np.percentile(kList, 50)
        k_eta_23_pos_er = np.percentile(kList, 84) - k_eta_23
        k_eta_23_neg_er = k_eta_23 - np.percentile(kList, 16)
        fitDict["k_eta_23"] = (k_eta_23, k_eta_23_pos_er, -k_eta_23_neg_er)
        print(f"{SOLVER_DICT[infoDict['solver']]}: {k_eta_23:.2f}^{{{k_eta_23_pos_er:.2f}}}_{{{-k_eta_23_neg_er:.2f}}}")
        # ylim = plCurObj.gca().get_ylim()
        # plCurObj.plot([k_nu_23, k_nu_23], [ylim[0], 2e-10], color=colorDict[infoDict['solver']], scaley=False)

    elif spectType == "vels":
        if verbose: print("Fitting for:", infoDict['solver'])
        if scp is False:
            kinFit = cfp.fit(Log10_P_kin, kFit, log10PTotFit, xerr=None, yerr=deltaLog10PTotFit, params=params, n_random_draws=10000)
            cfp.plot(fact * 10**Log10_P_kin(kFit, *(kinFit.popt)) * compensateFitFact, kFit, color="black")
            fitDict["A_kin"] = (kinFit.popt[0], kinFit.perr[0][1], kinFit.perr[0][0])
            fitDict["p_bn"] = (kinFit.popt[1], kinFit.perr[1][1], kinFit.perr[1][0])
            fitDict["k_bn"] = (kinFit.popt[2], kinFit.perr[2][1], kinFit.perr[2][0])
            fitDict["k_tilde_nu"] = (kinFit.popt[3], kinFit.perr[3][1], kinFit.perr[3][0])
        else:
            A_kin = params['A_kin']
            p_bn = params['p_bn']
            k_bn = params['k_bn']
            k_tilde_nu = params['k_tilde_nu']
            kinPopt, kinCor = curve_fit(Log10_P_kin, kFit, log10PTotFit, sigma=deltaLog10PTotFit, absolute_sigma=True, 
                                        p0=[A_kin[1], p_bn[1], k_bn[1], k_tilde_nu[1]], maxfev=10000,
                                        bounds=([A_kin[0], p_bn[0], k_bn[0], k_tilde_nu[0]], [A_kin[2], p_bn[2], k_bn[2], k_tilde_nu[2]]))
            cfp.plot(fact * 10**Log10_P_kin(kFit, *kinPopt) * compensateFitFact, kFit, color="black")
            print(f"A_kin: {kinPopt[0]}±{kinCor[0, 0]**0.5}")
            print(f"p_bn: {kinPopt[1]}±{kinCor[1, 1]**0.5}")
            print(f"k_bn: {kinPopt[2]}±{kinCor[2, 2]**0.5}")
            print(f"k_tilde_nu: {kinPopt[3]}±{kinCor[3, 3]**0.5}")
            fitDict["A_kin"] = (kinPopt[0], kinCor[0, 0]**0.5, -kinCor[0, 0]**0.5)
            fitDict["p_bn"] = (kinPopt[1], kinCor[1, 1]**0.5, -kinCor[1, 1]**0.5)
            fitDict["k_bn"] = (kinPopt[2], kinCor[2, 2]**0.5, -kinCor[2, 2]**0.5)
            fitDict["k_tilde_nu"] = (kinPopt[3], kinCor[3, 3]**0.5, -kinCor[3, 3]**0.5)

    return plObj, fitDict

def postPlot(plObj, spectType, compensated=False):
    ax = plObj.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    if spectType == "mags":
        ax.set_ylabel(r'$P_\mathrm{mag}$')
        ax.set_xlabel(r'$k$')
        # ax.legend(loc='best')
    elif spectType == "vels":
        if compensated:
            ax.set_ylabel(r'$k^2P_\mathrm{kin}$')
        else:
            ax.set_ylabel(r'$P_\mathrm{kin}$')
        ax.get_xaxis().set_ticks([])
        ax.legend(loc='best')
    elif spectType == "cur":
        ax.set_ylabel(r'$P_\mathrm{cur}$')
        ax.set_xlabel(r'$k$')
        # ax.legend(loc='best')

def plotScaleLoc(plObj, solverFit, type):
    ax = plObj.gca()
    ylim = ax.get_ylim()
    if type == "mags":
        maxKEta = 0
        for solver in solverFit:
            val = solverFit[solver]
            ax.plot([val["k_tilde_eta"][0]**(1/val["p_eta"][0]), val["k_tilde_eta"][0]**(1/val["p_eta"][0])], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
            if val["k_tilde_eta"][0]**(1/val["p_eta"][0]) > maxKEta:
                maxKEta = val["k_tilde_eta"][0]**(1/val["p_eta"][0])
        ax.text(maxKEta, 2.5*ylim[0], r"$k_\eta$", color="black")
    elif type == "vels":
        maxKNu = 0
        c_nu_22 = 0.025
        c_nu_23 = 0.116
        for solver in solverFit:
            val = solverFit[solver]
            k_nu_23 = val["k_tilde_nu"][0] * c_nu_23 / c_nu_22
            ax.plot([k_nu_23, k_nu_23], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
            if k_nu_23 > maxKNu:
                maxKNu = k_nu_23
        ax.text(maxKNu, 2.5*ylim[0], r"$k_\nu$", color="black")
    elif type == "cur":
        maxKEta = 0
        for solver in solverFit:
            val = solverFit[solver]
            ax.plot([val["k_eta_23"][0], val["k_eta_23"][0]], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
            if val["k_eta_23"][0] > maxKEta:
                maxKEta = val["k_eta_23"][0]
        ax.text(maxKEta, 2.5*ylim[0], r"$k_\eta$", color="black")

def main(args):

    simListTemp = [sim for sim in args.i if os.path.isdir(sim)]
    simList = ["" for i in simListTemp]
    # for i in simListTemp:
    #     infoDict = getInfoDict(i)
    #     simList[ORDER_DICT[infoDict['solver']]] = i
    simList=simListTemp

    if args.kin_spect:
        for sim in simList:
            generateSpectra(sim, args.v, "vels", getInfoDict(sim), args.lf, args.uf, args.stf, args.n)
    if args.cur_spect:
        for sim in simList:
            generateSpectra(sim, args.v, "cur", getInfoDict(sim), args.lf, args.uf, args.stf, args.n)
    if args.mag_spect:
        for sim in simList:
            generateSpectra(sim, args.v, "mags", getInfoDict(sim), args.lf, args.uf, args.stf, args.n)

    if args.kin_plot or args.table:
        solverKinFit = {}
        for simDir in simList:
            infoDict = getInfoDict(simDir)
            fact = FACT_DICT[infoDict['solver']]
            kinParams = {"A_kin": [0, 0.0015, np.inf], "p_bn": [0, 1, np.inf], "k_bn": [0, 4.0, 128], "k_tilde_nu": [0, 4.0, 128]}
            plKinObj, fitDict = plotSpectra(simDir, 1, "vels", fact, infoDict, kinParams, scp=False)
            solverKinFit[infoDict['solver']] = fitDict
            postPlot(plKinObj, "vels")
        plotScaleLoc(plKinObj, solverKinFit, "vels")
        plKinObj.savefig("Kinetic Spectra.pdf")
        plKinObj.clf(); plKinObj.cla(); plKinObj.close(); plKinObj = None

    
    if args.mag_plot or args.table:
        solverMagFit = {}
        for simDir in simList:
            infoDict = getInfoDict(simDir)
            fact = FACT_DICT[infoDict['solver']]
            magParams = {"A_mag": [0, 0.0001, np.inf], "p_mag": [0, 1, np.inf], "p_eta": [0, 1, np.inf], "k_tilde_eta": [0, 4.0, 128]}
            plMagObj, fitDict = plotSpectra(simDir, 1, "mags", fact, infoDict, magParams, scp=True)
            solverMagFit[infoDict['solver']] = fitDict
            postPlot(plMagObj, "mags")
        plotScaleLoc(plMagObj, solverMagFit, "mags")
        plMagObj.savefig("Magnetic Spectra.pdf")
        plMagObj.clf(); plMagObj.cla(); plMagObj.close(); plMagObj = None

    if args.cur_plot or args.table:
        solverCurFit = {}
        for simDir in simList:
            infoDict = getInfoDict(simDir)
            fact = FACT_DICT[infoDict['solver']]
            plCurObj, fitDict = plotSpectra(simDir, 1, "cur", fact, infoDict, None, scp=True)
            solverCurFit[infoDict['solver']] = fitDict
            postPlot(plCurObj, "cur")
        plotScaleLoc(plCurObj, solverCurFit, "cur")
        plCurObj.savefig("Current Spectra.pdf")
        plCurObj.clf(); plCurObj.cla(); plCurObj.close(); plCurObj = None

    if args.table:
        solverList = ["8wave", "bouchut-split", "Roe", "HLLD", "HLLC", "bk-usm"]
        table1Data = ''''''
        table2Data = ''''''

        c_nu = 0.025
        c_nu_pos_er = 0.005
        c_nu_neg_er = -0.006
        c_nu_pos = c_nu + c_nu_pos_er
        c_nu_neg = c_nu + c_nu_neg_er

        coef = 2.3
        coef_pos_er = 0.8
        coef_neg_er = -0.5
        coef_pos = coef + coef_pos_er
        coef_neg = coef + coef_neg_er

        for solver in solverList:
            kinFit = solverKinFit[solver]
            curFit = solverCurFit[solver]
            
            p_kin = -1.7
            p_nu = 1.0
            p_bn = kinFit["p_bn"][0]
            p_bn_pos_er = kinFit["p_bn"][1]
            p_bn_neg_er = kinFit["p_bn"][2]
            k_bn = kinFit["k_bn"][0]
            k_bn_pos_er = kinFit["k_bn"][1]
            k_bn_neg_er = kinFit["k_bn"][2]
            k_nu = kinFit["k_tilde_nu"][0]
            k_nu_pos_er = kinFit["k_tilde_nu"][1]
            k_nu_neg_er = kinFit["k_tilde_nu"][2]
            k_nu_pos = k_nu + k_nu_pos_er
            k_nu_neg = k_nu + k_nu_neg_er

            k_eta = curFit["k_eta_23"][0]
            k_eta_pos_er = curFit["k_eta_23"][1]
            k_eta_neg_er = curFit["k_eta_23"][2]
            k_eta_pos = k_eta + k_eta_pos_er
            k_eta_neg = k_eta + k_eta_neg_er

            Re = (k_nu / (c_nu * 2))**(4/3)
            Re_pos = ((k_nu_pos / (c_nu_neg * 2))**(4/3))
            Re_neg = ((k_nu_neg / (c_nu_pos * 2))**(4/3))
            Re_pos_er = Re_pos - Re
            Re_neg_er = Re_neg - Re
            Re_pos_er_disp = round(Re_pos_er / 10**getPwr(Re), 1)
            Re_neg_er_disp = round(Re_neg_er / 10**getPwr(Re), 1)

            Pm = (k_eta / (coef * k_nu))**2
            Pm_pos = ((k_eta_pos / (coef_neg * k_nu_neg))**2)
            Pm_neg = ((k_eta_neg / (coef_pos * k_nu_pos))**2)
            Pm_pos_er = Pm_pos - Pm
            Pm_neg_er = Pm_neg - Pm
            Pm_pos_er_disp = f"{Pm_pos_er:.1f}"
            Pm_neg_er_disp = f"{Pm_neg_er:.1f}"

            Rm = Pm * Re
            Rm_pos = Pm_pos * Re_pos
            Rm_neg = Pm_neg * Re_neg
            Rm_pos_er = Rm_pos - Rm
            Rm_neg_er = Rm_neg - Rm
            Rm_pos_er_disp = round(Rm_pos_er / 10**getPwr(Rm), 1)
            Rm_neg_er_disp = round(Rm_neg_er / 10**getPwr(Rm), 1)

            table1Data += f"{SOLVER_DICT[solver]} & ----- & ${p_bn:.2f}^{{+{p_bn_pos_er:.2f}}}_{{{p_bn_neg_er:.2f}}}$ & ${k_bn:.1f}^{{+{k_bn_pos_er:.1f}}}_{{{k_bn_neg_er:.1f}}}$ & ${k_nu:.2f}^{{+{k_nu_pos_er:.2f}}}_{{{k_nu_neg_er:.2f}}}$ &  ${k_eta:.1f}^{{+{k_eta_pos_er:.1f}}}_{{{k_eta_neg_er:.1f}}}$ \\\\ \n"
            table2Data += f"{SOLVER_DICT[solver]} & ${getNum(Re)}\\times10^{{{getPwr(Re)}}}$ & ${getNum(Rm)}\\times10^{{{getPwr(Rm)}}}$ "

        table1Str = f'''
        \\begin{{table*}}
        \\centering
        \\setlength{{\\tabcolsep}}{{1.8pt}}
        \\renewcommand{{\\arraystretch}}{{1.5}}
        \\begin{{tabular}}{{cccccccc}}
        \\hline
        Name & $\\tau$ & $p_\\mathrm{{bn}}$ & $k_\\mathrm{{bn}}$ & $\\tilde{{k}}_\\mathrm{{\\nu}}(=k_\\mathrm{{\\nu}})$ & $k_\\mathrm{{\\eta}}$\\\\
        (1) & (2) & (3) & (4) & (5) & (6)\\\\
        \\hline
        {table1Data}
        \\hline
        \\end{{tabular}}
        \\caption{{All parameters were measured/derived by averaging over the kinematic phase of the dynamo when $t>5t_\\mathrm{{turb}}$ and $10^{{-7}}\\le E_\\mathrm{{mag}}/E_\\mathrm{{kin}} \\le 10^{{-3}}$. Columns: \\textbf{{(1)}} Name of the numerical scheme as described in Table~\\ref{{tab:solvers}}. \\textbf{{(2)}} Exponent of the power law part of the kinetic spectra. \\textbf{{(3)}} Exponent of the bottleneck effect. \\textbf{{(4)}} Scaling wave number of the bottleneck effect. \\textbf{{(5)}} Viscous dissipation wave number. \\textbf{{(6)}} Resistive dissipation wave number.}}
        \\label{{tab:Turbulent dynamo fit parameters}}
        \\end{{table*}}
        '''
        table2Str = f'''
        \\begin{{table*}}
        \\centering
        \\setlength{{\\tabcolsep}}{{1.8pt}}
        \\renewcommand{{\\arraystretch}}{{1.5}}
        \\begin{{tabular}}{{ccc}}
        \\hline
        Name & $Re$ & $Rm$ \\\\
        (1) & (2) & (3) \\\\
        \\hline
        {table2Data}
        \\hline
        \\end{{tabular}}
        \\caption{{Columns: \\textbf{{(1)}} Name of the numerical scheme as described in Table~\\ref{{tab:solvers}}. \\textbf{{(2)}} Effective Hydrodynamic Reynolds number. \\textbf{{(3)}} Effective Magnetic Reynolds number.}}
        \\label{{tab:Turbulent dynamo effective flow numbers}}
        \\end{{table*}}
        '''
        table1File = open("Table1.txt", "w")
        table1File.write(table1Str)
        table1File.close()
        table2File = open("Table2.txt", "w")
        table2File.write(table2Str)
        table2File.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, nargs="*", help="Input Directories")
    parser.add_argument("-o", type=str, default="./", help="Output Directory")
    parser.add_argument("-kin_spect", type=int, default=0, help="Generate kinetic spectra")
    parser.add_argument("-cur_spect", type=int, default=0, help="Generate current spectra")
    parser.add_argument("-mag_spect", type=int, default=0, help="Generate magnetic spectra")
    parser.add_argument("-kin_plot", type=int, default=0, help="Plot kinetic spectra")
    parser.add_argument("-cur_plot", type=int, default=0, help="Plot current spectra")
    parser.add_argument("-mag_plot", type=int, default=0, help="Plot magnetic spectra")
    # parser.add_argument("-legend", type=int, default=0, help="Want legend? 0 will place legend on Kinetic spectra only")
    parser.add_argument("-table", type=int, default=0, help="Generate table of fit parameters")
    parser.add_argument("-v", type=int, default=0, help="Verbose")
    parser.add_argument("-lf", type=float, help="Lower bound for kinematic phase")
    parser.add_argument("-uf", type=float, help="Upper bound for kinematic phase")
    parser.add_argument("-stf", type=float, help="Start time for kinematic phase")
    parser.add_argument("-n", type=int, default=1, help="Number of processors for spectra generation")

    commonKeys = ["n", "o", "v", "table", "mag_spect", "kin_spect", "cur_spect", "mag_plot", "kin_plot", "cur_plot", "legend", "lf", "uf", "stf"]

    args = parseArgs(parser.parse_args(), commonKeys)
    print(args)

    main(args)