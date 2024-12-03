#!/usr/bin/env python3
import sys
import os
import glob
import socket
import configparser

import matplotlib
matplotlib.use('Agg')

import numpy as np  
import matplotlib.pyplot as pl
import pandas as pd
from scipy.special import kn
from scipy.optimize import fminbound
from scipy.interpolate import CubicSpline
import emcee

from constants import *
from datautils import *
from utils import *
import designParams
import flashlib as fl
import turblib as tl

from ipdb import set_trace as stop

p_kin = -1.7 # global hack

def Re_fun(k_nu, c_nu):
    return (k_nu/(c_nu * 2))**(4/3)


def k_eta_fun(k_tilde_eta, p_eta):
    return (k_tilde_eta)**(1/p_eta)


def Rm_fun(k_nu, c_nu, k_tilde_eta, p_eta, c_eta):
    k_eta = k_eta_fun(k_tilde_eta, p_eta)
    Pm = (k_eta/(c_eta*k_nu))**2
    return Pm * Re_fun(k_nu, c_nu)


def Log10_P_kin(k, A_kin, p_bn, k_bn, k_nu_tilde, p_nu):
    return np.log10(((k / k_bn)**p_kin + (k / k_bn)**p_bn) * np.exp(-(k / k_nu_tilde)**p_nu)) + A_kin


def Log10_P_mag(k, A_mag, p_mag, p_eta, k_tilde_eta):
    y = A_mag * (k**p_mag) * kn(0, (k / k_tilde_eta)**p_eta)
    return np.log10(y)


def generateSpectra(simDir, verbose, spectType, infoDict, lf, uf, stf, nProcs=1):
    '''
    verbose  : Could be 0 or 1 or 2
    spectType: Could be 'vels' or 'curr' or 'mags' or 'vort'
    '''

    if verbose == 2:
        datFileVerbose = 1
        devNull = ""
    else:
        datFileVerbose = 0
        devNull = " > /dev/null 2>&1"

    spectraDir = simDir + "/spectra"
    if not os.path.isdir(spectraDir):
        print("creating spectra dir: '"+spectraDir+"'")
        os.mkdir(spectraDir)

    files = sorted(glob.glob(simDir+"/Turb_hdf5_plt_cnt_????"))

    t_turb = 1 / (2 * infoDict['v'])
    if not os.path.exists(simDir+"/Turb.dat_cleaned"):
        if verbose > 0: print("Loading turb.dat for cleaning")
        datfile_obj = fl.datfile(simDir+"/Turb.dat", verbose=datFileVerbose)
        datfile_obj.write_cleaned()
        if verbose > 0: print("Turb.dat cleaned")
    else:
        print("Turb.dat_cleaned already exists")

    n = getNForTurbDat(simDir, res=5e-2, filename="Turb.dat_cleaned")
    if verbose > 0: print("Loading cleaned turb.dat")
    dat = loadFile(simDir + "/Turb.dat_cleaned", shift=10, n=n)
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
            print("In current spectra")
            if "nid" in socket.gethostname():
                derivCmd = f"srun -n {nProcs} /software/projects/pawsey0810/jwatt/flash-tools/tools/derivative_var/derivative_var {f} -current {devNull}"
            else:
                derivCmd = f"mpirun -np {nProcs} /home/100/jw5893/flash-tools/tools/derivative_var/derivative_var {f} -current {devNull}"
            os.system(derivCmd)
    elif spectType == "vort":
        types = "0 -dsets vorticity_x vorticity_y vorticity_z"
        for x in range(dumpNStart, dumpNEnd+1):
            if verbose > 0: print("Adding vorticity variable to dump:", x)
            f = simDir+"/Turb_hdf5_plt_cnt_{:04d}".format(x)
            if "nid" in socket.gethostname():
                derivCmd = f"srun -n {nProcs} /software/projects/pawsey0810/jwatt/flash-tools/tools/derivative_var/derivative_var {f} -vort {devNull}"
            else:
                derivCmd = f"mpirun -np {nProcs} /home/100/jw5893/flash-tools/tools/derivative_var/derivative_var {f} -vort {devNull}"
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
    elif spectType == "vort":
        spectFileString = "dset_vorticity_x_vorticity_y_vorticity_z"

    for x in range(dumpNStart, dumpNEnd+1):
        if verbose > 0: print("Processing dump to generate spectra:", x)
        f = simDir+"/Turb_hdf5_plt_cnt_{:04d}".format(x)
        if "nid" in socket.gethostname():
            spectrCmd = f"srun -n {nProcs} /software/projects/pawsey0810/jwatt/flash-tools/tools/spectra_mpi/spectra_mpi {f} -types {types} {devNull}"
        else:
            spectrCmd = f"mpirun -np {nProcs} /home/100/jw5893/flash-tools/tools/spectra_mpi/spectra_mpi {f} -types {types} {devNull}"
        print(f"Running: {spectrCmd}")
        os.system(spectrCmd)
    mvCmd = "mv "+simDir+"/Turb_hdf5_plt_cnt_????_spect_*.dat "+spectraDir
    os.system(mvCmd)

    spectFiles = [simDir+"/spectra/Turb_hdf5_plt_cnt_{:04d}".format(x)+"_spect_"+spectFileString+".dat" for x in range(dumpNStart, dumpNEnd+1)]
    if spectType == "vels": normalise = False
    if spectType == "mags" or spectType == "cur" or spectType == "vort": normalise = True # normalise magnetic spectra and current spectra before averaging, because of field growth with time
    averDat, headerAver = tl.aver_spect(spectFiles, normalise=normalise, verbose=0)
    outfile = simDir+"/spectra/aver_spect_"+spectFileString+".dat"
    tl.write_spect(outfile, averDat, headerAver, verbose=verbose)


def plot_with_error(ax, x, y, yerr, *args, **kwargs):
    ax.plot(x, y, *args, **kwargs)
    ax.fill_between(x, y-np.abs(yerr[0]), y+np.abs(yerr[1]), color=ax.lines[-1].get_color(), alpha=0.5, linewidth=0.0)


def plotSpectra(ax, simDir, verbose, spectType, fact, infoDict, params, outdir, compensate=False, fit=True, color="black", label=""):
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
    elif spectType == "vort":
        spectVar = "dset_vorticity_x_vorticity_y_vorticity_z"

    averFile = simDir + f"/spectra/aver_spect_{spectVar}.dat"
    averDf   = pd.read_csv(averFile, sep="\s+", header=0, skiprows=5)
    
    log10PTot = averDf['log10_Ptot']
    deltaLog10PTot = averDf['sigma_log10_Ptot']
    pTot = 10**log10PTot
    k = averDf['k']
    err = np.array([pTot * np.log(10) * deltaLog10PTot, pTot * np.log(10) * deltaLog10PTot])
    
    kFitMin = 4
    kFitMax = max(k) / 4
    kFit = np.array(k[(k >= kFitMin) & (k <= kFitMax)].tolist())
    log10PTotFit = np.array(log10PTot[(k >= kFitMin) & (k <= kFitMax)].tolist())
    deltaLog10PTotFit = np.array(deltaLog10PTot[(k >= kFitMin) & (k <= kFitMax)].tolist())

    if compensate and spectType == "vels":
        compensateFact = k**(-p_kin)
        compensateFitFact = kFit**(-p_kin)
        plot_with_error(ax, k, compensateFact * fact * pTot, fact * err * np.array([compensateFact, compensateFact]), label=label, color=color)
    else:
        compensateFact = 1
        compensateFitFact = 1
        plot_with_error(ax, k, fact * pTot, fact * err, label=label, color=color)

    if fit:
        fitDict = {}
        fitDict = fit_func(ax, spectType, simDir, kFit, log10PTotFit, deltaLog10PTotFit, params, compensateFitFact, fact, fitDict, verbose, infoDict)
    else:
        if spectType == "mags":
            fitDict = loadDict(f"{outdir}/magFitDict.pkl")[infoDict['solver']]
        elif spectType == "vels":
            fitDict = loadDict(f"{outdir}/kinFitDict.pkl")[infoDict['solver']]
            ax.plot(kFit, compensateFitFact * fact * 10**Log10_P_kin(kFit, fitDict["A_kin"][0], fitDict["p_bn"][0], fitDict["k_bn"][0], fitDict["k_nu_tilde"][0], fitDict["p_nu"][0]), color="black")
        elif spectType == "cur":
            fitDict = loadDict(f"{outdir}/curFitDict.pkl")[infoDict['solver']]
        elif spectType == "vort":
            fitDict = loadDict(f"{outdir}/vortFitDict.pkl")[infoDict['solver']]

    return fitDict


def fit_func(ax, spectType, simDir, kFit, log10PTotFit, deltaLog10PTotFit, params, compensateFitFact, fact, fitDict, verbose, infoDict):

    if spectType == "mags":
        fitDict["A_mag"] = (0,0,0)
        fitDict["p_mag"] = (0,0,0)
        fitDict["p_eta"] = (0,0,0)
        fitDict["k_tilde_eta"] = (0,0,0)

    elif spectType == "vort":
        fitDict["A_vort"] = (0,0,0)
        fitDict["p_vort"] = (0,0,0)
        fitDict["p_eta"] = (0,0,0)
        fitDict["k_tilde_eta"] = (0,0,0)

    elif spectType == "cur":
        if verbose: print("Current Spectra-> Working with:", infoDict['solver'])
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

        k_eta_23 = np.percentile(kList, 50)
        k_eta_23_pos_er = np.percentile(kList, 84) - k_eta_23
        k_eta_23_neg_er = np.percentile(kList, 16) - k_eta_23
        fitDict["k_eta_23"] = (k_eta_23, k_eta_23_neg_er, k_eta_23_pos_er)
        print(f"{infoDict['solver']}: {k_eta_23:.2f}^{{{k_eta_23_pos_er:.2f}}}_{{{k_eta_23_neg_er:.2f}}}")

    elif spectType == "vels":
        if verbose: print("Kin Spectra-> Fitting for:", infoDict['solver'])
        def log_likelihood(theta, x, y, yerr):
            A_kin, p_bn, k_bn, k_nu_tilde, p_nu = theta
            model = Log10_P_kin(x, A_kin, p_bn, k_bn, k_nu_tilde, p_nu)
            sigma2 = yerr**2
            return -0.5 * np.sum((y - model) ** 2 / sigma2)
        A_kin_min, A_kin_guess, A_kin_max = params["A_kin"]
        p_bn_min, p_bn_guess, p_bn_max = params["p_bn"]
        k_bn_min, k_bn_guess, k_bn_max = params["k_bn"]
        k_nu_tilde_min, k_nu_tilde_guess, k_nu_tilde_max = params["k_nu_tilde"]
        p_nu_min, p_nu_guess, p_nu_max = params["p_nu"]
        log_prior = lambda theta: 0.0 if A_kin_min < theta[0] < A_kin_max and p_bn_min < theta[1] < p_bn_max and k_bn_min < theta[2] < k_bn_max and k_nu_tilde_min < theta[3] < k_nu_tilde_max and p_nu_min < theta[4] < p_nu_max else -np.inf
        def log_probability(theta, x, y, yerr):
            lp = log_prior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + log_likelihood(theta, x, y, yerr)
        pos = np.array([[A_kin_guess, p_bn_guess, k_bn_guess, k_nu_tilde_guess, p_nu_guess] + 1e-4 * np.random.randn(5) for i in range(32)])
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(kFit, log10PTotFit, deltaLog10PTotFit))
        sampler.run_mcmc(pos, 10000, progress=True)
        labels = ["A_kin", "p_bn", "k_bn", "k_nu_tilde", "p_nu"]
        flat_samples = sampler.get_chain(discard=500, thin=50, flat=True)
        params = []
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            qrt = np.diff(mcmc)
            if compensateFitFact == 1:
                print(f"{labels[i]} = {mcmc[1]:.3f}_{{-{qrt[0]:.3f}}}^{{{qrt[1]:.3f}}}")
            fitDict[labels[i]] = (mcmc[1], -qrt[0], qrt[1])
            params.append(mcmc[1])
        ax.plot(kFit, fact * 10**Log10_P_kin(kFit, *(params)) * compensateFitFact, color="black")
    
    return fitDict


def postPlot(ax, spectType, showx=False, compensated=False):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ylabel=""
    xlabel=""
    if spectType == "mags":
        ylabel=r'$P_\mathrm{mag}$'
        ax.set_ylabel(ylabel)
        if not showx:
            ax.set_xticklabels([])
        # ax.set_xlabel(r'$k$')
        # ax.legend(loc='best')
    elif spectType == "vort":
        ylabel=r'$P_\mathrm{vort}$'
        ax.set_ylabel(ylabel)
        if not showx:
            ax.set_xticklabels([])
    elif spectType == "vels":
        if compensated:
            ylabel=r'$k^{('+str(-p_kin)+')}\,P_\mathrm{kin}$'
            ylabel=r'$k^{('+str(-p_kin)+')}\,P_\mathrm{kin}$'
        else:
            ylabel=r'$P_\mathrm{kin}$'
        if not showx:
            ax.set_xticklabels([])
        # ax.legend(loc='best')
    elif spectType == "cur":
        ylabel=r'$P_\mathrm{cur}$'
        xlabel=r'$k$'
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        # ax.legend(loc='best')
    return xlabel, ylabel

def plotScaleLoc(ax, solverFit, type, color_dict, printvals=True):
    ylim = ax.get_ylim()
    if type == "mags":
        ax.plot([10, 10], [ylim[0], 2*ylim[0]], color="white", scaley=False)
        ax.text(10, 2.5*ylim[0], r"$k_\nu$", color="white")
        # maxKEta = 0
        # for solver in solverFit:
        #     val = solverFit[solver]
        #     ax.plot([val["k_tilde_eta"][0]**(1/val["p_eta"][0]), val["k_tilde_eta"][0]**(1/val["p_eta"][0])], [ylim[0], 2*ylim[0]], color=COLOR_DICT[solver], scaley=False)
        #     if val["k_tilde_eta"][0]**(1/val["p_eta"][0]) > maxKEta:
        #         maxKEta = val["k_tilde_eta"][0]**(1/val["p_eta"][0])
        # ax.text(maxKEta, 2.5*ylim[0], r"$k_\eta$", color="black")
    elif type == "vort":
        ax.plot([10, 10], [ylim[0], 2*ylim[0]], color="white", scaley=False)
        ax.text(10, 2.5*ylim[0], r"$k_\nu$", color="white")
    elif type == "vels":
        maxKNu = 0
        for solver in solverFit:
            val = solverFit[solver]
            k_nu = val["k_nu_tilde"][0]**(1/val["p_nu"][0])
            k_nu = val["k_nu_tilde"][0]**(1/val["p_nu"][0])
            if printvals:
                print(f"Value of k_nu for {solver}: {k_nu}")
            ax.plot([k_nu, k_nu], [ylim[0], 2*ylim[0]], color=color_dict[solver], scaley=False)
            if k_nu > maxKNu:
                maxKNu = k_nu
        ax.text(maxKNu, 2.5*ylim[0], r"$k_\nu$", color="black")
    elif type == "cur":
        maxKEta = 0
        for solver in solverFit:
            val = solverFit[solver]
            ax.plot([val["k_eta_23"][0], val["k_eta_23"][0]], [ylim[0], 2*ylim[0]], color=color_dict[solver], scaley=False)
            if val["k_eta_23"][0] > maxKEta:
                maxKEta = val["k_eta_23"][0]
        ax.text(maxKEta, 2.5*ylim[0], r"$k_\eta$", color="black")


def distribute_parameter(heading, sub_heading, config, input_dirs, cast=str):
    param = [cast(i.strip()) for i in config[heading][sub_heading].split(",")]
    if len(param) == 1:  # Single value supplied
        return [param[0]] * len(input_dirs)
    elif len(param) == len(input_dirs):  # Match the length of input_dirs
        return param
    else:
        raise ValueError("The number of parameters does not match the number of directories.")


def main(config):
    base = config['common']['base']
    output_dir = config['common']['output_dir']
    kin_spect = config['common'].getboolean('kin_spect')
    cur_spect = config['common'].getboolean('cur_spect')
    mag_spect = config['common'].getboolean('mag_spect')
    vort_spect = config['common'].getboolean('vort_spect')
    kin_plot = config['common'].getboolean('kin_plot')
    cur_plot = config['common'].getboolean('cur_plot')
    mag_plot = config['common'].getboolean('mag_plot')
    vort_plot = config['common'].getboolean('vort_plot')
    verbose = config['common'].getint('verbose')
    lf = config['common'].getfloat('lf')
    uf = config['common'].getfloat('uf')
    stf = config['common'].getfloat('stf')
    n = config['common'].getint('n')
    refit = config['common'].getboolean('refit')
    no_shift = config['common'].getboolean('no_shift')
    showx = config['common'].getboolean('showx')
    oname = config['common']['oname']

    sim_list = [base+"/"+i.strip() for i in config['directories']['i'].split(",")]

    fit_file_list = distribute_parameter('args', 'fitFile', config, sim_list)
    factor_list = distribute_parameter('args', 'factor', config, sim_list, cast=float)
    color_list = distribute_parameter('args', 'color', config, sim_list)
    label_list = distribute_parameter('args', 'label', config, sim_list)

    color_dict = {}
    for index, sim in enumerate(sim_list):
        color_dict[getInfoDict(sim)['solver']] = color_list[index]

    if kin_spect:
        for sim in sim_list:
            generateSpectra(sim, verbose, "vels", getInfoDict(sim), lf, uf, stf, n)
    if cur_spect:
        for sim in sim_list:
            generateSpectra(sim, verbose, "cur", getInfoDict(sim), lf, uf, stf, n)
    if mag_spect:
        for sim in sim_list:
            generateSpectra(sim, verbose, "mags", getInfoDict(sim), lf, uf, stf, n)
    if vort_spect:
        for sim in sim_list:
            generateSpectra(sim, verbose, "vort", getInfoDict(sim), lf, uf, stf, n)


    if kin_plot:
        solverKinFit = {}
        fig, ax_kin_obj = pl.subplots()
        fig2, ax_kin_obj_comp = pl.subplots()
        for index, sim_dir in enumerate(sim_list):
            infoDict = getInfoDict(sim_dir)
            fit_file_path = fit_file_list[index]
            if no_shift:
                fact = 1
            else:
                fact = factor_list[index]
            fit = True
            if not refit and os.path.exists(f"{output_dir}/kinFitDict.pkl"):
                fit = False
            kinParams = {"A_kin": [-8, -5, -2], "p_bn": [0, 1, np.inf], "k_bn": [0.1, 4.0, 128], "k_nu_tilde": [0.1, 4.0, 128], "p_nu": [1, 1, 1+1e-6]}
            if os.path.exists(fit_file_path):
                print("Found fit file at ", fit_file_path)
                kinParams = txtToCfpDict(fit_file_path)
            elif os.path.exists(sim_dir+"/spectra/kinFitInit.txt"):
                print("Found kinFitInit.txt")
                kinParams = txtToCfpDict(sim_dir+"/spectra/kinFitInit.txt")
            else:
                print("Using default kin dict")
            print(f"Using kinParams: {kinParams}")
            fitDict = plotSpectra(ax_kin_obj, sim_dir, verbose, "vels", fact, infoDict, kinParams, output_dir, compensate=False, fit=fit, color=color_list[index], label=label_list[index])
            _ = plotSpectra(ax_kin_obj_comp, sim_dir, verbose, "vels", fact, infoDict, kinParams, output_dir, compensate=True, fit=False, color=color_list[index], label=label_list[index])
            solverKinFit[infoDict['solver']] = fitDict
        xlabel, ylabel = postPlot(ax_kin_obj, "vels", showx, False)
        xlabel_comp, ylabel_comp = postPlot(ax_kin_obj_comp, "vels", showx, True)
        dumpDict(solverKinFit, f"{output_dir}/kinFitDict.pkl")
        plotScaleLoc(ax_kin_obj, solverKinFit, "vels", color_dict=color_dict)
        plotScaleLoc(ax_kin_obj_comp, solverKinFit, "vels", color_dict=color_dict, printvals=False)
        ax_kin_obj.set_xlabel(xlabel)
        ax_kin_obj.set_ylabel(ylabel)
        ax_kin_obj_comp.set_xlabel(xlabel_comp)
        ax_kin_obj_comp.set_ylabel(ylabel_comp)
        ax_kin_obj_comp.legend(loc="best")
        print("Saving Kinetic Spectra")
        ax_kin_obj.figure.savefig(f"{output_dir}/Kinetic_Spectra{oname}.pdf")
        ax_kin_obj_comp.figure.savefig(f"{output_dir}/Compensated_Kinetic_Spectra{oname}.pdf")
        # plKinObj.ax().figure.clf(); plKinObj.ax().cla(); pl.close(); plKinObj = None


    if mag_plot:
        solverMagFit = {}
        fig, ax_mag_obj = pl.subplots()
        for index, sim_dir in enumerate(sim_list):
            infoDict = getInfoDict(sim_dir)
            if no_shift:
                fact = 1
            else:
                fact = factor_list[index]
            fit = True
            if not refit and os.path.exists(f"{output_dir}/magFitDict.pkl"):
                fit = False
            # Mag fit has been disabled for now
            magParams = {"A_mag": [0, 0.0001, np.inf], "p_mag": [0, 1, np.inf], "p_eta": [0, 1, np.inf], "k_tilde_eta": [0, 4.0, 128]}
            fitDict = plotSpectra(ax_mag_obj, sim_dir, 1, "mags", fact, infoDict, magParams, output_dir, compensate=False, fit=fit, color=color_list[index], label=label_list[index])
            solverMagFit[infoDict['solver']] = fitDict
        xlabel, ylabel = postPlot(ax_mag_obj, "mags", showx)
        dumpDict(solverMagFit, f"{output_dir}/magFitDict.pkl")
        plotScaleLoc(ax_mag_obj, solverMagFit, "mags", color_dict=color_dict)
        ax_mag_obj.set_xlabel(xlabel)
        ax_mag_obj.set_ylabel(ylabel)
        print("Saving Magnetic Spectra")
        ax_mag_obj.figure.savefig(f"{output_dir}/Magnetic_Spectra{oname}.pdf")


    if cur_plot:
        solverCurFit = {}
        fig, ax_cur_obj = pl.subplots()
        for index, sim_dir in enumerate(sim_list):
            infoDict = getInfoDict(sim_dir)
            if no_shift:
                fact = 1
            else:
                fact = factor_list[index]
            fit = True
            if not refit and os.path.exists(f"{output_dir}/curFitDict.pkl"):
                fit = False
            fitDict = plotSpectra(ax_cur_obj, sim_dir, 1, "cur", fact, infoDict, None, output_dir, compensate=False, fit=fit, color=color_list[index], label=label_list[index])
            solverCurFit[infoDict['solver']] = fitDict
        xlabel, ylabel = postPlot(ax_cur_obj, "cur", showx)
        dumpDict(solverCurFit, f"{output_dir}/curFitDict.pkl")
        plotScaleLoc(ax_cur_obj, solverCurFit, "cur", color_dict=color_dict)
        ax_cur_obj.set_xlabel(xlabel)
        ax_cur_obj.set_ylabel(ylabel)
        print("Saving Current Spectra")
        ax_cur_obj.figure.savefig(f"{output_dir}/Current_Spectra{oname}.pdf")

    
    if vort_plot:
        solverVortFit = {}
        fig, ax_vort_obj = pl.subplots()
        for index, sim_dir in enumerate(sim_list):
            infoDict = getInfoDict(sim_dir)
            if no_shift:
                fact = 1
            else:
                fact = factor_list[index]
            fit = True
            if not refit and os.path.exists(f"{output_dir}/vortFitDict.pkl"):
                fit = False
            vortParams = {"A_vort": [0, 0.0001, np.inf], "p_vort": [0, 1, np.inf], "p_eta": [0, 1, np.inf], "k_tilde_eta": [0, 4.0, 128]}
            fitDict = plotSpectra(ax_vort_obj, sim_dir, 1, "vort", fact, infoDict, vortParams, output_dir, compensate=False, fit=fit, color=color_list[index], label=label_list[index])
            solverVortFit[infoDict['solver']] = fitDict
        xlabel, ylabel = postPlot(ax_vort_obj, "vort", showx)
        dumpDict(solverVortFit, f"{output_dir}/vortFitDict.pkl")
        plotScaleLoc(ax_vort_obj, solverVortFit, "vort", color_dict=color_dict)
        ax_vort_obj.set_xlabel(xlabel)
        ax_vort_obj.set_ylabel(ylabel)
        print("Saving Vorticity Spectra")
        ax_vort_obj.figure.savefig(f"{output_dir}/Vorticity_Spectra{oname}.pdf")


if __name__ == "__main__":
    config = configparser.ConfigParser()
    if len(sys.argv) > 1:
        config.read(sys.argv[1])
    else:
        config.read("spectra.ini")

    main(config)