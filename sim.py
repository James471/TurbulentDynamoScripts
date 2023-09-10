#!/usr/bin/env python3
import os
import argparse
import pickle
import textwrap
import common


def argsToSolverSetupParams(args):
    solverSetupParams = ""
    if args.solver == "bk-usm" or args.solver == "Roe" or args.solver == "HLLC":
        solverSetupParams += "+usm"
    elif args.solver == "bouchut-split":
        solverSetupParams += "+bouchut"
    return solverSetupParams


def getFlashParSolverParams(args):
    if args.solver == "bk-usm":
        return f"""
        RiemannSolver = "HLL5R"
        hy_hll5r_low_mach_cut = {args.mcut}
        """
    elif args.solver == "Roe":
        return f"""
        RiemannSolver = "Roe"
        """
    elif args.solver == "HLLC":
        return f"""
        RiemannSolver = "HLLC"
        """
    elif args.solver == "HLLD":
        return f"""
        RiemannSolver = "HLLD"
        """


def getFlashParMiscParams(args):
    def getEUpwind(args):
        return ".true." if args.E_upwind else ".false."

    return f"""
    E_upwind = {getEUpwind(args)}
    """


def createOutputDirectory(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    outputDirName = common.argsToOutdirName(args)

    if os.path.exists(outputDirName):
        deleteExistingDir = input(
            "Output directory"
            + outputDirName
            + " already exists. Do you want to delete it? (y/n): "
        )
        if deleteExistingDir == "y":
            os.system("rm -rf " + outputDirName)
        else:
            print("Exiting...")
            exit(0)
    os.makedirs(outputDirName)


def createSimulationObjectDirectory(args):
    solverSetupParams = argsToSolverSetupParams(args)
    objDirectory = common.argsToSimulationObjectDirectory(args)

    createSimulationCmd = (
        args.flash_path
        + "/setup StirFromFile -auto -objdir="
        + objDirectory
        + " -3d -nxb="
        + str(args.nxb)
        + " -nyb="
        + str(args.nyb)
        + " -nzb="
        + str(args.nzb)
        + " +ug "
        + solverSetupParams
        + " +stir_ics +polytropic_eos -parfile=flash.par.james"
    )
    os.system(createSimulationCmd)


def createInfoDumpFile(args):
    infoDict = {}

    infoDict["simulation"] = "Turbulent Dynamo"
    infoDict["dimensions"] = "3d"
    infoDict["v"] = args.v
    infoDict["auto_adjust"] = args.auto_adjust
    infoDict["solver"] = args.solver
    if args.solver == "bk-usm":
        infoDict["mcut"] = args.mcut
    infoDict["cfl"] = args.cfl
    infoDict["E_upwind"] = args.E_upwind
    infoDict["nt"] = args.nt
    infoDict["dt"] = args.dt
    infoDict["grid"] = "uniform"
    infoDict["nxb"] = args.nxb
    infoDict["nyb"] = args.nyb
    infoDict["nzb"] = args.nzb
    infoDict["iprocs"] = args.iprocs
    infoDict["jprocs"] = args.jprocs
    infoDict["kprocs"] = args.kprocs
    infoDict["extra"] = args.extra
    infoDict["turnover_time"] = 1 / (2 * args.v)

    infoFilePath = common.argsToOutdirName(args) + "/info.pkl"
    with open(infoFilePath, "wb") as infoFile:
        pickle.dump(infoDict, infoFile)


def createFlashPar(args):
    turnOverTime = 1 / (2 * args.v)
    tmax = args.nt * turnOverTime
    checkpointFileIntervalTime = turnOverTime
    plotFileIntervalTime = args.dt * turnOverTime

    flashParFile = open(common.argsToOutdirName(args) + "/flash.par", "w")

    flashParFileContent = f"""
    {getFlashParSolverParams(args)}
    {getFlashParMiscParams(args)}
    # kinematic viscosity and magnetic resistivity
    useViscosity                 = .false.
    visc_whichCoefficientIsConst = 2 # then set diff_visc_nu
    diff_visc_nu                 = 1e-3 # L*Mach*cs/(2*Re) implies Re=1000

    useMagneticResistivity       = .false.
    resistivity                  = 1e-3

    useConductivity              = .false.
    useSTS                       = .false.
    useSTSforDiffusion           = .false.

    dt_diff_factor  = 0.2

    # numerical viscosity switches
    use_avisc = .false.
    cvisc     = 0.1

    # use the turbulence stirring module (stirring from file)
    useStir              =  .true.
    st_infilename        =  "TurbGen.par"
    rho_ambient          =  1.0
    c_ambient            =  1.0
    magnetic             =  .true.
    MagField_z           =  1e-10
    st_computeDt         =  .false.

    st_rmsMagneticField  = 0
    st_MagneticSpectForm = 1
    st_stirMagneticKMin  = 1.0
    st_stirMagneticKMax  = 3.0

    st_MPMagneticFluxConservation = .true.
    st_MPxBmeanTarget = 0.0
    st_MPyBmeanTarget = 0.0
    st_MPzBmeanTarget = 0.0

    # simulation box
    xmin            = -0.5
    xmax            = +0.5
    ymin            = -0.5
    ymax            = +0.5
    zmin            = -0.5
    zmax            = +0.5

    xl_boundary_type      = "periodic"
    xr_boundary_type      = "periodic"
    yl_boundary_type      = "periodic"
    yr_boundary_type      = "periodic"
    zl_boundary_type      = "periodic"
    zr_boundary_type      = "periodic"

    # Tracer particles
    useParticles  = .false.
    pt_maxPerProc = 1000
    pt_numX       = 10
    pt_numY       = 10
    pt_numZ       = 10
    pt_initialXMin = -0.5
    pt_initialXMax = +0.5
    pt_initialYMin = -0.5
    pt_initialYMax = +0.5
    pt_initialZMin = -0.5
    pt_initialZMax = +0.5

    basenm          = "Turb_"
    log_file        = "flash.log"
    stats_file      = "flash.dat"
    run_comment     = "FLASH4"

    # file numbers - if you restart you have to change the checkpointFileNumber
    restart                 = .false.
    checkpointFileNumber	= 0
    plotFileNumber		= 0

    # Movie
    use_movie = .true.
    movie_dt_dump =	0.05

    # set the time between dumps
    checkpointFileIntervalTime  = {checkpointFileIntervalTime}
    plotFileIntervalTime        = {plotFileIntervalTime}
    particleFileIntervalTime    = {plotFileIntervalTime}
    dtmax                       = 0.025
    tmax                        = {tmax}

    # set the number of steps between dumps
    checkpointFileIntervalStep  = 0
    plotFileIntervalStep	    = 0

    wall_clock_time_limit = 160000.0
    wall_clock_checkpoint = 86000.0
    wr_integrals_freq = 1

    dtinit = 1.e-4
    dtmin  = 1.e-99  # This parameter must be << minimum timestep
                    #  in order to avoid numerical instability
    smallt = 1.e-99
    smalle = 1.e-99
    smlrho = 1.e-99

    plot_var_1      = "dens"
    plot_var_2      = "velx"
    plot_var_3      = "vely"
    plot_var_4      = "velz"
    plot_var_5      = "magx"
    plot_var_6      = "magy"
    plot_var_7      = "magz"

    gamma           = 1.0001
    eintSwitch      = 0
    cfl             = {args.cfl}
    nend            = 1000000

    # magnetic fields
    flux_correct          = .true.
    killdivb              = .true.
    UnitSystem            = "CGS"

    # polytropic EOS
    usePolytrope         = .true.
    IsothermalKonst      = 1.0 # (cs^2)
    polytropeKonst       = 1.0 # (cs^2)
    polytropeGamma1      = 1.0
    polytropeGamma2      = 1.0
    polytropeGamma3      = 1.0
    polytropeGamma4      = 1.0
    polytropeGamma5      = 1.0
    polytropeDens1       = 1e-99
    polytropeDens2       = 1e99
    polytropeDens3       = 1e99
    polytropeDens4       = 1e99
    polytropeDens5       = 1e99

    #   AMR refinement parameters
    # lrefine_min = 1
    # lrefine_max = 1
    # refine_var_1 = "dens"
    # nblockx = 2
    # nblocky = 2
    # nblockz = 2

    #These parameters below are only necessary for the Uniform Grid

    iProcs = {args.iprocs}      #num procs in i direction
    jProcs = {args.jprocs}      #num procs in j direction
    kProcs = {args.kprocs}


    # When using UG, iProcs, jProcs and kProcs must be specified.
    # These are the processors along each of the dimensions

    #FIXEDBLOCKSIZE mode ::
    # When using fixed blocksize, iGridSize etc are redundant in
    # runtime parameters. These quantities are calculated as 
    # iGridSize = NXB*iprocs
    # jGridSize = NYB*jprocs
    # kGridSize = NZB*kprocs

    #NONFIXEDBLOCKSIZE mode ::
    # iGridSize etc must be specified. They constitute the global
    # number of grid points in the physical domain without taking 
    # the guard cell into account. The local blocksize is calculated
    # as iGridSize/iprocs  etc.
    """

    flashParFile.write(textwrap.dedent(flashParFileContent).strip())
    flashParFile.close()


def createTurbGenPar(args):
    turbGenParFile = open(common.argsToOutdirName(args) + "/TurbGen.par", "w")

    turbGenParFileContent = f"""
    # ********************************************************************************
    # *** Input parameter file for controlling turbulence driving                  ***
    # *** Please see Federrath et al. (2010, A&A 512, A81) for details and cite :) ***
    # ********************************************************************************
    ndim               = 3              # N-dimensional turbulence driving (1 or 2 or 3).
                                        # Note that ndim = 1.5 or 2.5 will create 2 or 3 vector field components, respectively.
    L                  = 1.0            # Length of simulation domain (box) in [x[y[z]] (can be comma-separated list to set each component).
    velocity           = {args.v}            # Target turbulent velocity dispersion.
                                        # The following parameters (ampl_factor) is used to adjust the driving amplitude in [x[y[z]],
                                        # to approach the desired target velocity dispersion. Note that it scales as velocity/velocity_measured,
                                        # e.g., given a target velocity dispersion of 'velocity' and a measured velocity dispersion of 'velocity_measured',
                                        # scale the current value of ampl_factor by velocity/velocity_measured, such that
                                        # ampl_factor(new) = ampl_factor(previous) * velocity / velocity_measured. This will need adjustment, because
                                        # different codes and numerical resolutions/schemes will result in somewhat different dissipation.
                                        # Further, if the driving range and/or solenoidal weight (see below) are changed, these parameters will also need
                                        # to be adjusted, so the target turbulent velocity dispersion is reached in x[y[z]].
    ampl_factor        = 1.0            # Adjust [x[y[z]] amplitude (can be comma-separated list to set each component).
    ampl_auto_adjust   = {args.auto_adjust}              # Automatic amplitude adjustment switch (0: off, 1: on).
    k_driv             = 2.0            # Characteristic driving scale in units of 2pi / L[x].
                                        # L[x], k_driv, and velocity are used to set the Ornstein-Uhlenbeck auto-correlation time.
    k_min              = 1.0            # Minimum driving wavenumber in units of 2pi / L[x].
    k_max              = 3.0            # Maximum driving wavenumber in units of 2pi / L[x].
                                        # Note that while this is set based on L[x] only, the driving is still going to be isotropic,
                                        # even if L[x] != L[y] != L[z], because the actual modes are set based on L[x], L[y], L[z] during initialisation.
    sol_weight         = {args.sol_weight}            # 1.0: solenoidal driving, 0.0: compressive driving, 0.5: natural mixture.
    spect_form         = 1              # Spectral form of the driving amplitude. 0: band/rectangle/constant, 1: paraboloid, 2: power law.
    power_law_exp      = -2.0           # If spect_form = 2, this sets the spectral power-law exponent (e.g., -5/3: Kolmogorov; -2: Burgers)
    angles_exp         = 1.0            # If spect_form = 2, this sets the number of modes (angles) in k-shell surface,
                                        # such that it increases as k^angles_exp.
                                        # For full sampling, angles_exp = 2.0; for healpix-type sampling, angles_exp = 0.0.
    random_seed        = 140281         # Random number seed for driving sequence.
    nsteps_per_t_turb  = 10             # Number of turbulence driving pattern updates per turnover time.
    # ***************************************************************************************************
    """

    turbGenParFile.write(textwrap.dedent(turbGenParFileContent).strip())
    turbGenParFile.close()


def createFlashExec(args):
    currentPath = os.getcwd()
    os.chdir(common.argsToSimulationObjectDirectory(args))
    os.system("make -j 20")
    os.chdir(currentPath)


def createFlashSymLink(args):
    # If the target of a symbolic link is a relative path, it's interpreted relative to the
    # directory containing the link, not the directory you were in when you created it.
    os.system(
        "ln -sv objStirFromFile/flash4 " + common.argsToOutdirName(args) + "/flash4"
    )


def runSimulation(args):
    currentPath = os.getcwd()
    os.chdir(common.argsToOutdirName(args))
    os.system("mpirun -np " + str(args.iprocs * args.jprocs * args.kprocs) + " flash4")
    os.chdir(currentPath)


def main(args):
    print("Creating output directory")
    createOutputDirectory(args)
    print("Creating simulation object directory")
    createSimulationObjectDirectory(args)
    print("Creating info dump file")
    createInfoDumpFile(args)
    print("Creating flash.par")
    createFlashPar(args)
    print("Creating TurbGen.par")
    createTurbGenPar(args)
    print("Creating flash executable")
    createFlashExec(args)
    print("Creating symbolic link to flash executable")
    createFlashSymLink(args)
    print("Running simulation")
    runSimulation(args)


def parseArgs(args):
    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run a simulation of a turbulent dynamo."
    )
    parser.add_argument("-v", default=1.0, type=float, help="Velocity amplitude")
    parser.add_argument(
        "-auto_adjust",
        default=0,
        type=int,
        help="Automatically adjust the velocity amplitude",
    )
    parser.add_argument(
        "-solver",
        default="bk-usm",
        choices=["HLLC", "Roe", "bk-usm", "bouchut-split", "HLLD"],
        type=str,
        help="Solver to use",
    )
    parser.add_argument("-mcut", default=0.1, type=float, help="Low mach cutoff")
    parser.add_argument("-outdir", default=".", type=str, help="Output directory")
    parser.add_argument(
        "-sol_weight",
        default=0.5,
        type=str,
        help="1.0: solenoidal driving, 0.0: compressive driving, 0.5: natural mixture",
    )
    parser.add_argument("-cfl", default=0.5, type=float, help="CFL number")
    parser.add_argument("-E_upwind", default=True, type=bool, help="Use E_upwind")
    parser.add_argument("-nt", default=100, type=float, help="Number of turnover times")
    parser.add_argument(
        "-dt",
        default=0.1,
        type=float,
        help="Fraction of turnover time for plotfile dump",
    )
    parser.add_argument(
        "-iprocs", default=4, type=int, help="Number of processors in i-direction"
    )
    parser.add_argument(
        "-jprocs", default=4, type=int, help="Number of processors in j-direction"
    )
    parser.add_argument(
        "-kprocs", default=4, type=int, help="Number of processors in k-direction"
    )
    parser.add_argument(
        "-nxb", default=32, type=int, help="Number of blocks in i-direction"
    )
    parser.add_argument(
        "-nyb", default=32, type=int, help="Number of blocks in j-direction"
    )
    parser.add_argument(
        "-nzb", default=32, type=int, help="Number of blocks in k-direction"
    )
    parser.add_argument(
        "-flash_path",
        default="/home/james471/Academics/Projects/MHD/Code/flash-rsaa",
        help="Path to flash directory",
    )
    parser.add_argument(
        "-extra", type=str, help="Extra arguments to pass to the simulation"
    )

    args = parseArgs(parser.parse_args())

    main(args)
