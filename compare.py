#!/usr/bin/env python3
import argparse
import sim
import plot
import common


def getArgsForDirPlot(dir, save, extra, title, outdir, ylim_mag, 
                      ylim_ratio, no_adj_mag, no_adj_ratio, show, fit, fit_range, skiprows):
    
    Object = lambda **kwargs: type("Object", (), kwargs)
    return Object(i=dir, outdir=outdir, t=1, save=save, ylim_mag=ylim_mag, 
                  ylim_ratio=ylim_ratio, no_adj_mag=no_adj_mag, 
                  no_adj_ratio=no_adj_ratio, show=show, e=extra, title=title,
                  fit=fit, fit_range=fit_range, skiprows=skiprows)


def createComparisionPlot(dirs, extraArgs, title, outdir, ylim_mag=None, 
                          ylim_ratio=None, no_adj_mag=False, no_adj_ratio=False, show=False, 
                          fit=False, fit_range=None, skiprows=0):
    
    for i in range(len(dirs)):
        if i == 0:
            save, fig, axes = False, None, None
            show_local = False
        elif i == len(dirs) - 1:
            save = True
            show_local = show
        else:
            save = False
            show_local = False

        fig, axes = plot.main(getArgsForDirPlot(dirs[i], save, extraArgs[i], title, outdir, 
                                                ylim_mag, ylim_ratio, no_adj_mag, no_adj_ratio, show_local,
                                                fit, fit_range, skiprows), fig, axes)


def createSim(simArgs):
    sim.main(simArgs)


def getArgsForSim(v, auto_adjust, solver, mcut, outdir, sol_weight, 
                  cfl, E_method, nt, dt, iprocs, jprocs, kprocs, nxb, nyb, nzb, flash_path, extra):
    
    Object = lambda **kwargs: type("Object", (), kwargs)
    return Object(v=v, auto_adjust=auto_adjust, solver=solver, mcut=mcut, outdir=outdir, 
                  sol_weight=sol_weight, cfl=cfl, E_method=E_method, nt=nt, dt=dt, 
                  iprocs=iprocs, jprocs=jprocs, kprocs=kprocs, nxb=nxb, nyb=nyb, nzb=nzb, 
                  flash_path=flash_path, extra=extra)


def main(args):
    if args.mode == "dir":
        dirList = args.d_list

    elif args.mode == "sim":
        numSim = len(args.v)
        dirList = []
        for i in range(numSim):
            simArgs = getArgsForSim(args.v[i], args.auto_adjust[i], args.solver[i], args.mcut[i], 
                                    args.outdir, args.sol_weight[i], args.cfl[i], args.E_method[i], 
                                    args.nt, args.dt, args.iprocs[i], args.jprocs[i], args.kprocs[i], 
                                    args.nxb[i], args.nyb[i], args.nzb[i], args.flash_path, args.extra[i])
            dirList.append(common.argsToOutdirName(simArgs))
            createSim(simArgs)
    
    createComparisionPlot(dirList, args.extra, args.title, args.outdir, 
                          args.ylim_mag, args.ylim_ratio, args.no_adj_mag, args.no_adj_ratio, args.show,
                          args.fit, args.fit_range, args.skiprows)


def checkArgs(args):
    if args.mode == "dir":
        if args.extra is not None and len(args.extra) != 1:
            assert len(args.extra) == len(args.d_list), "Extra arguments must be specified for each directory"

    elif args.mode == "sim":
        numSims = -1
        for key, value in vars(args).items():
            if key not in commonKeys:
                print(key)
                if value is not None and len(value) > 1:
                    if numSims == -1:
                        numSims = len(value)
                    else:
                        assert numSims == len(value), "All arguments must have the same number of values unless you have common arguments"


def parseArgs(args):
    checkArgs(args)

    if args.mode == "dir":
        if args.extra is None:
            args.extra = [None for i in range(len(args.d_list))]
        elif len(args.extra) == 1:
            args.extra = [args.extra[0] for i in range(len(args.d_list))]
    elif args.mode == "sim":
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

commonKeys = ["mode", "outdir", "title", "flash_path", "nt", "dt", "show", "ylim_mag", "ylim_ratio", "no_adj_mag", "no_adj_ratio"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to automate comparision of FLASH Turbulent Dynamo simulations"
    )
    subparsers = parser.add_subparsers(dest="mode")

    dirComp = subparsers.add_parser("dir", help="Compare two directories")
    dirComp.add_argument("-d_list", nargs="*", type=str, help="Directories")
    dirComp.add_argument("-extra", nargs="*", type=str, help="Extra arguments")

    simComp = subparsers.add_parser("sim", help="Compare two simulations")
    simComp.add_argument(
        "-v", nargs="*", default=[0.1], type=float, help="Velocity scales"
    )
    simComp.add_argument(
        "-auto_adjust",
        nargs="*",
        default=[0],
        type=int,
        help="Automatically adjust the velocity amplitude",
    )
    simComp.add_argument(
        "-solver",
        nargs="*",
        default=["bk-usm"],
        choices=["HLLC", "Roe", "bk-usm", "bouchut-split", "HLLD"],
        type=str,
        help="Solvers",
    )
    simComp.add_argument("-mcut", nargs="*", default=[0.1], type=float, help="Mass cut")
    simComp.add_argument(
        "-sol_weight",
        nargs="*",
        default=[1.0],
        type=str,
        help="1.0: solenoidal driving, 0.0: compressive driving, 0.5: natural mixture",
    )
    simComp.add_argument("-cfl", nargs="*", default=[0.5], type=float, help="CFL")
    simComp.add_argument(
        "-E_method", nargs="*", type=str, default=["GS"], help="Method for calculation of electric field"
    )
    simComp.add_argument(
        "-iprocs",
        nargs="*",
        type=int,
        default=[2],
        help="Number of processors in i-direction",
    )
    simComp.add_argument(
        "-jprocs",
        nargs="*",
        type=int,
        default=[2],
        help="Number of processors in j-direction",
    )
    simComp.add_argument(
        "-kprocs",
        nargs="*",
        type=int,
        default=[1],
        help="Number of processors in k-direction",
    )
    simComp.add_argument(
        "-nxb",
        nargs="*",
        type=int,
        default=[32],
        help="Number of blocks in i-direction",
    )
    simComp.add_argument(
        "-nyb",
        nargs="*",
        type=int,
        default=[32],
        help="Number of blocks in j-direction",
    )
    simComp.add_argument(
        "-nzb",
        nargs="*",
        type=int,
        default=[32],
        help="Number of blocks in k-direction",
    )
    simComp.add_argument("-extra", nargs="*", type=str, help="Extra arguments")

    simComp.add_argument(
        "-nt", type=float, default=100, help="Number of turnover times"
    )
    simComp.add_argument(
        "-dt",
        type=float,
        default=0.1,
        help="Fraction of turnover time for creating plot files",
    )
    simComp.add_argument(
        "-flash_path",
        default="/home/james471/Academics/Projects/MHD/Code/flash-rsaa",
        type=str,
        help="Path to flash directory",
    )

    parser.add_argument("-outdir", type=str, default="./", help="Output directory")
    parser.add_argument("-title", type=str, help="Title of the plot")
    parser.add_argument(
        "-ylim_mag", type=float, nargs=2, help="Y-axis limits for E_mag plot"
    )
    parser.add_argument(
        "-ylim_ratio", type=float, nargs=2, help="Y-axis limits for E_mag/ERot plot"
    )
    parser.add_argument(
        "-no_adj_ratio", action="store_true", help="Don't adjust Ekin/Emag axis"
    )
    parser.add_argument(
        "-no_adj_mag", action="store_true", help="Don't adjust Emag axis"
    )
    parser.add_argument("-show", action="store_true", help="Show the figure")
    parser.add_argument("-fit", action="store_true", help="Fit the data to a power law")
    parser.add_argument("-fit_range", type=float, nargs=2, help="Range of data to fit")
    parser.add_argument("-skiprows", type=int, default=0, help="Number of rows to skip in the data file")

    args = parseArgs(parser.parse_args())

    main(args)
