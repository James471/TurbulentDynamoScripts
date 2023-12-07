#!/usr/bin/env python3
import argparse
import sim
import plot
import common


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

        fig, axes = plot.main(common.Object(i=dirs[i], outdir=outdir, t=1, save=save, ylim_mag=ylim_mag, 
                                            ylim_ratio=ylim_ratio, no_adj_mag=no_adj_mag, 
                                            no_adj_ratio=no_adj_ratio, show=show_local, e=extraArgs[i], 
                                            title=title, fit=fit, fit_range=fit_range, skiprows=skiprows), 
                              fig=fig, axes=axes)


def createSim(simArgs):
    sim.main(simArgs)


def main(args):
    if args.mode == "dir":
        dirList = args.d_list

    elif args.mode == "sim":
        numSim = len(args.v)
        dirList = []
        for i in range(numSim):
            simArgs = common.Object(v=args.v[i], auto_adjust=args.auto_adjust[i], useVisc=args.useVisc[i], 
                                    Re=args.Re[i], useMgRes=args.useMgRes[i], Prm=args.Prm[i], L=args.L[i], 
                                    c=args.c[i], solver=args.solver[i], mcut=args.mcut[i], 
                                    outdir=args.outdir, sol_weight=args.sol_weight[i], cfl=args.cfl[i], 
                                    E_method=args.E_method[i], nt=args.nt, dt=args.dt, iprocs=args.iprocs[i], 
                                    jprocs=args.jprocs[i], kprocs=args.kprocs[i], nxb=args.nxb[i], 
                                    nyb=args.nyb[i], nzb=args.nzb[i], flash_path=args.flash_path, 
                                    extra=args.extra[i])
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

commonKeys = ["mode", "outdir", "title", "flash_path", "nt", "dt", "show", "ylim_mag", "ylim_ratio", "no_adj_mag", "no_adj_ratio", "fit", "fit_range", "skiprows"]

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
    simComp.add_argument("-useVisc", nargs="*", default=["false"], type=str, help="Use viscosity")
    simComp.add_argument("-Re", nargs="*", default=[1e3], type=float, help="Reynolds number")
    simComp.add_argument("-useMgRes", nargs="*", default=["false"], type=str, help="Use magnetic resistivity")
    simComp.add_argument("-Prm", nargs="*", default=[1], type=float, help="Prandtl number")
    simComp.add_argument("-L", nargs="*", default=[1], type=float, help="Length of the box")
    simComp.add_argument("-c", nargs="*", default=[1], type=float, help="Speed of sound")
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
        default=common.FLASH_PATH,
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
