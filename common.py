def argsToOutdirName(args):
    outputDirName = (
        args.outdir
        + "/Turb_v"
        + str(args.v)
        + "_auto-adj"
        + str(args.auto_adjust)
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
