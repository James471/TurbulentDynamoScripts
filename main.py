import argparse
import utils
import pathlib
import os


def main(args):

    scriptDir = pathlib.Path(__file__).parent.resolve()

    pathList = [args.i+_ for _ in args.i]
    simList = utils.getSolverSortedList(pathList)
    sims = simList.join(" ")

    if not args.no_growth:
        refit = ""
        if args.redo:
            refit = "-refit"
        cmd = f"python3 {scriptDir}/growth.py -i {sims} -o {args.o} -lf {args.lf} -uf {args.uf} -stf {args.stf} -sr {args.sr} -sbs {args.sbs} -ld {args.ld} -ud {args.ud} -fit {args.fit} {refit}"
        print("Generating growth plots")
        print("Running:", cmd)
        os.system(cmd)

    if not args.no_spectra:

        for sim in simList:
            if not os.path.exists(f"{sim}/spectra") or args.recreate_spectra:
                recreate = f"-kin_spect -mag_spect -cur_spect -n {args.n}"
                cmd = f"python3 {scriptDir}/spectra.py -i {sim} -o {args.o} {recreate}"
                print("Generating spectra files")
                print("Running:", cmd)
                os.system(cmd)

        shift = ""
        if args.no_shift:
            shift = "-no_shift"
        refit = ""
        if args.redo:
            refit = "-refit"
        cmd = f"python3 {scriptDir}/spectra.py -i {sim} -o {args.o} -kin_plot -mag_plot -cur_plot -lf {args.lf} -uf {args.uf} -stf {args.stf} {shift} {refit}"
        print("Generating spectra plots")
        print("Running:", cmd)
        os.system(cmd)

    if not args.no_snapshot:
        pass

    if not args.no_table:
        cmd = f"python3 {scriptDir}/table.py -i {args.i}"
        print("Generating tables")
        print("Running:", cmd)
        os.system(cmd)

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate production grade plots")
    parser.add_argument("-i", type=str, help="Input Directory containing all solver simulations")
    parser.add_argument("-o", type=str, default="./", help="Output Directory")
    parser.add_argument("-redo", action="store_true", help="Redo the analysis even if previous results are present")
    parser.add_argument("-no_growth", action="store_true", help="Do not plot growth rate")
    parser.add_argument("-no_spectra", action="store_true", help="Do not plot spectra")
    parser.add_argument("-no_snapshot", action="store_true", help="Do not plot snapshots")
    parser.add_argument("-no_table", action="store_true", help="Do not create tables")
    parser.add_argument("-lf", type=float, default=1e-7, help="Lower bound for kinematic phase")
    parser.add_argument("-uf", type=float, default=1e-3, help="Upper bound for kinematic phase")
    parser.add_argument("-stf", type=float, default=5.0, help="Start time for kinematic phase")
    parser.add_argument("-no_shift", action="store_true", help="Do not shift the spectra for different solvers")    

    # Growth plot arguments
    parser.add_argument("-sr", type=int, nargs="*", default=[12000], help="Skip rows")
    parser.add_argument("-sbs", type=float, default=[60.0], help="Skip turnover time before measuring saturation")
    parser.add_argument("-ld", type=float, nargs="?", default=1e-8, help="Low bound to display")
    parser.add_argument("-ud", type=float, default=5e0, help="Upper bound to display")
    parser.add_argument("-fit", type=str, default="syst", help="Fit method to use")

    # Spectra plot arguments
    parser.add_argument("-recreate_spectra", action="store_true", help="Recreate the spectra files")
    parser.add_argument("-n", type=int, default=4, help="Number of processors for spectra generation")

    # Snapshot plot arguments

    # Table arguments


    args = parser.parse_args()

    main(args)