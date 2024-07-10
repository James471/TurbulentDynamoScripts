## Turbulent Dynamo Scripts

### Running

**Snapshot functionality is currently broken**

In order to produce all the plots so far on gadi, do the following-

1. `cd` to the simulation directory located at `/scratch/ek9/jw5893/TurbulentDynamo/Mach0.01`.
2. Run `main.py` using `python3 ~/TurbulentDynamoScripts/main.py -i ./`.
Currently, this would do nothing since all fit parameters and spectra files are already present in the parent directory of the simulation. In order to reattempt the fits, use the flag `-redo`. In order to regenerate the spectra files, use `-recreate_spectra`. If you want to change the initial fit parameters passed for the kinetic fit, edit the file `simulationDirectory/spectra/kinFitInit.txt`.