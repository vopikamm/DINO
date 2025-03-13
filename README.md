# DINO [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15016824.svg)](https://doi.org/10.5281/zenodo.15016824)

A DIabatic Neverworld Ocean in NEMO 4.2.1. This repository contains the source code for the NEMO test configuration `DINO`. An archived version together with the reference experiments for 1°, 1/4° and 1/16° horizontal resolution are published at [zenodo](https://doi.org/10.5281/zenodo.15016824). A short summary of the configuration details can be found in [this Python module](https://github.com/vopikamm/dinostics/blob/7f776460128afb2af9153aab06af70e4aad152b7/dino_configuration.py) with some visualizations in [this jupyter notebook](https://github.com/vopikamm/dinostics/blob/7f776460128afb2af9153aab06af70e4aad152b7/DINO_config.ipynb). A more elaborate description of the configuration is currently in review.

## Content
* /EXPREF:  1°    spin-up from scratch (no forcing, restart files needed)
* /MY_SRC:  source files
* /cpp_DINO.fcm   keys for compilation       

## Installation
Please follow the [NEMO installation guide](https://sites.nemo-ocean.io/user-guide/install.html#essential-components) and install XIOS and NEMO [version 4.2.1](https://forge.nemo-ocean.eu/nemo/nemo/-/releases/4.2.1) (!) as described. If you are new to NEMO I suggest following the guide until you manage to run your own [test configuration](https://sites.nemo-ocean.io/user-guide/install.html#running-the-model) before continuing with DINO.

### Downloading necessary code/data

Go into the `tests` folder of your NEMO installation and clone the `DINO_4.2.1` branch of the repository:

```
cd <YOUR_NEMO_INSTALLATION_PATH>/tests
git clone -b DINO_4.2.1 https://github.com/vopikamm/DINO.git
```

### Compilation 

First you need to add `DINO` to the configuration list in `demo_cfgs.txt`:

```
echo "DINO  OCE" >> ../demo_cfgs.txt
```
If you chose to rename the configuration folder, replace `"DINO  OCE"` with `"<YOUR_NAME>  OCE"` and also rename the file [cpp_DINO.fcm](https://github.com/vopikamm/DINO/blob/4bbbe8337de59952eb2f42a20b2fd30250354dee/cpp_DINO.fcm) to `cpp_<YOUR_NAME>.fcm`. 
Now compile from the NEMO root directory, e.g.:

``` 
../makenemo -n 'DINO_compiled' -a 'DINO' -m '<YOUR_ARCH_FILE>' -j 32
```
This will compile DINO into a new folder `/DINO_compiled` add another new folder `/EXP00` inside, which contains a symbolic link to the excecutable `nemo.exe`. This experiment corresponds to running DINO at 1° from scratch.

### Build your experiment
If you plan to reproduce the reference experiments, you will need the restart files, evaporation-minus-precipitation (emp) fields, file definition files and namelists found in `/reference_experiments.zip` of the [**zenodo release**](https://doi.org/10.5281/zenodo.15016824). They are unfortunately too large to be provided via Github. Download and copy the folder into your compiled DINO configuration `/DINO_compiled`. Now you are ready to reproduce the 1° (R1), 1/4° (R4) and 1/16° (R16) experiments, both from the `initial state` (init) or just from the last few years of `data production` (prod), as indicated in the spin-up strategy:

![png](readme_figures/schematic.png)

Choose the experiment you want to run and if you want to run it from the `init` or `prod` restart file. Copy the forcing, restart, file definition and namelist files into the `/EXP00` folder (**overwriting** the namelist_cfg and the file_def_nemo-oce.xml!), e.g. for R1 from the production restart:

```
cd ./DINO_compiled/EXP00
cp -f ../Reference_experiments/EXP_R1/namelist_cfg_prod ./namelist_cfg
cp -f ../Reference_experiments/EXP_R1/file_def_nemo-oce.xml ./file_def_nemo-oce.xml
cp ../Reference_experiments/EXP_R1/restart_prod_1deg.nc ./
cp ../Reference_experiments/EXP_R1/emp_1deg.nc ./
```

### Run the experiment:
As for the [NEMO installation guide](https://sites.nemo-ocean.io/user-guide/install.html#running-the-model) you can run the experiment (Here in detached mode on 37 CPUs):

```
mpirun -np 36 ./nemo : -np 1 <PATH_TO_YOUR>/xios_server.exe
```

## Additional notes:
* NEMO stores the runtime in seconds. Therefore you cannot run DINO for more than 2.^31 / 24. /3600 ./ 360. = 69.04 to store such large numbers with integer 4. For this reason and not obtaining too large output files we suggest to run DINO in batches of 50 years (R1), 10 years (R4) and 1 year (R16) restarts.
* The computational ressources vary from machine to machine and depend on the level of parallelisation. But for reference, here are some rough estimates from my experiments on the Jean Zay HPC:
    - R1:  on 36 + 1   CPUs  ~ 2    hCPU / simulated year
    - R4:  on 585 + 15  CPUs  ~ 120  hCPU / simulated year
    - R16: on 4368 + 112 CPUs  ~ 9000 hCPU / simulated year
* All above was written for NEMO/XIOS in detached mode, where one CPU per computing node was reserved for XIOS. This is advised for efficient input/output reading/writing.
