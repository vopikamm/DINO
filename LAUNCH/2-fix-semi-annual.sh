#!/usr/bin/env bash
#SBATCH --job-name=2-fix-semi-annual-cycle-in-the-forcing
#SBATCH --hint=nomultithread       # 1 MPI process per physical core (no hyperthreading)
#SBATCH --time=02:00:00
#SBATCH --output=2-fix_semi_annual.out
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=david.kamm@locean.ipsl.fr
#SBATCH --account=omr@cpu # cpu accounting
#SBATCH --qos=qos_cpu-dev # QoS
######################

expdir='EXP00/'
expnam='2-fix-semi-annual-cycle-in-the-forcing'
cfgdir='DINO/'
cfgpath='/gpfswork/rech/omr/uym68qx/nemo_4.2.0/tests/'
savedir='/gpfsstore/rech/omr/uym68qx/nemo_output/'
# tasks_xios='1'
# tasks_nemo='16'

store=false

cd /gpfsscratch/rech/omr/uym68qx

if [[ ! -e $cfgdir ]]; then
    mkdir $cfgdir; cd $cfgdir
else
    cd $cfgdir
fi

if [[ ! -e $expnam ]]; then
    mkdir $expnam; cd $expnam
else
    cd $expnam
fi

pwd; hostname; date

module purge # purge modules inherited by default
# conda deactivate # deactivate environments inherited by default

# chargement des modules
intel_version=19.0.4
module load intel-compilers/${intel_version}
module load intel-mpi/${intel_version}
module load hdf5/1.10.5-mpi
module load netcdf/4.7.2-mpi
module load netcdf-fortran/4.5.2-mpi


cp ${cfgpath}${cfgdir}${expdir}/*.xml .
cp ${cfgpath}${cfgdir}${expdir}/namelist* ./
#cp ${cfgpath}${cfgdir}${expdir}/../BLD/bin/nemo.exe ./
cp ${cfgpath}${cfgdir}${expdir}/nemo ./
cp ${cfgpath}${cfgdir}${expdir}/mpmd.conf ./

source $I_MPI_ROOT/intel64/bin/mpivars.sh release_mt

time srun  --multi-prog ./mpmd.conf

# store results in the store directory
if [ "$store" = true ]; then
    cd ${savedir}

    if [[ ! -e $cfgdir ]]; then
        mkdir $cfgdir; cd $cfgdir
    else
        cd $cfgdir
    fi

    if [[ ! -e $expnam ]]; then
        mkdir $expnam
    fi

    cp *.nc ${savedir}${cfgdir}${expnam}
    cp output* ${savedir}${cfgdir}${expnam}
    cp communication_report.txt ${savedir}${cfgdir}${expnam}
    cp time.step ${savedir}${cfgdir}${expnam}
    cp run.stat ${savedir}${cfgdir}${expnam}
    cp ocean.output ${savedir}${cfgdir}${expnam}
fi

date

