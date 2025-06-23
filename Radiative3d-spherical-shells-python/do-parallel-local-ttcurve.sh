#!/bin/bash -ux

# Local testing config
BASE="/Users/balthazar/Downloads/Radiative3d-spherical-shells"
PYTHON_EXEC="python3"
ENV_PATH="$BASE/R3Denv/bin/activate"
PARAMFILE="$BASE/Params/source_params.txt"
NPROC=4   # Number of processes for local test

SOURCE_ID=$(cat $BASE/Params/source_id.txt)
echo "Running local test for Source : $SOURCE_ID"

RUNID="spherical-moon-R3D-parallel-ttcurve-$SOURCE_ID"
tstamp=`date +"%Y%m%d-%H%M%S"`
outdirname="$tstamp"-"$RUNID"
outdir="$BASE/data/Parallel-runs/$outdirname"
mkdir -p "$outdir"
echo $outdir
echo "$outdir" > $BASE/data/Parallel-runs/latest_run.txt

# Local MPI run (no PBS, no machinefile)
mpirun -np $NPROC $BASE/scripts/do-ttcurve-python.sh

# Time 10 secs to allow the run to complete
rm -rf $BASE/data/log
mkdir -p $BASE/data/log


## Figure processing for the BATCH RUN (paralllel run on terrawulf)

# go to base directory
cd $outdir

# Copy the do-figsonly.sh script to the output directory
cp $BASE/scripts/do-figsonly-python.sh $outdir

# Give the script execute permissions
chmod +x $outdir/do-figsonly-python.sh


#Get a copy of the stoud file
cp process_0/stdout.txt $outdir/stdout.txt

source $ENV_PATH

# Assemble the seifiles into a single metadata file
seispattern="process_0/seisfiles/seis_???.pkl" # Here using process_0 as the directory to identify the file pattern
echo "Assembling seismograms..."
for file1 in "$outdir"/$seispattern;do
    SeisName=$(basename $file1)
    echo "SeisName: $SeisName"
    "$PYTHON_EXEC" -q - <<EOF
import os, sys
sys.path.append("$BASE/Python/")
from Radiative_python_funcs import *
merge_traces(Directory="$outdir", SeisName="$SeisName", NCore=$NPROC)
EOF
done

# Make seismometer and model maps:
echo "Plotting the tt curve."

### To do : Set up the normcurve comparison
RAD=$(cat $BASE/Params/rad.txt)
# Individual Traveltime Curves:

produce_ttcurves() {  # $1: station code
                     # $2: start index
                     # $3: end index
                     # $4 gamma normalisation (optional, default 0.5)
  "$PYTHON_EXEC" -q - <<EOF
import os, sys
sys.path.append("$BASE/Python/")
from Radiative_python_funcs import *
current_path = os.getcwd()
station_code = "$1"  # Pass station code as a string

ibegin = $2 # Start index for the array
iend = $3 # End index for the array
array_dict = assembleArray("$outdir/seisfiles", ibegin, iend)
print(array_dict['BinSize'])
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
plot_ttcurve("$BASE",array_dict, ax, 0, station_code, gamma=$4, norm=0.3, theoretical=False,R=$RAD)
os.makedirs('Figures', exist_ok=True)
plt.savefig(f'Figures/traveltime-{station_code}-0-$3.pdf', dpi=200)
plt.close()
EOF
}

# TT Curves:
STA="STA"

# Get the maximum number of seismogram from the seisfiles
maxseis=$(($(ls $outdir/seisfiles/seis_???.pkl | wc -l) - 1))
echo "Maximum number of seismograms: $maxseis"
produce_ttcurves $STA 0 $maxseis 0.1   # (Out the P-wave principle axis)

## (Do not edit/remove this comment block.)
## ___FIG_GEN_END___
##
cd "$outdir"
echo '///'
echo "Figure generation complete. To re-run figures, run 'do-figsonly.sh' in"
echo "the output directory, which may be modified to customize figure output."
echo "Output has been placed in $outdir."

# Copy the files to the output 
## END
##
