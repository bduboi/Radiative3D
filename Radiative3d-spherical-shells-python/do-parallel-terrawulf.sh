#!/bin/bash -ux
#
#PBS -N Radiative3D
#PBS -M Balthazar.Dubois@anu.edu.au
#PBS -d /home/balthazar/Radiative3d-spherical-shells/data/log/
#PBS -o stdout.log
#PBS -e stderr.log
#PBS -l nodes=10:t3:ppn=12
#PBS -l walltime=999:00:00


# Configure the run name and output directory
ENV_PATH="$BASE/R3Denv/bin/activate"
BASE="/home/balthazar/Radiative3d-spherical-shells"
PYTHON_EXEC="/home/balthazar/R3Denv/bin/python"
PARAMFILE="$BASE/source_params.txt"
export OMPI_LOGDIR=$PBS_JOBID
echo "Running with $PBS_NP total processes."

SOURCE_ID=$(cat $BASE/source_id.txt)
echo "Running on Simulation for Source : $SOURCE_ID"

RUNID="spherical-moon-R3D-M2-tt-$SOURCE_ID"
tstamp=`date +"%Y%m%d-%H%M%S"`   # Time stamp for the run
outdirname="$tstamp"-"$RUNID"
outdir="$BASE/data/Parallel-runs/$outdirname"
mkdir -p "$outdir"
echo $outdir
# Save the directory name in a file for identification when running on other cores, possibly later after queuing.
echo "$outdir" > $BASE/data/Parallel-runs/latest_run.txt

mpirun -mca btl "^openib" -machinefile "$PBS_NODEFILE" --map-by ppr:12:node /usr/local/bin/logdir $BASE/do-moon-ttcurve-python.sh

rm -rf /home/balthazar/Radiative3d-spherical-shells/data/log
mkdir -p /home/balthazar/Radiative3d-spherical-shells/data/log


## Figure processing for the BATCH RUN (paralllel run on terrawulf)

# go to base directory
cd $outdir

# Copy the do-figsonly.sh script to the output directory
cp /home/balthazar/Radiative3d-spherical-shells/scripts/do-figsonly-python.sh $outdir

# Give the script execute permissions
chmod +x $outdir/do-figsonly-python.sh


#Get a copy of the stoud file
cp process_0/stdout.txt $outdir/stdout.txt

## Activate the Python environment:
source /home/balthazar/R3Denv/bin/activate


# Assemble the seifiles into a single metadata file
seispattern="process_0/seisfiles/seis_???_asc.dat" # Here using process_0 as the directory to identify the file pattern
echo "Assembling seismograms..."
for file1 in "$outdir"/$seispattern;do
    SeisName=$(basename $file1)
    echo "SeisName: $SeisName"
    "$PYTHON_EXEC" -q - <<EOF
import os, sys
sys.path.append("$BASE")
from Radiative_python_funcs import *
merge_traces(Directory="$outdir", SeisName="$SeisName", NCore=$PBS_NP)
EOF
done
# Make seismometer and model maps:
echo "Plotting model maps and cross-sections."

cat process_0/stdout.txt | \
    grep -A 10000 "#  R3D_GRID:" | \
    grep -B 10000 "#  END R3D_GRID" | \
    #revtail -n +2 | \
    #tail -n +9 | \
    sed 's/\*\*\*/000/g' > griddump.txt  # Dump grid data to file


"$PYTHON_EXEC" -q - <<EOF
import os,sys
sys.path.append("$BASE/Python/")
from Radiative_python_funcs import *
current_path = os.getcwd()
GG=read_gridgeom(f'{current_path}/griddump.txt')
fig,ax=plt.subplots(1,1,figsize=(10,5))
plot_earth_layers(GG,ax,colormap=None,zorder=1)
os.makedirs('Figures',exist_ok=True)
plt.savefig('Figures/gridGeometry.pdf',dpi=200)
EOF

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
#echo End Figure Generation: `date +"%Y.%m.%d-%H:%M:%S"`# >> "$LOGFILE"
echo '///'
echo "Figure generation complete. To re-run figures, run 'do-figsonly.sh' in"
echo "the output directory, which may be modified to customize figure output."
echo "Output has been placed in $outdir."
#echo "Logfile contents:"
#cat "$outdir"/logfile


# Copy the files to the output 
## END
##
