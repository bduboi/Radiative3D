#!/bin/bash -ux

# Local testing config
BASE=`pwd`  # Base directory for the run.
PYTHON_EXEC="python3"
ENV_PATH="$BASE/R3Denv/bin/activate"

NPROC=2   # Number of processes for local test

SOURCE_ID=$(cat $BASE/Params/source_id.txt)
echo "Running local test for Source : $SOURCE_ID"

RUNID="spherical-moon-R3D-parallel-$SOURCE_ID"
tstamp=`date +"%Y%m%d-%H%M%S"`
outdirname="$tstamp"-"$RUNID"
outdir="$BASE/data/Parallel-runs/$outdirname"
mkdir -p "$outdir"
echo $outdir
echo "$outdir" > $BASE/data/Parallel-runs/latest_run.txt

# Local MPI run (no PBS, no machinefile)
mpirun -np $NPROC $BASE/scripts/do-spherical-python.sh


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

## Activate the Python environment:
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
echo "Plotting model maps and cross-sections."

cat process_0/stdout.txt | \
    grep -A 10000 "#  R3D_GRID:" | \
    grep -B 10000 "#  END R3D_GRID" | \
    #revtail -n +2 | \
    #tail -n +9 | \
    sed 's/\*\*\*/000/g' > griddump.txt  # Dump grid data to file


"$PYTHON_EXEC" -q - <<EOF
import os,sys
sys.path.append("$BASE")
from Radiative_python_funcs import *
current_path = os.getcwd()
GG=read_gridgeom(f'{current_path}/griddump.txt')
fig,ax=plt.subplots(1,1,figsize=(10,5))
plot_earth_layers(GG,ax,$RAD,colormap=None,zorder=1)
os.makedirs('Figures',exist_ok=True)
plt.savefig('Figures/gridGeometry.pdf',dpi=200)
EOF



#if [ "$1" != "noseis" ]; then
if [ "seis" != "noseis" ]; then
    echo "Plotting seismograms..."
    # Print the directory
    "$PYTHON_EXEC" -q - <<EOF
import os, sys
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
sys.path.append("$BASE/Python/")
from Radiative_python_funcs import *

# Define the output folder
output_folder = 'Figures'
os.makedirs(output_folder, exist_ok=True)
listdir = os.listdir('$PWD/seisfiles')
n_files = len(listdir)
# Loop over the files and process them
for i in range(0,n_files,2):
    file =sprintf(i) + '.pkl'
    file_path = os.path.join('$PWD/seisfiles', file)
    ofile = f'{file.replace(".pkl", ".pdf")}'
    print(f"Processing {ofile}...")
    # Read and plot the seismogram
    seisplot(file_path)    
    
    # Save the figure
    plt.savefig(os.path.join(output_folder, ofile), dpi=200)
    plt.close()

EOF
fi

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
