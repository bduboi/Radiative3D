#!/bin/bash -ux
#PBS -N Radiative3D
#PBS -M Balthazar.Dubois@anu.edu.au
#PBS -m abe
#PBS -d /home/balthazar/Radiative/data/log
#PBS -o stdout.log
#PBS -e stderr.log
#PBS -l nodes=1:t3:ppn=4
#PBS -l walltime=999:00:00

# Load environment
BASE="/home/balthazar/Radiative"
pwd=`pwd`
echo $BASE > $pwd/BASE.txt

PYTHON_EXEC="$BASE/R3Denv/bin/python"
ENV_PATH="$BASE/R3Denv/bin/activate"
SCRIPT="$BASE/scripts/do-spherical-python.sh"
NPROC=$PBS_NP

# Output run config
SOURCE_ID=$(cat $BASE/Params/source_id.txt)
RUNID="spherical-moon-R3D-parallel-$SOURCE_ID"
tstamp=$(date +"%Y%m%d-%H%M%S")
outdirname="$tstamp"-"$RUNID"
outdir="$BASE/data/Parallel-runs/$outdirname"

mkdir -p "$outdir"
echo "$outdir" > "$BASE/data/Parallel-runs/latest_run.txt"
echo "Running on $NPROC cores. Source: $SOURCE_ID"
echo "Output directory: $outdir"

# Run main simulation
mpirun -np $PBS_NP -mca btl ^openib --map-by ppr:12:node --machinefile "$PBS_NODEFILE" $SCRIPT

# Prepare logs
#rm -rf "$BASE/data/log"
#mkdir -p "$BASE/data/log"

# ---- Postprocessing & Figures ----
cd "$outdir"

# Copy figure script
cp "$BASE/scripts/do-figsonly-python.sh" "$outdir"
chmod +x "$outdir/do-figsonly-python.sh"

# Copy stdout
cp process_0/stdout.txt "$outdir/stdout.txt"

# Activate Python
source "$ENV_PATH"

# Assemble seisfiles
seispattern="process_0/seisfiles/seis_???.pkl"
echo "Assembling seismograms..."
for file1 in "$outdir"/$seispattern; do
    SeisName=$(basename "$file1")
    echo "SeisName: $SeisName"
    "$PYTHON_EXEC" -q - <<EOF
import os, sys
sys.path.append("$BASE/Python/")
from Radiative_python_funcs import *
merge_traces(Directory="$outdir", SeisName="$SeisName", NCore=$NPROC)
EOF
done

echo "Plotting model maps and cross-sections."
RAD=$(cat $BASE/Params/rad.txt)
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
