#!/bin/bash

BASE="/home/balthazar/Radiative3d-spherical-shells"
PYTHON_EXEC="/home/balthazar/R3Denv/bin/python"


outdir=`pwd`

PBS_NP=$(ls -d process_* 2>/dev/null | wc -l) # Count the number of process directories

## Activate the Python environment:
source /home/balthazar/R3Denv/bin/activate

: '
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
'
# Make seismometer and model maps:
echo "Plotting model maps and cross-sections."

cat process_0/stdout.txt | \
    grep -A 10000 "#  R3D_GRID:" | \
    grep -B 10000 "#  END R3D_GRID" | \
    #revtail -n +2 | \
    #tail -n +9 | \
    sed 's/\*\*\*/000/g' > griddump.txt  # Dump grid data to file

: '
"$PYTHON_EXEC" -q - <<EOF
import os,sys
sys.path.append("$BASE")
from Radiative_python_funcs import *
current_path = os.getcwd()
GG=read_gridgeom(f'{current_path}/griddump.txt')
fig,ax=plt.subplots(1,1,figsize=(10,5))
plot_earth_layers(GG,ax,colormap=None,zorder=1)
os.makedirs('Figures',exist_ok=True)
plt.savefig('Figures/gridGeometry.pdf',dpi=200)
EOF
'

if [ -z "$BASE" ]; then
    echo "The BASE environment variable (path to current directory) is not set. Exiting."
    exit 1
fi

#if [ "$1" != "noseis" ]; then
if [ "seis" != "noseis" ]; then
    echo "Plotting seismograms..."
    # Print the directory
    "$PYTHON_EXEC" -q - <<EOF
import os, sys
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
sys.path.append("$BASE/")
from Radiative_python_funcs import *

# Define the output folder
output_folder = 'Figures'
os.makedirs(output_folder, exist_ok=True)

# Loop over the files and process them
for i in range(1,160,100):

    file =sprintf(i) + '_asc.dat.pkl'
    #file = f'seis_00{i}_asc.dat.pkl'
    file_path = os.path.join('$PWD/seisfiles', file)
    ofile = f'{file.replace("_asc.dat.pkl", ".pdf")}'
    print(f"Processing {ofile}...")
    # Read and plot the seismogram
    with open(file_path, 'rb') as f:
        metadata = pickle.load(f)
    seisplot(metadata)
    
    # Save the figure
    plt.savefig(os.path.join(output_folder, ofile), dpi=200)
    plt.close()

EOF
fi


### To do : Set up the normcurve comparison
# Individual ttcurves

produce_ttcurves() {  # $1: station code
                     # $2: start index
                     # $3: end index
                     # $4: caxis limit
  "$PYTHON_EXEC" -q - <<EOF
import os, sys
sys.path.append("$BASE")
from Radiative_python_funcs import *
current_path = os.getcwd()
station_code = "$1"  # Pass station code as a string
array = assembleArray(f'{current_path}/seisfiles', $2, $3,fsignature='asc.dat.pkl')
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
plot_ttcurve(array, ax, 0, station_code, gamma=2, norm=0.3, theoretical=True,basedir='$BASE')
plt.savefig(f'Figures/traveltime-{station_code}-xyz-individual.pdf', dpi=200)
plt.close()
EOF
}

# TT Curves:
STA="STA"
STB="STB"
STC="STC"

echo "Plotting individual tt curves..."
produce_ttcurves $STA 0 159    # (Out the P-wave principle axis)
#produce_ttcurves $STB 160 319    # (Out the S-wave principle axis)
#produce_ttcurves $STC 320 479    # (Out the equivalency axis)
## (Do not edit/remove this comment block.)
## ___FIG_GEN_END___
##
cd "$rwd"
echo End Figure Generation: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"
echo '///'
echo "Figure generation complete. To re-run figures, run 'do-figsonly.sh' in"
echo "the output directory, which may be modified to customize figure output."
echo "Output has been placed in $outdir."
echo "Logfile contents:"
cat "$outdir"/logfile

# Copy the files to the output 
## END
##

exit 0
