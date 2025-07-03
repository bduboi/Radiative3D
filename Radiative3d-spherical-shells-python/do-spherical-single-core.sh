#!/bin/bash
##
##  "DO"-script for a Radiative3D run:
##
##  Copy and edit this file to easily specify (and remember)
##  parameters used for a particular run.
##
##  If $1 is set, it is used as an extra identifier token in the
##  output directory name.
##
##  (Note: Do not edit line 3.  It is used by subsequent scripts to
##  identify this file as a "do script". (See wikify.sh.))



##
##  This DO-SCRIPT sets up a run based on a spherical user chosen model file.
##

## One-liner description: (Keep this BRIEF.)
##


BASE=`pwd`  # Base directory for the run.
PYTHON_EXEC="python3"
ENV_PATH="$BASE/R3Denv/bin/activate"

### Initialization of the simulation name and output directory:
RUNID="spherical-moon-single-core"
tstamp=`date +"%Y%m%d-%H%M%S"`   # Time stamp for the run
outdirname="$tstamp"-"$RUNID"
outdir="$BASE/data/Single-core-runs/$outdirname"
mkdir -p "$outdir"
echo "$outdir" > $BASE/data/Single-core-runs/latest_run.txt
# Here precise that we are running a single core simulation.
SINGLE_CORE=1
echo "$SINGLE_CORE" > $BASE/Params/single_core.txt

source $BASE/scripts/do-fundamentals.sh


INTENT="Spherical whole-Earth simulation."
CAMPAIGN="" ## (e.g. "For AGU Poster Dec 2016" or some such.)

SIMTARGET="waveform"             # Choice: 'waveform' or 'video'. Affects
                              # defaults not otherwise specified.
#MODIDX=16                   # Spherical Earth  # Model Index: Selects from custom coded models.
MODIDX=666                    # 666 : Use a model file.
                              # 1: Lop Nor (base or moho depends on COMPARGS),
                              # 5: North Sea Crust Pinch model,
                              # 8: Crust Upthrust model
                              # 16: HARD CODED Spherical Earth (PREM)


#SOURCETYP=$(EventParams "Selby" 0)     # Lookup source parameters
SOURCETYP=SDR,90,90,0                 # 
                              # Source param choices include 'expl', 'eq',
                              # 'Selby', etc.  Second paramter is usually
                              # isofrac (if range [-1.0,1.0]) or isoangle
                              # (if range [-90,90]).
FREQ=2.0                      # Phonon frequency to model
NUMPHONS=1M                   # Number of phonons to spray (Recommend: 50M)
RECTIME=3600                  # Recording duration of seismometers (seconds).
BINSIZE=0.5                  # Seismometer time-bin size in seconds
GATHER=500.0                  # Terminal gather radius, in kilometers.
#MODELPATH="/Users/balthazar/Downloads/Radiative3d-spherical-shells/Models/MoonModels/SimplifiedModels/SimplfiedISSI_MOON_M2.csv" #Path of the model to be used. Default is SimplifiedPrem.csv in Models/EarthModels/PREM
MODELPATH="$BASE/Models/EarthModels/PREM/Simplified.csv" #Path of the model to be used. Default is SimplifiedPrem.csv in Models/EarthModels/PREM
GRIDSOURCE=GRID_FROMFILE      # Specify that the grid is read from a file. Otherwise use GRID_UNSPEC.


# /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
# THE FOLLOWING FUNCTION WILL DIRECTLY PARSE THE SEISMIC STATION BASED ON THE SOURCE CHOSEN    

SOURCE_ID=A01 # Chose a DMQ nest ID to run the simulation with. This will automatically initialize the seismometers.
                          #nStation, Stationfile1, SourceID1, Stationfile2, SourceID2, ... 

# Path to txt files with source parameters.
AP12="$BASE/Params/Stations/AP12.txt"
AP14="$BASE/Params/Stations/AP14.txt"
AP15="$BASE/Params/Stations/AP15.txt"
AP16="$BASE/Params/Stations/AP16.txt"

# MAKE SURE THAT THE GATHER RADIUS IS INITIALISED 
getSourceParams 4 $SOURCE_ID $AP12 $AP14 $AP15 $AP16 # Get the source parameters from the files.
source "$BASE/Params/source_params.txt"

#FLATTEN="--flatten"          # If set to "--flatten", apply Earth-flattening
FLATTEN=""                    # transformation to depths and velocities.
                              # (Set to "" to disable.)

RAD=1737.1
ADDITIONAL+="--earthrad=$RAD"  # Additional params. (Such as --ocsraw,
                              # or --earthrad=xxxx) Radii: Earth: 6371 km,
                              # Mars: 3389, Mercury: 2440, Earth's
                              # moon: 1737



## RUN THE SIMULATION:
##
##  Set up output directory and run the sim:
##

PopDefaults $SIMTARGET        ## Defined in do-fundamentals.sh
#CheckBuildStatus              ##  ''
CreateOutputDirectory $@       ##  ''
RunSimulation                  ##  ''

#echo Begin Figure Generation: >> "$LOGFILE"

rwd=`pwd`       # Switch to output directory - the rest of our work
outdir=$(cat $BASE/data/Single-core-runs/latest_run.txt)
cd "$outdir"

source "$ENV_PATH"  # Activate the Python environment

## (Do not edit/remove this comment block.)
#!/bin/bash
## ___FIG_GEN_START___
##

##                                      # Note: this part of the do-script will
## GENERATE FIGURES AND GRAPHICS:       #   be carved off into 'do-figsonly.sh'
##                                      #   to facilitate re-running and/or
                                        #   customizing figure output after the
                                        #   initial run.

# Organize seisfiles
mkdir -p seisfiles
mv $outdir/seis_???* seisfiles/ 2>/dev/null


#Convert the seisfiles to pkl format
echo "Converting seismograms to pkl format..."
"$PYTHON_EXEC" -q - <<EOF
import os, sys
sys.path.append("$BASE/Python/")
from Radiative_python_funcs import *
current_path = os.getcwd()
seisfiles = os.listdir(os.path.join(current_path, 'seisfiles'))
for file in seisfiles:
    if file.endswith('_asc.dat'):
        file_path = os.path.join(current_path, 'seisfiles', file)
        dat_to_pickle(file_path)

EOF
echo "Seismograms converted to pkl format."        


# Text values for Wiki script:
cat > wikify.incl <<EOF
  STA="STA"
  STB="STB"
  STC="STC"
  MODELCODE="SPHERE"
EOF


# Make seismometer and model maps:
echo "Plotting model maps and cross-sections."
cat stdout.txt | \
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


if [ -z "$BASE" ]; then
    echo "The BASE environment variable (path to current directory) is not set. Exiting."
    exit 1
fi

if [ "$1" != "noseis" ]; then
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

cd "$rwd"
echo End Figure Generation: `date +"%Y.%m.%d-%H:%M:%S"`
echo '///'
echo "Figure generation complete. To re-run figures, run 'do-figsonly.sh' in"
echo "the output directory, which may be modified to customize figure output."
echo "Output has been placed in $outdir."
echo "Logfile contents:"
cat "$outdir"/logfile
## END
##
