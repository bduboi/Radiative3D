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

SINGLE_CORE=0
echo "$SINGLE_CORE" > $BASE/single_core.txt

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
SOURCETYP=SDR,90,90,0                 # Or specify directly
#SOURCELOC=0,0,-10                       # Source location
                              # Source param choices include 'expl', 'eq',
                              # 'Selby', etc.  Second paramter is usually
                              # isofrac (if range [-1.0,1.0]) or isoangle
                              # (if range [-90,90]).

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=100K                   # Number of phonons to spray (Recommend: 50M)
RECTIME=3600                  # Recording duration of seismometers (seconds).
BINSIZE=0.5                  # Seismometer time-bin size in seconds
GATHER=500.0                  # Terminal gather radius, in kilometers.
MODELPATH="/Users/balthazar/Downloads/Radiative3d-spherical-shells/Models/MoonModels/SimplifiedModels/SimplfiedISSI_MOON_M2.csv" #Path of the model to be used. Default is SimplifiedPrem.csv in Models/EarthModels/PREM
GRIDSOURCE=GRID_FROMFILE      # Specify that the grid is read from a file. Otherwise use GRID_UNSPEC.

# /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
# THE FOLLOWING FUNCTION WILL DIRECTLY PARSE THE SEISMIC STATION BASED ON THE SOURCE CHOSEN    

SOURCE_ID=A01 # Chose a DMQ nest ID to run the simulation with. This will automatically initialize the seismometers.
                          #nStation, Stationfile1, SourceID1, Stationfile2, SourceID2, ... 

# Path to txt files with source parameters.
AP12="$BASE/Stations/AP12.txt"
AP14="$BASE/Stations/AP14.txt"
AP15="$BASE/Stations/AP15.txt"
AP16="$BASE/Stations/AP16.txt"

# MAKE SURE THAT THE GATHER RADIUS IS INITIALISED 
getSourceParams 4 $SOURCE_ID $AP12 $AP14 $AP15 $AP16 # Get the source parameters from the files.
source "$BASE/source_params.txt"

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
#MoveTerrawulfLogs              ##  ''

#echo Begin Figure Generation: >> "$LOGFILE"

rwd=`pwd`       # Switch to output directory - the rest of our work
outdir=$(cat $BASE/data/Parallel-runs/latest_run.txt)
outdir="$outdir/process_$OMPI_COMM_WORLD_RANK"
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
sys.path.append("$BASE")
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

# Determine RunID from name of do-script:   (for labeling figs)
SHFILE=RUNID
for file in *.sh; do
    if grep -qx '##  "DO"-script for a Radiative3D run:' "$file"
    then                # Catch ONLY file with R3D do marker
        SHFILE="$file"  # Also catches no more than one do- file
    fi                  #
done
RUNID=${SHFILE%.sh}
RUNID=${RUNID#do-}
echo "Found RunID: $RUNID"

##
echo "Finished processing on core $OMPI_COMM_WORLD_RANK."
echo "Output has been placed in $outdir."
## END
##
