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
BASE=`pwd`
echo "Base directory: $BASE"
SINGLE_CORE=0 # Make sure this is set to 0 for parallel runs, and before sourcing.
echo "$SINGLE_CORE" > "$BASE/Params/single_core.txt"
source $BASE/scripts/do-fundamentals.sh
ENV_PATH="$BASE/R3Denv/bin/activate"
PYTHON_EXEC="python3"


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
#SOURCELOC=0,0,-10                       # Source location, defined in the file
                              # Source param choices include 'expl', 'eq',
                              # 'Selby', etc.  Second paramter is usually
                              # isofrac (if range [-1.0,1.0]) or isoangle
                              # (if range [-90,90]).

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=1M                   # Number of phonons to spray (Recommend: 50M)
RECTIME=1200                  # Recording duration of seismometers (seconds).
BINSIZE=0.5                  # Seismometer time-bin size in seconds
GATHER=500.0                  # Terminal gather radius, in kilometers.
MODELPATH="$BASE/Models/MoonModels/SimplifiedModels/SimplfiedISSI_MOON_M2.csv" #Path of the model to be used. Default is SimplifiedPrem.csv in Models/EarthModels/PREM
GRIDSOURCE=GRID_FROMFILE      # Specify that the grid is read from a file. Otherwise use GRID_UNSPEC.



# THE FOLLOWING FUNCTION WILL DIRECTLY PARSE THE SEISMIC STATION BASED ON THE SOURCE CHOSEN                            



FLATTEN=""                    # transformation to depths and velocities.
                              # (Set to "" to disable.)

#Sourcing this will initialise the source and the Apollo stations
source "$BASE/source_params.txt"

#Undefine the Apollo stations because we are not interested.
SEIS1=""
SEIS2=""
SEIS3=""
SEIS4=""


#### SETUP THE SEISMOMETER ARRAY(s):

nAzi=1                      # Number of azimuthal resolution points
                              # (default is 9).  This is the number of
                              # azimuthal arrays to be generated.
nDist=10                    # Number of distance arrays (default is 50).
                              # This is the number of seismometers in each
                              # distance array.
dist1=0                   # Minimum distance (in km) of the first
                              # distance array.                              
dist2=5000                # Maximum distance (in km) of the last. Here I chose approx half circumference of moon.
                              # distance array.  (Note: this is the
                              # distance of the last seismometer in the
                              # last distance array.)
# Function can be found in fundamentals.sh
# This will generate the seismometer array to cover the planet and ehance the number of phonons detected
generateSeisArray $GATHER $nAzi $nDist $dist1 $dist2 # 400 km gather radius, 
                                                     #   9 azimuthal resolution,
                                                     # 50 seismometers in each distance array, 
                                                     #0 to 12000 km distance


RAD=1737.1
echo $RAD > "$BASE/Params/rad.txt"
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


# Convert the seisfiles to pkl format
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


## ___FIG_GEN_END___
##
cd "$rwd"
#echo End Figure Generation: `date +"%Y.%m.%d-%H:%M:%S"`# >> "$LOGFILE"
echo '///'
echo "Finished processing on core $OMPI_COMM_WORLD_RANK."
echo "Figure generation complete. To re-run figures, run 'do-figsonly.sh' in"
echo "the output directory, which may be modified to customize figure output."
echo "Output has been placed in $outdir."
## END
