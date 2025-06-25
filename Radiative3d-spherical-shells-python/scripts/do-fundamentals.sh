#!/bin/bash
##
##  DO-SCRIPT FUNDAMENTALS
##
##  (Library of routines for Radiative3D Do-Scripts.)
##
##  Include (source) this file in your own do-scripts, with a line like this:
##
##    source scripts/do-fundamentals.sh
##
##
R3DBASE=`pwd`  # Base directory for the run.
echo "R3D BASE DIRECTORY: $R3DBASE"
export R3DBASE

SINGLE_CORE=$(cat $R3DBASE/Params/single_core.txt)
echo "SINGLE_CORE: $SINGLE_CORE"

if [[ "$SINGLE_CORE" -eq 1 ]]; then
    echo "Running in single-core mode."
    OUTPARALLEL=$(cat "$R3DBASE/data/Single-core-runs/latest_run.txt")
    mkdir -p "$OUTPARALLEL"
    echo "Storing results in $OUTPARALLEL"
else
    PROCNB=$OMPI_COMM_WORLD_RANK
    if [ -z "$PROCNB" ]; then
        echo "OMPI_COMM_WORLD_RANK is not set. Running in single-core mode. Exiting."
        exit 1  
    fi
    
    echo "Running in parallel mode."
    echo " Hello from process $PROCNB !"
    OUTPARALLEL=$(cat $R3DBASE/data/Parallel-runs/latest_run.txt)
    
    OUTPARALLEL="$OUTPARALLEL/process_$PROCNB"
    mkdir -p "$OUTPARALLEL/process_$PROCNB"
    echo "Storing results in $OUTPARALLEL/process_$PROCNB"
fi

######
## Helpers:
##
function SecondsToHours {
    printf %.2f\\n "$((10**9 * $1/3600))e-9"
}
function SecondsDeltaToHours {
    local Delta=$(($1 - $2))
    echo $(SecondsToHours $Delta)
}
function VCSVersion {
    if [ -d .svn ]; then
        svnversion
    elif [ -d .git ]; then
        git rev-parse --short HEAD
    else
        echo "Unversioned directory"
    fi
}
function VCSStatusInfo {
    if [ -d .svn ]; then
        echo "Revision number (svnversion): $(svnversion)"
        svn info
        svn st
        svn diff
    elif [ -d .git ]; then
        echo "Revision number (git): $(git rev-parse --short HEAD)"
        git status
        git diff
    else
        echo "Unversioned directory"
    fi
}
function RealPath {
    # realpath not a base install on Mac (although can be added with
    # `brew install coreutils`). We fall back to just echoing the path if
    # not installed.
    if which realpath > /dev/null; then
        realpath "$1"
    else
        echo "$1"
    fi
}

######
## FUNCTION:  PopDefaults()
##
##   Populates default values for control veriables that are not
##   already specified, based on a selected "target" simulation type.
##   E.g., certain command line option and their values are determined
##   mainly by whether you are producing envelopes or videos, but do
##   not otherwise change from one simulation to another. These are
##   packed into here and chosen via the $1 argument.
##
##   SIDE AFFECTS:
##
##   Sets (but will not override if already exists) $R3D_EXE,
##   $OUT_BASE, $TOA_DEGREE, $REPORTS, (and others, tbd)
##
##
function PopDefaults {

    [ -z "$R3D_EXE" ]  && R3D_EXE=$R3DBASE/main  # Executable to use
    [ -z "$OUT_BASE" ] && OUT_BASE=$OUTPARALLEL  # Base directory for output
                                                                         # to put output
                                                                         # subdirectories

    case "$1" in
        waveform)
            [ -z "$TOA_DEGREE" ] && TOA_DEGREE=8
            [ -z "$REPORTS" ]    && REPORTS=INV
            [ -z "$FINTAG" ]     && FINTAG=R3D
            ;;
        video)
            [ -z "$TOA_DEGREE" ] && TOA_DEGREE=8
            [ -z "$REPORTS" ]    && REPORTS=ALL_ON
            [ -z "$FINTAG" ]     && FINTAG=R3VID
            ;;
        template)
            [ -z "$FINTAG" ]     && FINTAG=R3BATCH
            # No further defualts if we're just templating a batch.
            ;;
        *)
            echo Simulation target not recognized.
            exit
            ;;
    esac

}


######
## FUNCTION:  EventParams()
##
##   A selector for named event types.  First arg is event tag.
##   Subsequent args depend on selection.
##
##   USAGE:
##
##     SOURCETYP=$(EventParams [args])
##
##
function EventParams {

    local result="EQ"
    local isofrac=0

    case "$1" in
        expl)   # Generic explosion
            SOURCETYP=EXPL                    #
            ;;
        Selby|selby)
            if [ $# -gt 1 ]; then
                isofrac=$2
            fi
            result=SDR,125,40,90,$isofrac
            ;;
        *)
            echo Event code not recognized.
            exit
            ;;
    esac
    echo "$result"
}


######
## FUNCTION:  CheckBuildStatus()
##
##   Checks whether Make indicates that $R3D_EXE is up-to-date.  Exits
##   script with a message if it is not.
##
function CheckBuildStatus {
    if ! make -q
    then
        echo Make indicates $R3D_EXE is not up to date.
        echo Please run make and re-run this script.
        exit
    fi
}


######
## FUNCTION:  CreateOutputDirectory()
##
##   Sets up directory in which to put simulation output, graphical
##   results, and and in which to store input records and scripts,
##   etc.  Directory will be named with a timestamp and an extra token
##   if $1 is set.  Useage:
##
##     CreateOutputDirectory $@
##
##   (The $@ is needed to send script args to this function, the main
##   one is a tag id in the output directory name.)
##
##   ASSUMPTIONS:
##
##   Function assumes that $OUT_BASE and $R3D_EXE are already defined.
##
##   SIDE EFFECTS:
##
##   A directory will be created and populated. $R3D_EXE will be
##   modified.  $outdir will be set.  A logfile will be started and
##   $LOGFILE will be set.
##
function CreateOutputDirectory {

    #local SHFILE="$(basename "$0")"     # Guess a default RunID from
    #local RUNID=${SHFILE%.sh}           # the shell file name,
    #      RUNID=${RUNID#do-}            #
    #[ -z "$1" ] || RUNID="$1"           # but use $1 if it was provided.

    #[ -z "$FINTAG" ] && FINTAG=R3D
    #local tstamp=`date +"%Y%m%d-%H%M%S"`
    #local outdirname="$tstamp"-"$RUNID"-"$FINTAG"
    outdir="$OUT_BASE"
    #mkdir -p "$outdir"
    #local rcptfile="`basename ${0%.sh}`.$$.runrcpt"
    #echo `hostname`:`RealPath "$outdir"` >> $rcptfile

    #if [ -L "$OUT_BASE/previous" ]; then    # Generate "latest" and "previous"
    #    rm "$OUT_BASE/previous"             # links to make it easier to cd
    #fi                                      # into new output
    #if [ -L "$OUT_BASE/latest" ]; then
    #    mv "$OUT_BASE/latest" "$OUT_BASE/previous"
    #fi
    #ln -sf "$outdirname" "$OUT_BASE/latest"
    #

    cp "$0" "$outdir"/                # Copy do-script
    cp $R3D_EXE "$outdir"/            # Keep the exe for replicability
    R3D_EXE="$outdir"/$(basename $R3D_EXE)  # Run the copy instead of the orig
    echo "(BRIEF description of run goes here.)" > "$outdir"/description.txt
    [ -z "$INTENT" ] || echo "$INTENT" >> "$outdir"/description.txt
    [ -z "$CAMPAIGN" ] || echo "$CAMPAIGN" >> "$outdir"/description.txt
    VCSStatusInfo > "$outdir"/svn-info.txt # Record detailed revision info

    LOGFILE="$outdir"/logfile
    > "$LOGFILE"
    echo "Machine: " `hostname` >> "$LOGFILE"
    echo "User: " `whoami` >> "$LOGFILE"
    echo "Revision: " `VCSVersion` >> "$LOGFILE"

    # Carve-off and save figure-generation part of do-script in case
    # we want to re-run the figures later (as might happen if we
    # improve or tweak the plotting routines after a run).
    cat "$0" | grep -A 10000 -B 1 -E "^## ___FIG_GEN_START___$" | \
               grep -A 0 -B 10000 -E "^## ___FIG_GEN_END___$" \
               > "$outdir"/do-figsonly.sh
    chmod u+x "$outdir"/do-figsonly.sh

}


######
## FUNCTION:  RunSimulation()
##
##   Runs the actual Radiative3D simulation. Logs the begin time, run
##   time, and command line to the log file.
##
##   ASSUMPTIONS:
##
##   Assumes that $R3D_EXE and $LOGFILE have been defined and properly
##   set, and that all variables used in defining the command line
##   have been defined.  Typically, these are defined in the body of
##   the main do-script.
##
##   SIDE EFFECTS:
##
##   Begins simulation, appends to >>$LOGFILE
##
##
function RunSimulation {
    echo Begin Radiative3D Run: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"
    local binsizearg=""
    [ -z "$BINSIZE" ] || binsizearg="--binsize=$BINSIZE"
    local cylradarg=""
    [ -z "$CYLRAD" ] || cylradarg="--range=$CYLRAD"
    local modelarg=""
    [ -z "$MODIDX" ] || modelarg="--grid-compiled=$MODIDX" 
    local TIMEBEGIN=`date +"%s"`
    local R3D_CMDLN="\
$R3D_EXE --reports=$REPORTS \
        --output-dir=\"$outdir\" \
        --report-file=reports.dat \
        --output-envelope=1 \
        --mparams-outfile=out_mparams.octv \
        --num-phonons=$NUMPHONS \
        --toa-degree=$TOA_DEGREE \
        --source=$SOURCETYP \
        --source-loc=$SOURCELOC \
        --frequency=$FREQ \
        --timetolive=$RECTIME \
        --model-path=$MODELPATH\
        $binsizearg \
        $modelarg \
        $cylradarg \
        $FLATTEN \
        $MFPOVERRIDE \
        $NODEFLECT \
        --model-args=$COMPARGS \
        --dump-grid \
        $ADDITIONAL \
"
    echo $R3D_CMDLN >> "$LOGFILE"
    eval "$R3D_CMDLN" | tee "$outdir"/stdout.txt
    echo End__ Radiative3D Run: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"
    local TIMEEND=`date +"%s"`
    local TIMEELAPSED=$(SecondsDeltaToHours $TIMEEND $TIMEBEGIN)
    echo "($TIMEELAPSED Hours Run-Time)" >> "$LOGFILE"

}

######
## FUNCTION:  generateSeisArray()
##
##   Generates a set of seismometer parameters for the simulation.
##   The parameters are generated based on the input arguments.
##   The function creates a set of variables SEIS1, SEIS2, ..., SEISn
##   where n is the number of seismometers. The parameters are
##   generated based on the input arguments. The function also
##   generates an ADDITIONAL string that contains the parameters for
##   the seismometers.
##
function generateSeisArray() { #1 : gather radius, #2 : Azimuthal resolution, #3 : distance resolution, #4 : dist1, #5 : dist2
    local gather=$1

    local range1=$4
    local range2=$5
    local max_index=$2
    local step=$(echo "180 / $max_index" | bc)
    local base_params="20.0,20.0,$gather,$3"
#    echo "$base_params" >> $BASE/array_params.txt
    
    for i in $(seq 1 $max_index); do #seq is used to generate a sequence of numbers
        local azimuth=$(echo "($i - 1) * $step" | bc) #bc is used for floating point arithmetic. Converts the echo output to a number
        local seisorig="$range1,$azimuth,0"
        local seisdest="$range2,$azimuth,0"
        local seis_var="--seis-p2p=$seisorig,$seisdest,$base_params"
        eval "SEIS$i='$seis_var'"
    done

    # Build ADDITIONAL using the defined SEIS variables
    ADDITIONAL=""
    for i in $(seq 1 $max_index); do
        ADDITIONAL+="\${SEIS$i} "

    done
}

######
## FUNCTION:  get_multiple_source_params()
##
##   Retrieves parameters for multiple sources from specified files.
##   The function reads the source parameters from the files and generates the corresponding SEIS variables. 
##   It also sets the SOURCELOC variable for the first source.
##  USAGE:
## Need to construct for each station file, e.g., Appolo 12, need to have 4 columns : ID, Depth, Distance, Azimuth

## If you have 4 stations, you need 4 files.

getSourceParams() { #$1 : number of stations, #2 SourceID #$3 : file1 $4 : file2, $5 : file3, $6 : file4 ...
  # This function retrieves parameters for multiple sources from specified files.
  local n=$1
  local source=$2
  shift 2
  local PARAMFILE="$R3DBASE/source_params.txt"
 ADDITIONAL=""
  for ((i=1; i<=n; i++)); do
    local file=$1
    shift 

    if [[ ! -f "$file" ]]; then
      echo "File '$file' not found."
      return 1
    fi

    local line
    line=$(grep "^$source " "$file")

    if [[ -z "$line" ]]; then
      echo "Source '$source' not found in '$file'."
      return 1
    fi

    local depth distance azimuth
    read -r _ depth distance azimuth <<< "$line"

    eval "export SEISORIG${i}=\"$distance,$azimuth,0\""
    eval "export SEISDEST${i}=\"$((distance+1)),$azimuth,0\""
    eval "export SEIS${i}=\"--seis-p2p=\$SEISORIG${i},\$SEISDEST${i},1.0,2.0,\$GATHER,2\""
    varname="SEIS${i}"
    echo "$varname=${!varname}" >> "$PARAMFILE"
    ADDITIONAL+="\${SEIS$i} "
    if [[ $i -eq 1 ]]; then
      echo "SOURCELOC=0,0,-$depth" >> "$PARAMFILE"
    fi
done
}
######
## FUNCTION:  SpiralIndex n=$1 i0=$2 j0=$3
##
##  Given an index n as a counter in an i,j index space, return "i j" after
##  tracing a spiral path from centerpoint index pair (i0,j0). This is used
##  to iterate a 2-D parameter space "from the center outward" so that we
##  can view early results clustered around the center and the later results
##  will paint the perimeter.  n=0 returns "i0 j0". n>0 traces path like:
##
##      13 - 14 - 15 ...         +i -->
##       |                    +j
##      12    2 -- 3 -- 4      |
##       |    |         |     \ /
##      11    1 -- 0    5      V
##       |              |
##      10 -- 8 -- 7 -- 6
##
function SpiralIndex {

    local n="$1"
    [ "$n" -ge "0" ] || return  # Integer check

    local i="$2"
    local j="$3"
    [ "$i" -ge "0" ] || return  # Integer check
    [ "$j" -ge "0" ] || return  # Integer check

    if [ "$n" -eq "0" ]; then
        echo "$i $j"
        return
    fi

    local p=1
    local step=-1
    local idx=0
    while [ "$p" -lt "40" ]; do  # (Set large sanity limit)
        local q=1
        local qmax="$p"; ((qmax*=2))
        while [ "$q" -le "$qmax" ]; do
            if [ "$q" -le "$p" ]; then
                ((j+=step))
            else
                ((i+=step))
            fi
            ((idx+=1))
            if [ "$idx" -eq "$n" ]; then
                echo $i $j
                return
            fi
            ((q+=1))
        done
        ((p+=1))
        ((step*=-1))
    done
    # shouldn't get here unless n too big

}
## END
##