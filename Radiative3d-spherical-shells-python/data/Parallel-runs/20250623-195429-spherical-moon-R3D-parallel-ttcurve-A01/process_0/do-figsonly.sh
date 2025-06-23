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


## ___FIG_GEN_END___
