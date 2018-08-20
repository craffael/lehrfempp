#!/bin/bash

# Run this script in the LehrFEM++ root folder

TESTPRG=`find . -name "lf.refinement.rrt_test" -print | head -1`;
if [[ -z $TESTPRG ]]; then
    echo "This script must be run in the LehrFEM ++ root directory";
else
    TESTDIR="`dirname $TESTPRG`";
    PLOTLFMESH=`find . -name "plot_lf_mesh.m" -print | head -1`;
    PLOTLFMESHDIR="`dirname $PLOTLFMESH`";
    PLOTMULTIREF=`find . -name "plot_multiref.m" -print | head -1`;
    PLOTMULTIREFDIR="`dirname $PLOTMULTIREF`";
    cd $TESTDIR || { echo "$TESTDIR inaccessible"; exit -1; }
    # Run the LehrFEM++ test code
    echo "<< Running testing code in `pwd`>>";
    # ./lf.refinement.rrt_test && { echo "<<< Test run successful >>>"; }
    # Run the matlab postprocessing
    echo "<< Running matlab for postprocessing >>";
    # echo "$PLOTLFMESH, $PLOTMULTIREF";
    # matlab -nosplash -nodesktop -r "addpath($PLOTLFMESHDIR,$PLOTMULTIREFDIR); $PLOTMULTIREF"; 
fi
