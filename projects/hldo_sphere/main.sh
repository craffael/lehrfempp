#!/bin/bash

#
# script to run experiment
#

help()
{
    echo -e "Usage: main.sh [-s | --short] [-a | --all] [-h | --help]\n\
\tWithout options the method runs all the experiments described in 
\tthe experiment section of the thesis.\n\n\
\t-s | --short \trun only short version of the experiments using\n\t\t\tless refinement levels resulting in significant speedup\n\
\t-a | --all  \trun all experiments including the ones in the debugging \n\t\t\tsection\n\
\t-p | --plot  \tonly plot sample plots for the hodge laplace experiments\n\t\t\tthis option only works together with -s and requires \n\
\t\t\tpython, click, pandas, numpy and matplotlib to be installed.\n
\t-h | --help \tshow this help"
    exit 2
}

SHORT=s,a,h,p
LONG=short,all,help,plot
OPTS=$(getopt --name main --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

s=0
a=0
p=0

while :
do
    case "$1" in
        -s | --short)
            s=1
            shift 1
            ;;
        -a | --all)
            a=1
            shift 1
            ;;
        -h | --help)
            help
            shift 1
            ;;
        -p | --plot)
            p=1
            shift 1
            ;;
        --)
            shift;
            break
            ;;
        *)
            echo "Unexpected option $1"
            help
            ;;
    esac
done


# The following contains the calls to experiments in the section experiments
if [[ $s -eq 0 ]]; then

    echo -e "\nExecute stability experiments (section 6.2)\n"

    cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
    # arguments: max_refinement_level min_k^2 max_k^2 step_size_in_k
    ./projects.hldo_sphere.experiments.hodge_laplacians.test 4 0.001 16.0 0.1
    cd -

    cd ../../build/projects/hldo_sphere/experiments/dirac_operators
    # arguments: max_refinement_level min_k max_k step_size_in_k
    ./projects.hldo_sphere.experiments.dirac_operators.test 4 0.01 4.0 0.1

    echo -e "\nExecute convergence experiments (section 6.3)\n"

    cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
    # arguments: max_refinement_level min_k^2 max_k^2
    ./projects.hldo_sphere.experiments.hodge_laplacians.test 6 0.25 0.25
    cd -

    cd ../../build/projects/hldo_sphere/experiments/dirac_operators
    # arguments: max_refinement_level min_k max_k
    ./projects.hldo_sphere.experiments.dirac_operators.test 6 0.5 0.5
    cd -


    echo -e "\nExecute Hodge and Dirac experiments (section 6.4)\n"

    cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
    # arguments: max_refinement_level min_k^2 max_k^2
    ./projects.hldo_sphere.experiments.hodge_laplacians.testdirac 6 0.25 0.25
    cd -

    cd ../../build/projects/hldo_sphere/experiments/hodge_and_dirac
    # arguments: max_refinement_level min_k max_k
    ./projects.hldo_sphere.experiments.hodge_and_dirac.test 5 0.5 0.5
    cd -

    if [[ $a -eq 1 ]]; then

        echo -e "\nExecute debugging experiment convergence of the solution\
(section 5.3)\n"

        cd ../../build/projects/hldo_sphere/debugging/experiments
        # arguments: max_refinement_level k
        ./projects.hldo_sphere.debugging.experiments.convergence_test 6 0.5
        cd -

        echo -e "\nExecute debugging experiment bilinear form convergence\
(section 5.4)\n"

        cd ../../build/projects/hldo_sphere/debugging/experiments
        # arguments: max_refinement_level
        ./projects.hldo_sphere.debugging.experiments.basis_expansion_test 6
        cd -

        cd ../../build/projects/hldo_sphere/debugging/experiments
        # arguments: max_refinement_level
        ./projects.hldo_sphere.debugging.experiments.whitney_one_curl_test 6
        cd -

    fi

else

    echo -e "\nExecute SHORTENED stability experiments (section 6.2)\n"

    cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
    # arguments: max_refinement_level min_k^2 max_k^2 step_size_in_k
    ./projects.hldo_sphere.experiments.hodge_laplacians.test 3 0.001 16.0 0.5
    cd -

    cd ../../build/projects/hldo_sphere/experiments/dirac_operators
    # arguments: max_refinement_level min_k max_k step_size_in_k
    ./projects.hldo_sphere.experiments.dirac_operators.test 3 0.01 4.0 0.5
    cd -

    echo -e "\nExecute SHORTENED convergence experiments (section 6.3)\n"

    cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
    # arguments: max_refinement_level min_k^2 max_k^2
    ./projects.hldo_sphere.experiments.hodge_laplacians.test 4 0.25 0.25
    cd -

    cd ../../build/projects/hldo_sphere/experiments/dirac_operators
    # arguments: max_refinement_level min_k max_k
    ./projects.hldo_sphere.experiments.dirac_operators.test 4 0.5 0.5
    cd -

    echo -e "\nExecute SHORTENED Hodge and Dirac experiments (section 6.4)\n"

    cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
    # arguments: max_refinement_level min_k^2 max_k^2
    ./projects.hldo_sphere.experiments.hodge_laplacians.testdirac 4 0.25 0.25
    cd -

    cd ../../build/projects/hldo_sphere/experiments/hodge_and_dirac
    # arguments: max_refinement_level min_k max_k
    ./projects.hldo_sphere.experiments.hodge_and_dirac.test 4 0.5 0.5
    cd -

    if [[ $a -eq 1 ]]; then

        echo -e "\nExecute SHORTENED debugging experiment convergence of\
 the solution (section 5.3)\n"

        cd ../../build/projects/hldo_sphere/debugging/experiments
        # arguments: max_refinement_level k
        ./projects.hldo_sphere.debugging.experiments.convergence_test 3 0.5
        cd -

        echo -e "\nExecute SHORTENED debugging experiment bilinear\
 form convergence (section 5.4)\n"

        cd ../../build/projects/hldo_sphere/debugging/experiments
        # arguments: max_refinement_level
        ./projects.hldo_sphere.debugging.experiments.basis_expansion_test 3
        cd -

        cd ../../build/projects/hldo_sphere/debugging/experiments
        # arguments: max_refinement_level
        ./projects.hldo_sphere.debugging.experiments.whitney_one_curl_test 3
        cd -

    fi

    # If the plot parameter is chosen just plot selected plots of the thesis just example calls
    if [[ $p -eq 1 ]]; then

        cd ../../build/projects/hldo_sphere/experiments/hodge_laplacians
        python ./ploteigenvalues.py --file ./results/result_test_0_0316228-3_93713.csv --stat sd --noshow
        echo "Created image ../../build/projects/hldo_sphere/experiments/hodge_laplacians/eigenvalslap.png"
        python ./plotlaplace.py --file ./results/result_test_0_5-0_5.csv --log --noshow
        echo "Created image ../../build/projects/hldo_sphere/experiments/hodge_laplacians/hodgelaplacians.png"
        cd -

    fi
fi
