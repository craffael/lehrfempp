
set terminal pdf
set logscale

set xlabel 'N'
set ylabel 'error'

set output 'convergence_l2_regular.pdf'
set title 'L2 convergence'
plot 'output_regular.dat' using 2:3 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_regular.pdf'
set title 'DG-norm convergence'
plot 'output_regular.dat' using 2:4 with lines title 'DG-norm error'


set output 'convergence_l2_regular_modified.pdf'
set title 'L2 convergence'
plot 'output_regular.dat' using 2:5 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_regular_modified.pdf'
set title 'DG-norm convergence'
plot 'output_regular.dat' using 2:6 with lines title 'DG-norm error'


set output 'convergence_l2_irregular.pdf'
set title 'L2 convergence'
plot 'output_irregular.dat' using 2:3 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_irregular.pdf'
set title 'DG-norm convergence'
plot 'output_irregular.dat' using 2:4 with lines title 'DG-norm error'


set output 'convergence_l2_irregular_modified.pdf'
set title 'L2 convergence'
plot 'output_irregular.dat' using 2:5 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_irregular_modified.pdf'
set title 'DG-norm convergence'
plot 'output_irregular.dat' using 2:6 with lines title 'DG-norm error'
