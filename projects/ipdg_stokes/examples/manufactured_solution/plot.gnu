
set terminal pdf
set logscale

set xlabel 'N'
set ylabel 'error'

set output 'convergence_l2.pdf'
set title 'L2 convergence'
plot 'output.dat' using 2:3 with lines title 'L2 error'

set output 'convergence_dg.pdf'
set title 'DG-norm convergence'
plot 'output.dat' using 2:4 with lines title 'DG-norm error'

set output 'convergence_l2f.pdf'
set title 'L2 convergence'
plot 'output.dat' using 2:5 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dgf.pdf'
set title 'DG-norm convergence'
plot 'output.dat' using 2:6 with lines title 'DG-norm error'


set output 'convergence_l2_modified.pdf'
set title 'L2 convergence'
plot 'output.dat' using 2:7 with lines title 'L2 error'

set output 'convergence_dg_modified.pdf'
set title 'DG-norm convergence'
plot 'output.dat' using 2:8 with lines title 'DG-norm error'

set output 'convergence_l2f_modified.pdf'
set title 'L2 convergence'
plot 'output.dat' using 2:9 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dgf_modified.pdf'
set title 'DG-norm convergence'
plot 'output.dat' using 2:10 with lines title 'DG-norm error'
