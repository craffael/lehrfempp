set terminal pdf
set logscale

set xlabel 'N'
set ylabel 'error'

set output 'convergence_l2_zero.pdf'
set title 'L2 convergence'
plot 'output_builder.dat' using 1:2 with lines title 'L2 error'

set output 'convergence_dg_zero.pdf'
set title 'DG-norm convergence'
plot 'output_builder.dat' using 1:3 with lines title 'DG-norm error'

set output 'convergence_l2_nonzero.pdf'
set title 'L2 convergence'
plot 'output_builder.dat' using 1:4 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_nonzero.pdf'
set title 'DG-norm convergence'
plot 'output_builder.dat' using 1:5 with lines title 'DG-norm error'


set output 'convergence_l2_zero_modified.pdf'
set title 'L2 convergence'
plot 'output_builder.dat' using 1:6 with lines title 'L2 error'

set output 'convergence_dg_zero_modified.pdf'
set title 'DG-norm convergence'
plot 'output_builder.dat' using 1:7 with lines title 'DG-norm error'

set output 'convergence_l2_nonzero_modified.pdf'
set title 'L2 convergence'
plot 'output_builder.dat' using 1:8 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_nonzero_modified.pdf'
set title 'DG-norm convergence'
plot 'output_builder.dat' using 1:9 with lines title 'DG-norm error'



set output 'convergence_l2_zero_files.pdf'
set title 'L2 convergence'
plot 'output_files.dat' using 1:2 with lines title 'L2 error'

set output 'convergence_dg_zero_files.pdf'
set title 'DG-norm convergence'
plot 'output_files.dat' using 1:3 with lines title 'DG-norm error'

set output 'convergence_l2_nonzero_files.pdf'
set title 'L2 convergence'
plot 'output_files.dat' using 1:4 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_nonzero_files.pdf'
set title 'DG-norm convergence'
plot 'output_files.dat' using 1:5 with lines title 'DG-norm error'


set output 'convergence_l2_zero_modified_files.pdf'
set title 'L2 convergence'
plot 'output_files.dat' using 1:6 with lines title 'L2 error'

set output 'convergence_dg_zero_modified_files.pdf'
set title 'DG-norm convergence'
plot 'output_files.dat' using 1:7 with lines title 'DG-norm error'

set output 'convergence_l2_nonzero_modified_files.pdf'
set title 'L2 convergence'
plot 'output_files.dat' using 1:8 with lines title 'L2 error', \
     x**(-0.5) title 'N^{-1/2}' dt 2

set output 'convergence_dg_nonzero_modified_files.pdf'
set title 'DG-norm convergence'
plot 'output_files.dat' using 1:9 with lines title 'DG-norm error'
