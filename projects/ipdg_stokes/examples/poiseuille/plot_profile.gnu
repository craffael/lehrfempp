
set terminal svg
set datafile separator ','

set xlabel 'y'
set ylabel 'u'

set output 'velocity_profile.svg'
set title 'Velocity Profile'
plot for [i=1:10] \
    'magnitude.csv' using 'Points:1':(column(3*i+1)) with lines title 'refinement level = '.i
