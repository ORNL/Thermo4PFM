set terminal png
set pointsize 2
set size 0.7,0.7

set output "FvsC.png"
set xlabel "Composition"
set ylabel "Energy (J/mol)"
set datafile separator ","
set key autotitle columnhead

plot 'FvsC.csv' u 1:2 t 'phase L' w lines lw 2, \
     'FvsC.csv' u 1:3 t 'phase A' w lines lw 2, \
     'FvsC.csv' u 1:4 t 'phase B' w lines lw 2

