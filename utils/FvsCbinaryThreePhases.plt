set terminal png
set pointsize 2
set size 0.7,0.7

set output "FvsC.png"
set xlabel "Composition"
set ylabel "Energy (J/mol)"

plot 'FvsC.dat' i 0 u 1:2 t 'phase 1' w lines lw 2, \
     'FvsC.dat' i 1 u 1:2 t 'phase 2' w lines lw 2, \
     'FvsC.dat' i 2 u 1:2 t 'phase 3' w lines lw 2

