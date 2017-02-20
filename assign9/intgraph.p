set terminal wxt persist
set multiplot layout 1, 2

set xlabel "x"
set ylabel "y"
set grid
plot 'runge1.dat' with lines title "y1(x)" lw 2, 'runge2.dat' with lines title "y2(x)" lw 2

set grid
set xlabel "y1"
set ylabel "y2"
plot 'runge3.dat' with lines notitle lw 2

unset multiplot


set terminal png enhanced
set output "runge12.png"
set xlabel "x"
set ylabel "y"
set grid
plot 'runge1.dat' with lines title "y1(x)" lw 2, 'runge2.dat' with lines title "y2(x)" lw 2

set output "runge3.png"
set grid
set xlabel "y1"
set ylabel "y2"
plot 'runge3.dat' with lines notitle lw 2


