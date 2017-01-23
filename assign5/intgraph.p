set terminal wxt persist
set multiplot layout 1, 2

set xlabel "x"
set ylabel "y"
set grid
plot 'runge1.dat' with linespoints title "y1(x)", 'runge2.dat' with linespoints title "y2(x)"

set grid
set xlabel "y1"
set ylabel "y2"
plot 'runge3.dat' with linespoints notitle

unset multiplot


set terminal png enhanced
set output "runge12.png"
set xlabel "x"
set ylabel "y"
set grid
plot 'runge1.dat' with linespoints title "y1(x)", 'runge2.dat' with linespoints title "y2(x)"

set output "runge3.png"
set grid
set xlabel "y1"
set ylabel "y2"
plot 'runge3.dat' with linespoints notitle


