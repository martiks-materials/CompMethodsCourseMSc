set terminal wxt persist
set multiplot layout 1, 2

set xlabel "Iterations (Sample doublings)"
set ylabel "Value of erf(2)"
set grid
plot 'numtrack1.dat' with lp title "Uniform" pt 5, 'numtrack2.dat' with lp title "Importance (Pdf1)" pt 6, 'numtrack3.dat'with lp title "Importance (Pdf2)" pt 7

set grid
set xlabel "Iterations (Sample doublings)"
set ylabel "log(Standard Error)"
plot 'ertrack1.dat' with lp title "Uniform" pt 5, 'ertrack2.dat' with lp title "Importance (Pdf1)" pt 6, 'ertrack3.dat'with lp title "Importance (Pdf2)" pt 7

unset multiplot


set terminal png enhanced
set output "values.png"
set xlabel "Iterations (Sample doublings)"
set ylabel "Value of erf(2)"
set grid
plot 'numtrack1.dat' with lp title "Uniform" pt 5, 'numtrack2.dat' with lp title "Importance (Pdf1)" pt 6, 'numtrack3.dat'with lp title "Importance (Pdf2)" pt 7

set output "errors.png"
set grid
set xlabel "Iterations (Sample doublings)"
set ylabel "log(Standard Error)"
plot 'ertrack1.dat' with lp title "Uniform" pt 5, 'ertrack2.dat' with lp title "Importance (Pdf1)" pt 6, 'ertrack3.dat'with lp title "Importance (Pdf2)" pt 7
