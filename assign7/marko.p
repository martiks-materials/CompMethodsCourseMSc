set terminal png enhanced
set output "marko.png"
set xlabel "x_1"
set ylabel "x_2"
set grid
plot "markovdata.dat" using 1:2 with p title "Function Sampling"

set output "burnin.png"
set xlabel "x_1"
set ylabel "x_2"
plot "burndata.dat" using 1:2 with p title "Burn-in Period"

set output "multivar.png"
set xlabel "x_1"
set ylabel "x_2"
plot "multivar.dat" using 1:2 with p title "Burn-in Period"

set output "multivarburn.png"
set xlabel "x_1"
set ylabel "x_2"
plot "multivarburn.dat" using 1:2 with p title "Burn-in Period"
