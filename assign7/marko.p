set terminal png enhanced
set output "marko.png"
set xlabel "x_1"
set ylabel "x_2"
set grid
set title "MCMC of Rosenbrock for 12 Markov chains"
plot "markovdata.dat" using 1:2 with p title "Function Sampling"

set output "burnin.png"
set xlabel "x_1"
set ylabel "x_2"
set title "Burn-in data for MCMC of the Rosenbrock for 12 Markov chains"
plot "burndata.dat" using 1:2 with p title "Burn-in Period"

set output "multivar.png"
set xlabel "x_1"
set ylabel "x_2"
set title "MCMC of Rosenbrock (12 chains) with Multivariate Gaussian"
plot "multivar.dat" using 1:2 with p title "Function Sampling"

set output "multivarburn.png"
set xlabel "x_1"
set ylabel "x_2"
set title "Burn-in data for Rosenbrock MCMC (12 chains), Multivariate Gaussian"
plot "multivarburn.dat" using 1:2 with p title "Burn-in Period"
