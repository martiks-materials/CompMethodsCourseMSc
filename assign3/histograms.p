set terminal wxt persist
set multiplot layout 3, 1

set yzeroaxis
set xrange [0:1]
set xlabel "Uniform Deviate x bins"
set ylabel "Frequency of each bin"
set grid
set boxwidth 0.01 absolute
set style fill solid 1.0 noborder
bin_width = 0.01;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
f(x) = 1000
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Uniform deviate", f(x) lw 2.5 lt -1  title "10^5*dx"

set yzeroaxis
set xrange [0:3.14159]
set grid
set xlabel "Transformed variable Y bins"
set ylabel "Frequency of each bin"
set boxwidth 0.0314159 absolute
set style fill solid 1.0 noborder
bin_width = 0.0314159;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
g(x) = 50000*sin(x)*0.0314159
plot 'transform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Transformed Dist.", g(x) lw 2.5 lt -1 title "0.5sin(x)dx"


set yzeroaxis
set xrange [0:3.14159]
set grid
set xlabel "Random Variable Y Obtained Via Rejection Method"
set ylabel "Frequency of each bin"
set boxwidth 0.0314159 absolute
set style fill solid 1.0 noborder
bin_width = 0.0314159;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
h(x)= (200000/3.141592)*sin(x)*sin(x)*0.0314159
plot 'reject.dat' using (rounded($1)):(1) smooth frequency with boxes title "Rejection Method Dist", h(x) lw 2.5 lt -1 title "(2/pi)sin^2(x)dx"


unset multiplot

set terminal png enhanced
set output "histo1.png"
set border 3
set title "Histogram of 10^5 Uniform Deviates (bin width dx=0.01"
set yzeroaxis
set xrange [0:1]
set xlabel "Uniform Deviate x bins"
set ylabel "Frequency of each bin"
set grid
set boxwidth 0.01 absolute
set style fill solid 1.0 noborder
bin_width = 0.01;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Uniform deviate", f(x) lw 2.5 lt -1 title "10^5*dx"

set output "histo2.png"
set border 3
set yzeroaxis
set xrange [0:3.14159]
set grid
set title "Histogram of 10^5 random numbers via transformation method (dx=pi/100)"
set xlabel "Transformed variable Y bins"
set ylabel "Frequency of each bin"
set boxwidth 0.0314159 absolute
set style fill solid 1.0 noborder
bin_width = 0.0314159;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'transform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Transformed Dist.", g(x) lw 2.5 lt -1 title "0.5sin(x)dx"


set output "histo3.png"
set border 3
set yzeroaxis
set xrange [0:3.14159]
set grid
set title "Histogram of 10^5 random numbers via rejection method (dx=pi/100)"
set xlabel "Random Variable Y Obtained Via Rejection Method"
set ylabel "Frequency of each bin"
set boxwidth 0.0314159 absolute
set style fill solid 1.0 noborder
bin_width = 0.0314159;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'reject.dat' using (rounded($1)):(1) smooth frequency with boxes title "Rejection Method Dist", h(x) lw 2.5 lt -1 title "(2/pi)sin^2(x)dx"

