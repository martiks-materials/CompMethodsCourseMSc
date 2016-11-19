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
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Uniform deviate"

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
plot 'transform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Transformed Dist."

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
plot 'reject.dat' using (rounded($1)):(1) smooth frequency with boxes title "Rejection Method Dist"

unset multiplot

set terminal png enhanced
set output "histo1.png"
set border 3
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
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Uniform deviate"

set output "histo2.png"
set border 3
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
plot 'transform.dat' using (rounded($1)):(1) smooth frequency with boxes title "Transformed Dist."

set output "histo3.png"
set border 3
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
plot 'reject.dat' using (rounded($1)):(1) smooth frequency with boxes title "Rejection Method Dist"

