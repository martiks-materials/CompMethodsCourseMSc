set terminal png enhanced
set output "histo1.png"
set border 3
set yzeroaxis
set xrange [0:1]
set grid
set boxwidth 0.001 absolute
set style fill solid 1.0 noborder
bin_width = 0.001;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes

set output "histo2.png"
set border 3
set yzeroaxis
set xrange [0:3.14159]
set grid
set boxwidth 0.00314159 absolute
set style fill solid 1.0 noborder
bin_width = 0.00314159;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'transform.dat' using (rounded($1)):(1) smooth frequency with boxes

set output "histo3.png"
set border 3
set yzeroaxis
set xrange [0:3.14159]
set grid
set boxwidth 0.00314159 absolute
set style fill solid 1.0 noborder
bin_width = 0.00314159;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'reject.dat' using (rounded($1)):(1) smooth frequency with boxes

