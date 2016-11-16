set terminal wxt persist
set multiplot layout 3, 1
set yzeroaxis
set xrange [0:1]
set grid
set boxwidth 0.001 absolute
set style fill solid 1.0 noborder
bin_width = 0.001;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes

plot 'transform.dat' using (rounded($1)):(1) smooth frequency with boxes

plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes


