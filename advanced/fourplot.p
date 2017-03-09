set grid
set term png
set output 'im_gauss.png'
unset key
set multiplot layout 1, 2
set xrange [-5:5]
set yrange [0:2.2]
set xlabel "Position"                               
plot 'space_gauss.dat' using 1:2 with lines lw 1.5
set yrange [0:1.8]
set xrange [-5:5]
set xlabel "Momentum"
plot 'freq_gauss.dat' using 1:2 with lines lw 1.5
unset multiplot
set output 'im_triangle.png'
unset key                               
set multiplot layout 1, 2
set yrange [0:1.1]
set xrange [-5:5]  
set xlabel "Position"                             
plot 'space_tri.dat' using 1:2 with lines lw 1.5
set yrange [0:1.1]
set xrange [-5:5]
set xlabel "Momentum"
plot 'freq_tri.dat' using 1:2 with lines lw 1.5
unset multiplot
set output 'im_square.png'
unset key                               
set multiplot layout 1, 2
set yrange [0:2.2]
set xrange [-20:20]                              
set xlabel "Position" 
plot 'space_squ.dat' using 1:2 with lines lw 1.5
set yrange [-0.5:2]
set xrange [-5:5]
set xlabel "Momentum"
plot 'freq_squ.dat' using 1:2 with lines lw 1.5
unset multiplot
set output 'im_double.png'
unset key                               
set multiplot layout 1, 2
set yrange [0:2.2]
set xrange [-10:10]  
set xlabel "Position"                             
plot 'space_doub.dat' using 1:2 with lines lw 1.5
set yrange [-2.5:2.5]
set xrange [-5:5]
set xlabel "Momentum"
plot 'freq_doub.dat' using 1:2 with lines lw 1.5
unset multiplot

