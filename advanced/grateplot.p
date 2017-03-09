set grid
do for [i=1:9] {
  set term png
  set output 'gratimag/im_000'.i.'.png'
  unset key
  set multiplot layout 1, 2
  set xrange [-20:20]
  set yrange [0:2.2]                               
  plot 'gratspace/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:1.8]
  set xrange [-5:5]
  plot 'gratmom/mom_'.i.'.data' using 1:2 with lines lw 1.5
  unset multiplot
}
do for [i=10:99] {
  set term png
  set output 'gratimag/im_00'.i.'.png'
  unset key                               
  set multiplot layout 1, 2
  set yrange [0:2.2]
  set xrange [-20:20]                               
  plot 'gratspace/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:1.8]
  set xrange [-5:5]
  plot 'gratmom/mom_'.i.'.data' using 1:2 with lines lw 1.5
  unset multiplot
}
do for [i=100:599] { 
  set term png
  set output 'gratimag/im_0'.i.'.png'
  unset key                               
  set multiplot layout 1, 2
  set yrange [0:2.2]
  set xrange [-20:20]                               
  plot 'gratspace/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:1.8]
  set xrange [-5:5]
  plot 'gratmom/mom_'.i.'.data' using 1:2 with lines lw 1.5
  unset multiplot
}
