do for [i=1:9] {
  set term png
  set output 'seepimag/im_0000'.i.'.png'
  unset key
  set multiplot layout 3, 1 
  set xrange [-256:256]
  set yrange [0:0.5]                               
  plot 'seepdat/space/seepspace1/space_'.i.'.data' using 1:2 with lines lw 1.5, 'seepdat/space/seepspace2/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:60]
  set xrange [0:0.5]
  set style fill solid
  plot 'seepdat/mom/seepmom1/mom_'.i.'.data' using 1:2 with boxes
  set yrange [0:60]
  set xrange [0:0.5]
  plot 'seepdat/mom/seepmom2/mom_'.i.'.data' using 1:2 with boxes
  unset multiplot
}

do for [i=10:99] {
  set term png
  set output 'seepimag/im_000'.i.'.png'
  unset key
  set multiplot layout 3, 1
  set xrange [-256:256]
  set yrange [0:0.5]
  unset style fill
  plot 'seepdat/space/seepspace1/space_'.i.'.data' using 1:2 with lines lw 1.5, 'seepdat/space/seepspace2/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:60]
  set xrange [0:0.5]
  set style fill solid
  plot 'seepdat/mom/seepmom1/mom_'.i.'.data' using 1:2 with boxes
  set yrange [0:60]
  set xrange [0:0.5]
  plot 'seepdat/mom/seepmom2/mom_'.i.'.data' using 1:2 with boxes
  unset multiplot
}

do for [i=100:999] { 
  set term png
  set output 'seepimag/im_00'.i.'.png'
  unset key                               
  set multiplot layout 3, 1 
  set xrange [-256:256]
  set yrange [0:0.5]                               
  unset style fill
  plot 'seepdat/space/seepspace1/space_'.i.'.data' using 1:2 with lines lw 1.5, 'seepdat/space/seepspace2/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:60]
  set xrange [0:0.5]
  set style fill solid
  plot 'seepdat/mom/seepmom1/mom_'.i.'.data' using 1:2 with boxes
  set yrange [0:60]
  set xrange [0:0.5]
  plot 'seepdat/mom/seepmom2/mom_'.i.'.data' using 1:2  with boxes
  unset multiplot
}
do for [i=1000:9999] { 
  set term png
  set output 'seepimag/im_0'.i.'.png'
  unset key                               
  set multiplot layout 3, 1
  set xrange [-256:256]
  set yrange [0:0.5]                               
  unset style fill
  plot 'seepdat/space/seepspace1/space_'.i.'.data' using 1:2 with lines lw 1.5, 'seepdat/space/seepspace2/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:60]
  set xrange [0:0.5]
  set style fill solid
  plot 'seepdat/mom/seepmom1/mom_'.i.'.data' using 1:2  with boxes
  set yrange [0:60]
  set xrange [0:0.5]
  plot 'seepdat/mom/seepmom2/mom_'.i.'.data' using 1:2 with boxes
  unset multiplot
}
do for [i=10000:19523] { 
  set term png
  set output 'seepimag/im_'.i.'.png'
  unset key                               
  set multiplot layout 3, 1 
  set xrange [-256:256]
  set yrange [0:0.5]                               
  unset style fill
  plot 'seepdat/space/seepspace1/space_'.i.'.data' using 1:2 with lines lw 1.5, 'seepdat/space/seepspace2/space_'.i.'.data' using 1:2 with lines lw 1.5
  set yrange [0:60]
  set xrange [0:0.5]
  set style fill solid
  plot 'seepdat/mom/seepmom1/mom_'.i.'.data' using 1:2 with boxes
  set yrange [0:60]
  set xrange [0:0.5]
  plot 'seepdat/mom/seepmom2/mom_'.i.'.data' using 1:2 with boxes
  unset multiplot
}
