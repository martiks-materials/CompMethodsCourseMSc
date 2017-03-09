
do for [i=0:9] {
  set term png
  set output 'images/im_000'.i.'.png'
  unset key
  set yrange [0:0.1]                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines lw 1.8, 'potential.data' with lines lw 1
}
do for [i=10:99] {
  set term png
  set output 'images/im_00'.i.'.png'
  unset key                               
  set yrange [0:0.10]                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines lw 1.8, 'potential.data' with lines lw 1
}
do for [i=100:999] {
  set term png
  set output 'images/im_0'.i.'.png'
  unset key                               
  set yrange [0:0.10]                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines lw 1.8, 'potential.data' with lines lw 1
}
do for [i=1000:2499] {
  set term png
  set output 'images/im_'.i.'.png'
  unset key                               
  set yrange [0:0.10]                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines lw 1.8, 'potential.data' with lines lw 1
}
