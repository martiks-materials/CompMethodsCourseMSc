
do for [i=0:9] {
  set term png
  set output 'images/im_00'.i.'.png'
  unset key                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines
}
do for [i=10:99] {
  set term png
  set output 'images/im_0'.i.'.png'
  unset key                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines
}
do for [i=100:999] {
  set term png
  set output 'images/im_'.i.'.png'
  unset key                               
  plot 'Data/psi_'.i.'.data' using 1:2 with lines
}
