set terminal gif animate delay 2
set output 'test2.gif'
stats filename nooutput
set xrange[-3.0:3.0]
set yrange[-1.0:1.5]
N=int(STATS_blocks)
do for[i=1:N]{
  plot filename index (i-1) with lp
}
