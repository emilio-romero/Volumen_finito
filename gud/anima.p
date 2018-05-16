set terminal gif animate delay 6 
set output 'test.gif'
stats 'datos.dat' nooutput
set xrange[-1.2:1.2]
set yrange[-15.1:15.1]
N=int(STATS_blocks)
do for[i=1:N]{
  plot 'datos.dat' index (i-1) with lp
}
