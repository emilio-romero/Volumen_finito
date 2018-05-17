set terminal gif animate delay 5 
set output 'test3.gif'
stats 'datos2d.dat' nooutput
set xrange[-5.0:5.0]
set yrange[-5.0:5.0]
set zrange[-1.0:16.0]
N=int(STATS_blocks)
set dgrid3d
set hidden3d
set view equal xyz
do for[i=1:N]{
  splot 'datos2d.dat' index (i-1) u 1:2:3 with lp
}
