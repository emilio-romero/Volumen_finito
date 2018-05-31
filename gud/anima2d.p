reset 
set terminal gif animate delay 5 
set output 'test3.gif'
stats 'datos2d.dat' nooutput
set xrange[-5:5]
set yrange[-5:5]
set zrange[-5:5]
set view equal xyz
N=int(STATS_blocks)
set dgrid3d 
set pm3d at b
set ticslevel 0.8
set isosample 40,40
set palette rgbformulae 22,13,-31
set cbrange[-2:3]
do for[i=1:N]{
  splot 'datos2d.dat' index (i-1) using 1:2:3 w l
}
