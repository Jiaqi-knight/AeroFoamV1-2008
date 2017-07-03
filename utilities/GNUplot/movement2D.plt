reset

# Set terminal and output
set term post eps enh "Times-Italic" 12
set output "../Fig/movement2D.eps";

# Set line style
set palette gray
set line style 1 lt -1 lw 2

# Set margin
NX=2; NY=3;
DX=0.15; DY=0.1; DY2 = 0.02; SX=0.8; SY=0.6;
set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY;
set size SX*NX+DX*3.2,SY*NY+DY*2.0+2*DY2

set grid
set nokey;

set multiplot;                    
set size SX,SY;

# alpha
set origin DX,DY+SY*2+2*DY2
set ylabel "{/Symbol a} [{/Symbol \260}]";
set format y "%4.2f"
set format x ""
plot './Log/Displacements.txt' u 1:3 w l notitle ls 1;

# h
set origin DX,DY+SY+DY2
set ylabel 'h [m]'
set format y "%4.2f"
set format x ""
plot './Log/Displacements.txt' u 1:2 w l notitle ls 1;

# alpha dot
set origin 2.2*DX+SX,DY+SY*2+2*DY2
set ylabel 'd{/Symbol a}/dt [{/Symbol \260}/s]'
set format y "%4.1f"
set format x ""  
plot './Log/Displacements.txt' u 1:6 w l notitle ls 1;

# h dot
set origin 2.2*DX+SX,DY+SY+DY2
set ylabel 'dh/dt [m/s]'
set format y "%4.1f"
set format x "" 
plot './Log/Displacements.txt' u 1:5 w l notitle ls 1;

# delta aileron
set origin DX,DY
set ylabel "{/Symbol d}_A [{/Symbol \260}]";
set format y "%4.2f"
set format x "%6.3f"
set label 1 "t [s]" at graph 0.5, graph -0.1 center
plot './Log/Displacements.txt' u 1:4 w l notitle ls 1;

# delta dot aileron
set origin 2.2*DX+SX,DY
set ylabel "d{/Symbol d}_A/dt [{/Symbol \260}/s]";
set format y "%4.1f"
set format x "%6.3f"
set label 1 "t [s]" at graph 0.5, graph -0.1 center
plot './Log/Displacements.txt' u 1:7 w l notitle ls 1;

unset multiplot 
