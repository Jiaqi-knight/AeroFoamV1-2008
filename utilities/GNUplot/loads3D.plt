reset

# Set terminal and output
set term post eps enh "Times-Italic" 14 #enhanced color 
set output "../Fig/loads3D.eps";

# Set line style
set palette gray
set line style 1 lt 1 lw 2

# Set margin
NX=2; NY=3;
DX=0.15; DY=0.1; DY2 = 0.02; SX=0.8; SY=0.6;
set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY;
set size SX*NX+DX*3.2,SY*NY+DY*2.0+2*DY2

set grid
set nokey;

set multiplot;                    
set size SX,SY;

# CFX
set origin DX,DY+SY*2+2*DY2
set ylabel 'C_{F,X}'
set format x ""
plot './Log/AerodynamicLoads.txt' u 1:2 w l notitle ls 1;

# CFY
set origin DX,DY+SY+DY2
set ylabel 'C_{F,Y}'
set format x ""
plot './Log/AerodynamicLoads.txt' u 1:3 w l notitle ls 1;

# CMX
set origin 2.2*DX+SX,DY+SY*2+2*DY2
set ylabel 'C_{M,X}'
set format x "" 
plot './Log/AerodynamicLoads.txt' u 1:5 w l notitle ls 1;

# CMY
set origin 2.2*DX+SX,DY+SY+DY2
set ylabel 'C_{M,Y}'
set format x "" 
plot './Log/AerodynamicLoads.txt' u 1:6 w l notitle ls 1;

# CFZ
set origin DX,DY
set ylabel 'C_{F,Z}'
set format x "%6.3f"  
set label 1 "t [s]" at graph 0.5, graph -0.1 center
plot './Log/AerodynamicLoads.txt' u 1:4 w l notitle ls 1;

# CMZ
set origin 2.2*DX+SX,DY
set ylabel 'C_{M,Z}'
set format x "%6.3f" 
set label 1 "t [s]" at graph 0.5, graph -0.1 center
plot './Log/AerodynamicLoads.txt' u 1:7 w l notitle ls 1;

unset multiplot 
