reset

# Set terminal and output
set term post eps enh "Times-Italic" 14 #solid color
set output "../Fig/loads2D.eps";

# Set line style
set line style 1 lt 1 lw 2

# Set margin
NX=2; NY=2;
DX=0.1; DY=0.1; DX2 = 0.02; DY2 = 0.02; SX=0.8; SY=0.8;
set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY;
set size SX*NX+DX*3.2,SY*NY+DY*2.0+DY2*2

set grid
set nokey;

set multiplot;                    
set size SX,SY;

# CL
set origin DX2+DX,DY+SY+2*DY2
set ylabel 'C_L'
set format x "" 
plot './Log/AerodynamicLoads.txt' u 1:2 w l notitle ls 1;

# CM
set origin 2*DX2+2.2*DX+SX,DY+SY+2*DY2
set ylabel 'C_{Mx}'
set format x ""
plot './Log/AerodynamicLoads.txt' u 1:4 w l notitle ls 1;

# CD
set origin DX2+DX,DY
set ylabel 'C_D'
set label 1 "t [s]" at graph 0.5, graph -0.1 center
set format x "%6.3f"
plot './Log/AerodynamicLoads.txt' u 1:3 w l notitle ls 1;

# CH
set origin 2*DX2+2.2*DX+SX,DY
set ylabel 'C_{H}'
set label 1 "t [s]" at graph 0.5, graph -0.1 center
set format x "%6.3f" 
plot './Log/AerodynamicLoads.txt' u 1:5 w l notitle ls 1;

unset multiplot 
