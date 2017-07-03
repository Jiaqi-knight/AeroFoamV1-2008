reset

# Set terminal and output
set term post eps enh "Times-Italic" 14 solid color 
set output "../Fig/modal3D.eps";

# Set line style
set line style 1 lt 1 lw 3 lc rgb "#000000"
set line style 2 lt 1 lw 3 lc rgb "#646464"
set line style 3 lt 1 lw 3 lc rgb "#969696"
set line style 4 lt 1 lw 3 lc rgb "#C7C7C7"

# Set margin
NX=1; NY=2;
DX=0.1; DY=0.1; DY2 = 0.02; SX=1.0; SY=0.6;
set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY;
set size SX*NX+DX*2.5,SY*NY+DY*2.0+DY2

set grid
set key left;

set multiplot;                    
set size SX,SY;
# q
set origin 1.5*DX,DY+SY+DY2
set ylabel "q [m]"
set format x "" 
plot './Log/qqs.txt' u 1:2 w l title " q_1" ls 1, \
     './Log/qqs.txt' u 1:3 w l title " q_2" ls 2, \
     './Log/qqs.txt' u 1:4 w l title " q_3" ls 3, \
     './Log/qqs.txt' u 1:5 w l title " q_4" ls 4
# Q
set origin 1.5*DX,DY
set ylabel "Q [m^2]"
set format x "%6.3f" 
set label 1 "{/Symbol t} [{/Symbol -}]" at graph 0.5, graph -0.1 center
plot './Log/QQa.txt' u 1:2 w l title " Q_1" ls 1, \
     './Log/QQa.txt' u 1:3 w l title " Q_2" ls 2, \
     './Log/QQa.txt' u 1:4 w l title " Q_3" ls 3, \
     './Log/QQa.txt' u 1:5 w l title " Q_4" ls 4
unset multiplot 
