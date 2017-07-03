#!/usr/bin/python

# ==========================================================
#                *** MODULES ***
# ==========================================================
# Built-in module
import os

# ==========================================================
#                   *** MAIN PROGRAM ***
# ==========================================================
# Size data
f = open('./Log/CCp_up.txt', 'r')
Nup = 0
for line in f.readlines():
   Nup = Nup + 1
f.close()

f = open('./Log/CCp_lo.txt', 'r')
Nlow = 0
for line in f.readlines():
   Nlow = Nlow + 1
f.close()

# Initialization
xup   = [-1.0]*Nup
xNup  = [-1.0]*Nup
xlow  = [-1.0]*Nlow
xNlow = [-1.0]*Nlow
Cpup  = [-1.0]*Nup
Cplow = [-1.0]*Nlow

# Read data
f = open('./Log/CCp_up.txt', 'r')
k = 0
for line in f.readlines():
    
   xup[k]  = float( line.split()[0] )
   xNup[k] = float( line.split()[0] )
   Cpup[k] = float( line.split()[1] )
   k = k + 1
    
f.close()   

f = open('./Log/CCp_lo.txt', 'r')
k = 0
for line in f.readlines():
    
   xlow[k]  = float( line.split()[0] )
   xNlow[k] = float( line.split()[0] )
   Cplow[k] = float( line.split()[1] )
   k = k + 1
    
f.close()  

# Sort upper data
xNup.sort()
CpNup = [-1.0]*Nup
for k in range(0, Nup):
   CpNup[k] = Cpup[ xup.index(xNup[k]) ]

xNlow.sort()
CpNlow = [-1.0]*Nlow
for k in range(0, Nlow):
   CpNlow[k] = Cplow[ xlow.index(xNlow[k]) ]
   
# Write on file
f = open( 'CCp_up.txt', 'w' )
for k in range(0, Nup):
   str_line = str(xNup[k]) + " " + str(CpNup[k]) + "\n"
   f.write( str_line )
   
f.close()

f = open( 'CCp_lo.txt', 'w' )
for k in range(0, Nlow):
   str_line = str(xNlow[k]) + " " + str(CpNlow[k]) + "\n"
   f.write( str_line )
   
f.close()

# Generate .plt file
f = open('xCp2D.plt', 'w')

f.write('reset\n')
f.write('set term post eps enh "Times-Italic" 14 solid color\n') 
f.write('set output "../Fig/xCp2D.eps";\n')
f.write('set nokey;\n')
f.write('set line style 1 lt 1 lw 2 pt 9 lc rgb "red"\n')
f.write('set line style 2 lt 1 lw 2 pt 11 lc rgb "blue"\n')
f.write('set grid\n')
f.write('set xlabel "x/c"\n')
f.write('set ylabel "{/Symbol -}C_p"\n')
f.write('set key right;\n')

f.write('plot "CCp_up.txt" w lp ls 1 title "Upper Surface", \\\n')
f.write('     "CCp_lo.txt" w lp ls 2 title "Lower Surface"\n')
     
f.close()

# Plot
os.system('gnuplot xCp2D.plt')
os.system('rm *.txt')
os.system('rm xCp2D.plt')
