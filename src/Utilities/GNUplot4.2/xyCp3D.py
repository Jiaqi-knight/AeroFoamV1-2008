#!/usr/bin/python

# ==========================================================
# ==========================================================
#                     *** INPUT ***
# ==========================================================
# ==========================================================

# Data related to read file
spanAxis = 'y';       # Right wing along spanAxis ['y'; 'z']
zsym     = 0.0        # Set upper and lower wing surface

# Interpolation mesh to plot Cp
xDistribution = 'c'   # Points distribution along x: 'l'=linear; 'c'=cosine
scaleAirfoil  = 2.5   # Scale factor of airfoil width 
meshDeform    = -10.0 # Plot deformed mesh: [0; K] with K scale factor of the deformation
                      # If negative the mesh of section Cp will not be deformed
Npx    = 40;          # Number of point along x-direction
Npy    = 15;          # Number of point along y-direction
Nlsi   = 8;           # Number of least squares interpolation points
fix_Cp = 0.0;         # Set LE Cp to the value fix_Cp if fix_Cp > 0.0

# Span sections
SpanSection = (0.01, 0.25, 0.50, 0.75, 0.99)   # % of span
expData     = 0
# ONERA DATA:
#SpanSection = (0.2 0.44 0.65 0.8 0.9 0.95 0.99)   # % of span
#expData     = 0       # WARNING: Check the experimental data and the plot
#nameExp     = 'OM6_340k_M084_a3' 

# ==========================================================
# ==========================================================
#                      *** MODULES ***
# ==========================================================
# ==========================================================
# Built-in module
import os
# Mathematical module
import math

# ==========================================================
# ==========================================================
#                *** FUNCTION DEFINITION ***
# ==========================================================
# ==========================================================
# ------------------------------
# *** MATMUL ***
# ------------------------------
# mat1 [m1 x n1]; mat2 [m2 x n2]
def matmul( mat1, mat2, m1, n1, m2, n2 ):
   
   if ( n1 != m2 ):
      print '*** ERROR in matmul: wrong matrix dimensions!'
      return 0
      
   mat12 = [0.0]*m1*n2
   for ki1 in range(0, m1):
      for ki2 in range(0, n2):
         for kj in range(0, n1):
            
            aa = mat1[kj + ki1*n1]
            bb = mat2[ki2 + kj*n2]
            mat12[ki2 + ki1*n2] = mat12[ki2 + ki1*n2] + aa*bb
            
   return mat12      

# ------------------------------
# *** FIND_NEAREST ***
# ------------------------------
def find_nearest( xx, yy, fxy, x0, y0, N ):
   
   # Compute the distance
   nn = len(xx)
   dd = [None]*nn
   for k in range(0, nn):
      dd[k] = math.sqrt( ( xx[k] - x0 )**2 + ( yy[k] - y0 )**2 )
   
   # Find the nearest
   vec = [None]*N*3
   for k in range(0, N):
      vmin = min(dd)
      id_v = dd.index(vmin)
      vec[k]       =  xx[id_v]
      vec[k + N]   =  yy[id_v]
      vec[k + 2*N] = fxy[id_v]
      dd[id_v]     = 100.0
   
   return vec

# ------------------------------
# *** TRANSPOSE ***
# ------------------------------
def transpose( A, m, n ):
   
   AT = [None]*n*m
   for ki in range(0, m):
      for kj in range(0, n):
         AT[ki + kj*m] = A[kj + ki*n]
   
   return AT

# ------------------------------
# *** INVERSE_3 ***
# ------------------------------
def inverse_3(A):
   
   if len(A) != 9:
      print '*** ERROR in inverse_3: wrong matrix dimensions!'
      return 0
   
   # Input
   a11 = A[0]; a12 = A[1]; a13 = A[2]
   a21 = A[3]; a22 = A[4]; a23 = A[5]
   a31 = A[6]; a32 = A[7]; a33 = A[8]
   
   # Determinant
   deta = a11*( a22*a33 - a32*a23 ) - \
          a12*( a21*a33 - a31*a23 ) + \
          a13*( a21*a32 - a31*a22 )
   
   if abs(deta) < 1e-15:
      print '*** ERROR in inverse_3: singular matrix!'
      deta = 1e-5
   
   # Inverse matrix
   A_1 = [None]*3*3
   A_1[0] =    a22*a33 - a23*a32;   A_1[1] = -( a12*a33 - a13*a32 ); A_1[2] =    a12*a23 - a13*a22
   A_1[3] = -( a21*a33 - a23*a31 ); A_1[4] =    a11*a33 - a13*a31;   A_1[5] = -( a11*a23 - a13*a21 )
   A_1[6] =    a21*a32 - a22*a31;   A_1[7] = -( a11*a32 - a12*a31 ); A_1[8] =    a11*a22 - a12*a21
   
   A_1 = [x/deta for x in A_1]
   
   return A_1

# ------------------------------
# *** INTERPOLATION_3 ***
# ------------------------------
# f(x,y) = a0 + a1*x + a2*y
def interpolation_3( xP, yP, xor, yor, f_xyor, N ):

   # Number of coefficient
   Ncoeff = 3
   
   # Find near points
   xyNearP = find_nearest( xor, yor, f_xyor, xP, yP, N )
   
   # Build A matrix [N x Ncoeff]
   A   = [None]*N*Ncoeff
   RHS = [None]*N 
   for ki in range(0, N):
      xi = xyNearP[ki]
      yi = xyNearP[ki + N]
      fi = xyNearP[ki + 2*N]
      A[ki*Ncoeff:(ki + 1)*Ncoeff] = (1.0, xi, yi) 
      RHS[ki] = fi
   
   # Compute A transpose [Ncoeff x N]
   AT = transpose( A, N, Ncoeff )
   
   # Invert [AT]*[A]
   AT_A   = matmul( AT, A, Ncoeff, N, N, Ncoeff )
   AT_A_1 = inverse_3( AT_A )
   
   # Compute coefficient (least squares)
   AT_RHS = matmul( AT, RHS, Ncoeff, N, N, 1 )
   X = matmul( AT_A_1, AT_RHS, Ncoeff, Ncoeff, Ncoeff, 1 )
   
   # Compute interpolated value
   Af = (1.0, xP, yP)
   tmp = matmul( Af, X, 1, Ncoeff, Ncoeff, 1 )
   f_xyP = tmp[0]

   return f_xyP
   
# ------------------------------
# *** INTERPOLATION_2 ***
# ------------------------------
# y = y0 + (y1 - y0)/(x1 - x0)*(x - x0) -> y = a0 + a1*x
def interpolation_2( xP, xor, yor ):

   Nor = len( xor )

   # Find right point
   flag = 1
   diff = [ x-xP for x in xor]
   absdiff = [abs(x) for x in diff]
   k = 0
   while flag == 1:
      dmin = min( absdiff )
      ind = absdiff.index( dmin )
      absdiff[ind] = 100.0
      if diff[ind] >= 0:
         flag = 0
         x1 = xor[ind]; y1 = yor[ind]
      if k >= Nor:
         flag = 0
      k = k + 1
   
   # Find left point
   flag = 1
   diff = [ x-xP for x in xor]
   absdiff = [abs(x) for x in diff]
   k = 0
   while flag == 1:
      dmin = min( absdiff )
      ind = absdiff.index( dmin )
      absdiff[ind] = 100.0
      if diff[ind] < 0:
         flag = 0
         x0 = xor[ind]; y0 = yor[ind]
      if k >= Nor:
         flag = 0
      k = k + 1   
   
   # Out of the range
   xmin = min(xor)
   xmax = max(xor)
   if xP <= xmin: 
      x0 = xmin
      y0 = yor[xor.index(xmin)]
   if xP >= xmax: 
      x1 = xmax
      y1 = yor[xor.index(xmax)]

   # Linear interpolation
   yP = y0 + (y1 - y0)/(x1 - x0 + 1e-15)*(xP - x0)
   
   return yP

# ==========================================================
# ==========================================================
#                    *** MAIN PROGRAM ***
# ==========================================================
# ==========================================================
# Assumption:
# x-axis oriented as the stream flow
# y-axis oriented as the right wing
# z-axis according to the right hand rule

print '-----------------------------------------'
print '  *** 3D Linear interpolation of Cp *** '
print ' Reading data ...'
# -------------------------------------------------------
# *** READ CP FILE ***
# -------------------------------------------------------
# Data size
f = open('../Log/CCp.txt', 'r')
N_n = 0
for line in f.readlines():
   N_n = N_n + 1
f.close()

# Initialization
xx = [None]*N_n
yy = [None]*N_n
zz = [None]*N_n
Cp = [None]*N_n

# Read data
f = open('../Log/CCp.txt', 'r')
k = 0
for line in f.readlines():
    
   xx[k] = float( line.split()[0] )
   Cp[k] = float( line.split()[3] )
   if spanAxis == 'y':
      yy[k] = float( line.split()[1] )
      zz[k] = float( line.split()[2] )
   elif spanAxis == 'z':
      zz[k] = float( line.split()[1] )
      yy[k] = float( line.split()[2] )
   k = k + 1
    
f.close()    

# Wing span
ymin = min( yy )
ymax = max( yy )

# Count upper data
Nup = 0
for k in range(0,N_n):

   if zz[k] >= zsym:
      Nup = Nup + 1  
       
Nlow = N_n - Nup
   
# Initialization
xxup  = [None]*Nup;  yyup  = [None]*Nup;  zzup  = [None]*Nup;  Cpup  = [None]*Nup
xxlow = [None]*Nlow; yylow = [None]*Nlow; zzlow = [None]*Nlow; Cplow = [None]*Nlow

# Upper and lower data
idup  = 0; idlow = 0
Nupr = 0; Nlowr = 0
Nupt = 0; Nlowt = 0
for k in range(0,N_n):
    
   if zz[k] >= zsym:
      # Upper surface
      xxup[idup] = xx[k]
      yyup[idup] = yy[k]
      zzup[idup] = zz[k]
      Cpup[idup] = Cp[k]
      idup = idup + 1
      
   else:
      # Lower surface
      xxlow[idlow] = xx[k]
      yylow[idlow] = yy[k]
      zzlow[idlow] = zz[k]
      Cplow[idlow] = Cp[k]
      idlow = idlow + 1
      
   # Count root points
   if abs( yy[k] - ymin ) <= 1.0e-8:
      if zz[k] >= zsym:
         Nupr  = Nupr + 1
      if zz[k] <= zsym:
         Nlowr = Nlowr + 1
      
   # Count tip points
   if abs( yy[k] - ymax ) <= 1.0e-8:
      if zz[k] >= zsym:
         Nupt  = Nupt + 1
      if zz[k] <= zsym:
         Nlowt = Nlowt + 1
   
# Root chord
Xuproot  = [None]*Nupr
Zuproot  = [None]*Nupr
Xlowroot = [None]*Nlowr
Zlowroot = [None]*Nlowr
# Tip chord
Xuptip  = [None]*Nupt
Zuptip  = [None]*Nupt
Xlowtip = [None]*Nlowt
Zlowtip = [None]*Nlowt

idupr = 0; idlowr = 0
idupt = 0; idlowt = 0
for k in range(0, N_n):

   # Find root points
   if abs( yy[k] - ymin ) <= 1.0e-8:
      if zz[k] >= zsym:
         Xuproot[idupr] = xx[k]
         Zuproot[idupr] = zz[k]
         idupr = idupr + 1
      if zz[k] <= zsym:
         Xlowroot[idlowr] = xx[k]
         Zlowroot[idlowr] = zz[k]
         idlowr = idlowr + 1
      
   # Find tip points
   if abs( yy[k] - ymax ) <= 1.0e-8:
      if zz[k] >= zsym:
         Xuptip[idupt] = xx[k]
         Zuptip[idupt] = zz[k]
         idupt = idupt + 1
      if zz[k] <= zsym:
         Xlowtip[idlowt] = xx[k]
         Zlowtip[idlowt] = zz[k]
         idlowt = idlowt + 1

xrLE = min( Xuproot )     
xrTE = max( Xuproot )  
xtLE = min( Xuptip )     
xtTE = max( Xuptip )
cr = xrTE - xrLE
ct = xtTE - xtLE

# Scale
Zuproot  = [ scaleAirfoil*x for x in Zuproot ]
Zlowroot = [ scaleAirfoil*x for x in Zlowroot ]
Zuptip  = [ scaleAirfoil*x for x in Zuptip ]
Zlowtip = [ scaleAirfoil*x for x in Zlowtip ]

# -------------------------------------------------------
# *** READ DEFORMATION FILE ***
# ------------------------------------------------------
if abs(meshDeform) > 0.0:
   print ' > Plot mesh deformation amplified by K = ' + str(abs(meshDeform))
f = open('../Log/StructuralModel.deform', 'r')
N_n = 0
id_loc = 1
for line in f.readlines():
   if id_loc == 4:
      N_n = N_n + 1
      id_loc = id_loc + 1
   elif id_loc == 2:
      N_n = N_n + 1
      id_loc = id_loc + 1
   elif id_loc == 1:
      N_n = N_n + 1
      id_loc = id_loc + 1
   elif id_loc == 7:
      id_loc = 1
   else:
      id_loc = id_loc + 1
f.close()

# Initialization
N_n = N_n/4
xU = [100.0]*N_n
yU = [100.0]*N_n
zU = [100.0]*N_n

f = open('../Log/StructuralModel.deform', 'r')
id_loc = 1
id_U = 0
for line in f.readlines():
   # Read point
   xtmp = -100; ytmp = -100; ztmp = -100;
   if id_loc == 4:
      xtmp = float( line.split()[0] )
      ytmp = float( line.split()[1] )
      ztmp = float( line.split()[2] )
      id_loc = id_loc + 1
   elif id_loc == 2:
      xtmp = float( line.split()[0] )
      ytmp = float( line.split()[1] )
      ztmp = float( line.split()[2] )
      id_loc = id_loc + 1
   elif id_loc == 1:
      xtmp = float( line.split()[0] )
      ytmp = float( line.split()[1] )
      ztmp = float( line.split()[2] )
      id_loc = id_loc + 1
   elif id_loc == 7:
      id_loc = 1
   else:
      id_loc = id_loc + 1
   
   # Add point only if it is new
   xdiff = [ abs(x - xtmp) for x in xU ]
   ydiff = [ abs(y - ytmp) for y in yU ]
   zdiff = [ abs(z - ztmp) for z in zU ]
   xdm = min(xdiff); ydm = min(ydiff); zdm = min(zdiff)
   if ( xdm > 0.0 )or( ydm > 0.0 )or( zdm > 0.0 ):
      xU[id_U] = xtmp
      yU[id_U] = ytmp
      zU[id_U] = ztmp
      id_U = id_U + 1
   
f.close()
N_n = id_U - 1
xU = xU[:N_n]
yU = yU[:N_n]
zU = zU[:N_n]

# -------------------------------------------------------
# *** WING MESH ***
# -------------------------------------------------------
# Build Mesh interpolation yMesh
My = 40; Mx = 40;
print ' Building wing surface   [' + str(Mx) + 'x' + str(My) + '] ...'
Dy = ( ymax - ymin )/( My - 1.0 )
yMesh = [None]*Mx*My
yvec = range(0, My)
yvec = [Dy*x for x in yvec]
yMesh[0:My] = yvec
for k in range(0, Mx):
   yMesh[k*My:(k+1)*My] = yvec

# Build Mesh interpolation xMesh
xMesh = [None]*Mx*My
for k in range(0, My):
   
   xTE = xrTE + ( xtTE - xrTE )/( ymax - ymin )*( yMesh[k] - ymin )
   xLE = xrLE + ( xtLE - xrLE )/( ymax - ymin )*( yMesh[k] - ymin )
   if xDistribution == 'l':
      # Linear distribution
      Dx = ( xTE - xLE )/( Mx - 1.0 )
      xvec = [Dx*kk for kk in range(0, Mx)]
      
   elif xDistribution == 'c':
      # Cosine distribution
      cloc = xTE - xLE
      Dx = ( xTE - xLE )/( Mx - 1.0 )
      xtmp = [ Dx*kk/cloc for kk in range(0, Mx) ]
      xvec = [ 0.5*(1.0 - math.cos(math.pi*x))*cloc for x in xtmp ]

   for kk in range(0, Mx):
      xMesh[k+kk*My] = xvec[kk] + xLE

# Build Mesh interpolation zMesh
zupMesh  = [None]*Mx*My
zlowMesh = [None]*Mx*My
for ky in range(0, My):

   xTE = xrTE + ( xtTE - xrTE )/( ymax - ymin )*( yMesh[ky] - ymin )
   xLE = xrLE + ( xtLE - xrLE )/( ymax - ymin )*( yMesh[ky] - ymin )
   c = xTE - xLE
   for kx in range(0, Mx):
      
      xP = xMesh[ky+kx*My]
      xP0 = ( xP - xLE )/c*cr + xrLE
      xP1 = ( xP - xLE )/c*ct + xtLE
      # Upper surface
      zP0up = interpolation_2( xP0, Xuproot, Zuproot )
      zP1up = interpolation_2( xP1, Xuptip,  Zuptip )
      zupMesh[ky+kx*My] = zP0up + ( zP1up - zP0up )/( ymax - ymin )*( yMesh[ky] - ymin )
      # Lower surface
      zP0low = interpolation_2( xP0, Xlowroot, Zlowroot )
      zP1low = interpolation_2( xP1, Xlowtip,  Zlowtip )
      zlowMesh[ky+kx*My] = zP0low + ( zP1low - zP0low )/( ymax - ymin )*( yMesh[ky] - ymin )

# Scale
Zuproot  = [ x/scaleAirfoil for x in Zuproot ]
Zlowroot = [ x/scaleAirfoil for x in Zlowroot ]
Zuptip  = [ x/scaleAirfoil for x in Zuptip ]
Zlowtip = [ x/scaleAirfoil for x in Zlowtip ]

# Mesh deformation
zupMeshDef  = zupMesh[:]
zlowMeshDef = zlowMesh[:]
if abs(meshDeform) > 0.0:
   for k in range(0, Mx*My):
      tmp = interpolation_3( xMesh[k], yMesh[k], xU, yU, zU, Nlsi )
      zupMeshDef[k]  = zupMeshDef[k]  + tmp*abs(meshDeform)
      zlowMeshDef[k] = zlowMeshDef[k] + tmp*abs(meshDeform)
      if meshDeform > 0.0:
         zupMesh[k]  = zupMesh[k]  + tmp*meshDeform
         zlowMesh[k] = zlowMesh[k] + tmp*meshDeform
 
# -------------------------------------------------------
# *** CP INTERPOLATION ON SECTIONS ***
# -------------------------------------------------------
# Build matrix interpolation yMatrix
Nsec = len( SpanSection )
print ' Creating sections data  [' + str(Npx) + 'x' + str(Nsec) + '] ...'
ySection = [None]*Npx*Nsec
yvec = [ ( ymax - ymin )*x for x in SpanSection ]
ySection[0:Nsec] = yvec
for k in range(0, Npx):
   ySection[k*Nsec:(k+1)*Nsec] = yvec

# Build matrix interpolation xSection
xSection = [None]*Npx*Nsec
for k in range(0, Nsec):
   
   xTE = xrTE + ( xtTE - xrTE )/( ymax - ymin )*( ySection[k] - ymin )
   xLE = xrLE + ( xtLE - xrLE )/( ymax - ymin )*( ySection[k] - ymin )
   if xDistribution == 'l':
      # Linear distribution
      Dx = ( xTE - xLE )/( Npx - 1.0 )
      xvec = [Dx*kk for kk in range(0, Npx)]
      
   elif xDistribution == 'c':
      # Cosine distribution
      cloc = xTE - xLE
      Dx = ( xTE - xLE )/( Npx - 1.0 )
      xtmp = [ Dx*kk/cloc for kk in range(0, Npx) ]
      xvec = [ 0.5*(1.0 - math.cos(math.pi*x))*cloc for x in xtmp ]

   for kk in range(0, Npx):
      xSection[k+kk*Nsec] = xvec[kk] + xLE

# Build matrix interpolation zSection
zupSection  = [None]*Npx*Nsec
zlowSection = [None]*Npx*Nsec
for ky in range(0, Nsec):

   xTE = xrTE + ( xtTE - xrTE )/( ymax - ymin )*( ySection[ky] - ymin )
   xLE = xrLE + ( xtLE - xrLE )/( ymax - ymin )*( ySection[ky] - ymin )
   c = xTE - xLE
   for kx in range(0, Npx):
      
      xP = xSection[ky+kx*Nsec]
      xP0 = ( xP - xLE )/c*cr + xrLE
      xP1 = ( xP - xLE )/c*ct + xtLE
      # Upper surface
      zP0up = interpolation_2( xP0, Xuproot, Zuproot )
      zP1up = interpolation_2( xP1, Xuptip,  Zuptip )
      zupSection[ky+kx*Nsec] = zP0up + ( zP1up - zP0up )/( ymax - ymin )*( ySection[ky] - ymin )
      # Lower surface
      zP0low = interpolation_2( xP0, Xlowroot, Zlowroot )
      zP1low = interpolation_2( xP1, Xlowtip,  Zlowtip )
      zlowSection[ky+kx*Nsec] = zP0low + ( zP1low - zP0low )/( ymax - ymin )*( ySection[ky] - ymin )
      #zlowSection[ky+kx*Nsec] = zupSection[ky+kx*Nsec]

# Mesh deformation
if meshDeform > 0.0:
   for k in range(0, Npx*Nsec):
      tmp = interpolation_3( xSection[k], ySection[k], xU, yU, zU, Nlsi )
      zupSection[k]  = zupSection[k]  + tmp*meshDeform
      zlowSection[k] = zlowSection[k] + tmp*meshDeform

# Build matrix interpolation CpSection
CpupSection  = [None]*Npx*Nsec
CplowSection = [None]*Npx*Nsec
for k in range(0, Npx*Nsec):
   
   CpupSection[k]  = interpolation_3( xSection[k], ySection[k], xxup, yyup, Cpup, Nlsi )
   CplowSection[k] = interpolation_3( xSection[k], ySection[k], xxlow, yylow, Cplow, Nlsi )

# Fix the minimun Cp of each section
if fix_Cp > 0.0:
   for k in range(0, Nsec):
      CpupSection[k]  = fix_Cp
      CplowSection[k] = fix_Cp

# -------------------------------------------------------
# *** MESH OFFSET TO MAXIMUN CP ***
# -------------------------------------------------------
zoffsetup  = max( CpupSection )
zoffsetlow = max( CplowSection )
zoffset = zoffsetup
if zoffsetlow > zoffset:
   zoffset = zoffsetlow
zupMesh     = [ x - zoffset for x in zupMesh ]
zlowMesh    = [ x - zoffset for x in zlowMesh ]
zupMeshDef  = [ x - zoffset for x in zupMeshDef ]
zlowMeshDef = [ x - zoffset for x in zlowMeshDef ]
zupSection  = [ x - zoffset for x in zupSection ]
zlowSection = [ x - zoffset for x in zlowSection ]

# -------------------------------------------------------
# *** WRITE ON FILE ***
# -------------------------------------------------------
print ' Creating temporary files ...'

# Write mesh on file
fup  = open( 'mesh_up.txt', 'w' )
flow = open( 'mesh_low.txt', 'w' )
for ky in range(0,My-1):
   
   for kx in range(0,Mx-1):
      
      id1 = ky +     (kx    )*My
      id2 = ky +     (kx + 1)*My
      id3 = ky + 1 + (kx + 1)*My
      id4 = ky + 1 + (kx    )*My
      
      # p1;p2; <space> p4;p3; <2x space>
      # Upper surface
      str_line = str(  xMesh[id1]) + " " + str(yMesh[id1]) + " " + \
                 str(zupMesh[id1]) + " " + str(zupMesh[id1]) + "\n"
      fup.write(str_line)
      str_line = str(  xMesh[id2]) + " " + str(yMesh[id2]) + " " + \
                 str(zupMesh[id2]) + " " + str(zupMesh[id2]) + "\n\n"
      fup.write(str_line)
      str_line = str(  xMesh[id4]) + " " + str(yMesh[id4]) + " " + \
                 str(zupMesh[id4]) + " " + str(zupMesh[id4]) + "\n"
      fup.write(str_line)
      str_line = str(  xMesh[id3]) + " " + str(yMesh[id3]) + " " + \
                 str(zupMesh[id3]) + " " + str(zupMesh[id3]) + "\n\n\n"
      fup.write(str_line)
      
      # Lower surface
      str_line = str(   xMesh[id1]) + " " + str(yMesh[id1]) + " " + \
                 str(zlowMesh[id1]) + " " + str(zlowMesh[id1]) + "\n"
      flow.write(str_line)
      str_line = str(   xMesh[id2]) + " " + str(yMesh[id2]) + " " + \
                 str(zlowMesh[id2]) + " " + str(zlowMesh[id2]) + "\n\n"
      flow.write(str_line)
      str_line = str(   xMesh[id4]) + " " + str(yMesh[id4]) + " " + \
                 str(zlowMesh[id4]) + " " + str(zlowMesh[id4]) + "\n"
      flow.write(str_line)
      str_line = str(   xMesh[id3]) + " " + str(yMesh[id3]) + " " + \
                 str(zlowMesh[id3]) + " " + str(zlowMesh[id3]) + "\n\n\n"
      flow.write(str_line)

fup.close()
flow.close()

# Write deformed mesh on file
fup  = open( 'meshDef_up.txt', 'w' )
flow = open( 'meshDef_low.txt', 'w' )
for ky in range(0,My-1):
   
   for kx in range(0,Mx-1):
      
      id1 = ky +     (kx    )*My
      id2 = ky +     (kx + 1)*My
      id3 = ky + 1 + (kx + 1)*My
      id4 = ky + 1 + (kx    )*My
      
      # p1;p2; <space> p4;p3; <2x space>
      # Upper surface
      str_line = str(     xMesh[id1]) + " " + str(     yMesh[id1]) + " " + \
                 str(zupMeshDef[id1]) + " " + str(zupMeshDef[id1]) + "\n"
      fup.write(str_line)
      str_line = str(     xMesh[id2]) + " " + str(     yMesh[id2]) + " " + \
                 str(zupMeshDef[id2]) + " " + str(zupMeshDef[id2]) + "\n\n"
      fup.write(str_line)
      str_line = str(     xMesh[id4]) + " " + str(     yMesh[id4]) + " " + \
                 str(zupMeshDef[id4]) + " " + str(zupMeshDef[id4]) + "\n"
      fup.write(str_line)
      str_line = str(     xMesh[id3]) + " " + str(     yMesh[id3]) + " " + \
                 str(zupMeshDef[id3]) + " " + str(zupMeshDef[id3]) + "\n\n\n"
      fup.write(str_line)
      
      # Lower surface
      str_line = str(      xMesh[id1]) + " " + str(     yMesh[id1]) + " " + \
                 str(zlowMeshDef[id1]) + " " + str(zlowMeshDef[id1]) + "\n"
      flow.write(str_line)
      str_line = str(      xMesh[id2]) + " " + str(     yMesh[id2]) + " " + \
                 str(zlowMeshDef[id2]) + " " + str(zlowMeshDef[id2]) + "\n\n"
      flow.write(str_line)
      str_line = str(      xMesh[id4]) + " " + str(     yMesh[id4]) + " " + \
                 str(zlowMeshDef[id4]) + " " + str(zlowMeshDef[id4]) + "\n"
      flow.write(str_line)
      str_line = str(      xMesh[id3]) + " " + str(     yMesh[id3]) + " " + \
                 str(zlowMeshDef[id3]) + " " + str(zlowMeshDef[id3]) + "\n\n\n"
      flow.write(str_line)

fup.close()
flow.close()

# Write sections of file
f = open( 'mesh_section.txt', 'w' )
for ky in range(0,Nsec-1):
   
   for kx in range(0,Npx-1):
      
      id1 = ky +     (kx    )*Nsec
      id2 = ky +     (kx + 1)*Nsec
      id3 = ky + 1 + (kx + 1)*Nsec
      id4 = ky + 1 + (kx    )*Nsec
      
      # p1;p2; <space> p4;p3; <2x space>
      # Upper surface
      str_line = str(  xSection[id1]) + " " + str(ySection[id1]) + " " + \
                 str(zupSection[id1]) + " 0.0\n"
      f.write(str_line)
      str_line = str(  xSection[id2]) + " " + str(ySection[id2]) + " " + \
                 str(zupSection[id2]) + " 0.0\n\n\n"
      f.write(str_line)
      str_line = str(  xSection[id4]) + " " + str(ySection[id4]) + " " + \
                 str(zupSection[id4]) + " 0.0\n"
      f.write(str_line)
      str_line = str(  xSection[id3]) + " " + str(ySection[id3]) + " " + \
                 str(zupSection[id3]) + " 0.0\n\n\n"
      f.write(str_line)
      
      # Lower surface
      #str_line = str(   xSection[id1]) + " " + str(ySection[id1]) + " " + \
      #           str(zlowSection[id1]) + " 0.0\n"
      #f.write(str_line)
      #str_line = str(   xSection[id2]) + " " + str(ySection[id2]) + " " + \
      #           str(zlowSection[id2]) + " 0.0\n\n"
      #f.write(str_line)
      #str_line = str(   xSection[id4]) + " " + str(ySection[id4]) + " " + \
      #           str(zlowSection[id4]) + " 0.0\n"
      #f.write(str_line)
      #str_line = str(   xSection[id3]) + " " + str(ySection[id3]) + " " + \
      #           str(zlowSection[id3]) + " 0.0\n\n\n"
      #f.write(str_line)

f.close()

# Write local frame of reference
f = open( 'sys_section.txt', 'w' )
for ky in range(0,Nsec):
      
   ind = ky
   #locCp = -CpupSection[ind]
   locCp = [abs(x) for x in Cp]
   locCp = max(locCp)
   
   str_line = str(xSection[ind]) + " " + str(ySection[ind]) + " " + str(zupSection[ind]) + " 0.0\n"
   f.write(str_line)
   str_line = str(xSection[ind]) + " " + str(ySection[ind]) + " " + \
              str(locCp + zupSection[ind] + zoffset) + " 0.0\n\n\n"
   f.write(str_line)
   ind = ky + (Npx - 1)*Nsec
   str_line = str(xSection[ind]) + " " + str(ySection[ind]) + " " + str(zupSection[ind]) + " 0.0\n"
   f.write(str_line)
   str_line = str(xSection[ind]) + " " + str(ySection[ind]) + " " + \
              str(-CpupSection[ind] + zupSection[ind] + zoffset) + " 0.0\n\n\n"
   f.write(str_line)
   
   # Horizontal line
   ind = ky
   str_line = str(xSection[ind]) + " " + str(ySection[ind]) + " " + \
              str(zupSection[ind] + zoffset) + " 0.0\n"
   f.write(str_line)
   ind = ky + (Npx - 1)*Nsec
   cloc = ( xSection[ind] - xSection[ky] )*1.15
   str_line = str(cloc+xSection[ky]) + " " + str(ySection[ind]) + " " + \
              str(zupSection[ind] + zoffset) + " 0.0\n\n\n"
   f.write(str_line)

f.close()

# 
f = open('wing_border.txt', 'w')
# LE
for k in range(0,My-1):
   str_line = str(xMesh[k]) + " " + str(yMesh[k]) + " " + str(zupMesh[k]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[k+1]) + " " + str(yMesh[k+1]) + " " + str(zupMesh[k+1]) + "\n\n\n"
   f.write(str_line)
# TE
k = 0
str_line = str(xMesh[k]) + " " + str(yMesh[k]) + " " + str(zupMesh[k]) + "\n"
f.write(str_line)
str_line = str(xMesh[k+My*(Mx-1)]) + " " + str(yMesh[k+My*(Mx-1)]) + " " + \
           str(zupMesh[k+My*(Mx-1)]) + "\n"
f.write(str_line)
for k in range(0,My-1):
   str_line = str(xMesh[k+My*(Mx-1)]) + " " + str(yMesh[k+My*(Mx-1)]) + " " + \
              str(zupMesh[k+My*(Mx-1)]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[k+1+My*(Mx-1)]) + " " + str(yMesh[k+1+My*(Mx-1)]) + " " + \
              str(zupMesh[k+1+My*(Mx-1)]) + "\n\n\n"
   f.write(str_line)
# Root chord
for k in range(0, Mx-1):
   str_line = str(xMesh[k*My]) + " " + str(yMesh[k*My]) + " " + str(zupMesh[k*My]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[(k+1)*My]) + " " + str(yMesh[(k+1)*My]) + " " + str(zupMesh[(k+1)*My]) + "\n\n\n"
   f.write(str_line)
   str_line = str(xMesh[k*My]) + " " + str(yMesh[k*My]) + " " + str(zlowMesh[k*My]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[(k+1)*My]) + " " + str(yMesh[(k+1)*My]) + " " + str(zlowMesh[(k+1)*My]) + "\n\n\n"
   f.write(str_line)
# Tip chord
for k in range(0, Mx-1):
   str_line = str(xMesh[(k+1)*My-1]) + " " + str(yMesh[(k+1)*My-1]) + " " + \
              str(zupMesh[(k+1)*My-1]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[(k+2)*My-1]) + " " + str(yMesh[(k+2)*My-1]) + " " + \
              str(zupMesh[(k+2)*My-1]) + "\n"
   f.write(str_line)
f.close()

# 
f = open('wingDef_border.txt', 'w')
# LE
for k in range(0,My-1):
   str_line = str(xMesh[k]) + " " + str(yMesh[k]) + " " + str(zupMeshDef[k]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[k+1]) + " " + str(yMesh[k+1]) + " " + str(zupMeshDef[k+1]) + "\n\n\n"
   f.write(str_line)
# TE
k = 0
str_line = str(xMesh[k]) + " " + str(yMesh[k]) + " " + str(zupMeshDef[k]) + "\n"
f.write(str_line)
str_line = str(xMesh[k+My*(Mx-1)]) + " " + str(yMesh[k+My*(Mx-1)]) + " " + \
           str(zupMeshDef[k+My*(Mx-1)]) + "\n"
f.write(str_line)
for k in range(0,My-1):
   str_line = str(xMesh[k+My*(Mx-1)]) + " " + str(yMesh[k+My*(Mx-1)]) + " " + \
              str(zupMeshDef[k+My*(Mx-1)]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[k+1+My*(Mx-1)]) + " " + str(yMesh[k+1+My*(Mx-1)]) + " " + \
              str(zupMeshDef[k+1+My*(Mx-1)]) + "\n\n\n"
   f.write(str_line)
# Root chord
for k in range(0, Mx-1):
   str_line = str(xMesh[k*My]) + " " + str(yMesh[k*My]) + " " + str(zupMeshDef[k*My]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[(k+1)*My]) + " " + str(yMesh[(k+1)*My]) + " " + str(zupMeshDef[(k+1)*My]) + "\n\n\n"
   f.write(str_line)
   str_line = str(xMesh[k*My]) + " " + str(yMesh[k*My]) + " " + str(zlowMesh[k*My]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[(k+1)*My]) + " " + str(yMesh[(k+1)*My]) + " " + str(zlowMesh[(k+1)*My]) + "\n\n\n"
   f.write(str_line)
# Tip chord
for k in range(0, Mx-1):
   str_line = str(xMesh[(k+1)*My-1]) + " " + str(yMesh[(k+1)*My-1]) + " " + \
              str(zupMeshDef[(k+1)*My-1]) + "\n"
   f.write(str_line)
   str_line = str(xMesh[(k+2)*My-1]) + " " + str(yMesh[(k+2)*My-1]) + " " + \
              str(zupMeshDef[(k+2)*My-1]) + "\n"
   f.write(str_line)
f.close()

# Wing sections Cp
for ky in range(0, Nsec):
   
   f1 = open('xCp_up' + str(ky+1) + '.txt', 'w')
   f2 = open('xCp_low' + str(ky+1) + '.txt', 'w')
   f3 = open('xfilled' + str(ky+1) + '.txt', 'w')
   for kx in range(0, Npx):
      
      ind = ky+kx*Nsec
      str_line = str(    xSection[ind]) + " " + str(ySection[ind]) + " " + \
                 str(-CpupSection[ind]+zupSection[ind]+zoffset) + "\n"
      f1.write(str_line)
      str_line = str(     xSection[ind]) + " " + str(ySection[ind]) + " " + \
                 str(-CplowSection[ind]+zupSection[ind]+zoffset) + "\n"
      f2.write(str_line)
      str_line = str(     xSection[ind]) + " " + str(ySection[ind]) + " " + \
                 str(-CplowSection[ind]+zupSection[ind]+zoffset) + "\n"
      f3.write(str_line)
      str_line = str(    xSection[ind]) + " " + str(ySection[ind]) + " " + \
                 str(-CpupSection[ind]+zupSection[ind]+zoffset) + "\n\n\n"
      f3.write(str_line)

   f1.close() 
   f2.close()
   f3.close() 

# -------------------------------------------------------
# *** WRITE SECTION FILE TO GNUPLOT ***
# -------------------------------------------------------
f = open('sec3D.plt', 'w')

f.write('reset\n')
f.write('set term X11 enh\n')
f.write('set grid\n')
f.write('set xlabel "x"\n')
f.write('set ylabel "y"\n')
f.write('set zlabel "{/Symbol -}C_p"\n')
f.write('set view 45, 345\n')
f.write('set size 1.0 ,1.0\n')
f.write('set autoscale\n')
#f.write('set palette gray\n')
f.write('set palette defined ( 0 0.9 0.9 0.9, 1 0.9 0.9 0.9 )\n')
#f.write('set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )\n')
f.write('set pm3d at s explicit\n')
f.write('unset colorbox\n\n')

f.write('set style line 1 lt 1 lw 2.0 lc rgb "blue"\n')
f.write('set style line 2 lt 1 lw 2.0 lc rgb "red"\n')
f.write('set style line 3 lt 1 lw 0.01 lc rgb "black"\n')
f.write('set style line 4 lt 2 lw 1.0 lc rgb "black"\n')
f.write('set style line 5 lt 1 lw 0.1 lc rgb "black"\n')
f.write('set style line 6 lt 1 lw 0.8 lc rgb "black"\n')
f.write('set style line 7 pt 7 ps 1 lc rgb "black"\n\n')

f.write('set style arrow 1 head filled size screen 0.008,20 ls 6\n\n')

for ky in range(0,Nsec):
      
   ind = ky
   #locCp = -CpupSection[ind]
   locCp = [abs(x) for x in Cp]
   locCp = max(locCp)
   P1 = str(xSection[ind]) + "," + str(ySection[ind]) + "," + str(zupSection[ind])
   P2 = str(xSection[ind]) + "," + str(ySection[ind]) + "," + str(locCp + zupSection[ind] + zoffset)
   f.write('set arrow from ' + P1 + ' to ' + P2 + ' as 1\n')
      
   # Horizontal line
   ind = ky
   P1 = str(xSection[ind]) + "," + str(ySection[ind]) + "," + str(zupSection[ind] + zoffset)
   ind = ky + (Npx - 1)*Nsec
   cloc = ( xSection[ind] - xSection[ky] )*1.15
   P2 = str(cloc+xSection[ky]) + "," + str(ySection[ind]) + "," + str(zupSection[ind] + zoffset)
   f.write('set arrow from ' + P1 + ' to ' + P2 + ' as 1\n')

f.write('\n')

f.write('splot "mesh_low.txt" using 1:2:3:4 w pm3d notitle, \\\n')
f.write('      "mesh_up.txt" using 1:2:3:4 w pm3d notitle, \\\n')
f.write('      "mesh_section.txt" using 1:2:3 w l ls 4 notitle, \\\n')
f.write('      "wing_border.txt" using 1:2:3 w l ls 5 notitle, \\\n')
f.write('      "sys_section.txt" using 1:2:3 w l ls 6 notitle, \\\n')

# Experimental data
if expData == 1:
   for k in range(0, Nsec):
      f.write('      "./Exp/gnu_' + nameExp + '_' + str(k+1) + '_up.txt" using 1:2:3 w p ls 6 notitle, \\\n')
      f.write('      "./Exp/gnu_' + nameExp + '_' + str(k+1) + '_low.txt" using 1:2:3 w p ls 6 notitle, \\\n')
######

for k in range(0, Nsec):
   f.write('      "xfilled' + str(k+1) + '.txt" using 1:2:3 w l ls 3 notitle, \\\n')
for k in range(0, Nsec-1):
   f.write('      "xCp_up' + str(k+1) + '.txt" using 1:2:3 w l ls 1 notitle, \\\n')
   f.write('      "xCp_low' + str(k+1) + '.txt" using 1:2:3 w l ls 2 notitle, \\\n')
f.write('      "xCp_up' + str(Nsec) + '.txt" using 1:2:3 w l ls 1 notitle, \\\n')
f.write('      "xCp_low' + str(Nsec) + '.txt" using 1:2:3 w l ls 2 notitle\n\n')

f.write('print "Press return to export current image into secCp3D.eps file...";\n')
f.write('pause -1\n\n')

f.write('set term post eps enh "Times-Italic" 10 color solid\n')
f.write('set output "../Fig/secCp3D.eps"\n')
f.write('replot\n')

f.close()

# -------------------------------------------------------
# *** WRITE DEFORMATION FILE TO GNUPLOT ***
# -------------------------------------------------------
f = open('def3D.plt', 'w')

f.write('reset\n')
f.write('set term X11 enh\n')
f.write('set grid\n')
f.write('set xlabel "x"\n')
f.write('set ylabel "y"\n')
f.write('set zlabel "z"\n')
f.write('set view 45, 345\n')
f.write('set size 1.0 ,1.0\n')
f.write('set autoscale\n')
#f.write('set palette gray\n')
#f.write('set palette defined ( 0 0.9 0.9 0.9, 1 0.9 0.9 0.9 )\n')
f.write('set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )\n')
f.write('set pm3d at s explicit\n')
#f.write('unset colorbox\n\n')

f.write('set line style 1 lt 1 lw 0.01 lc rgb "black"\n')
f.write('set line style 2 lt 1 lw 1.0 lc rgb "black"\n\n')

f.write('splot "meshDef_low.txt" using 1:2:3:4 w pm3d notitle, \\\n')
f.write('      "meshDef_low.txt" using 1:2:3:4 w l ls 1 notitle, \\\n')
f.write('      "meshDef_up.txt" using 1:2:3:4 w pm3d notitle, \\\n')
f.write('      "meshDef_up.txt" using 1:2:3:4 w l ls 1 notitle, \\\n')
f.write('      "wingDef_border.txt" using 1:2:3 w l ls 2 notitle\n\n')

f.write('print "Press return to export current image into deformAero3D.eps file...";\n')
f.write('pause -1\n\n')

f.write('set term post eps enh "Times-Italic" 10 color solid\n')
f.write('set output "../Fig/deformAero3D.eps"\n')
f.write('replot\n')

f.close()

print ' OK!'
print '----------------------------------------'

# -------------------------------------------------------
# *** PLOT ***
# -------------------------------------------------------
os.system('gnuplot sec3D.plt')
os.system('gnuplot def3D.plt')
os.system('rm *.txt')
os.system('rm sec3D.plt')
os.system('rm def3D.plt')
