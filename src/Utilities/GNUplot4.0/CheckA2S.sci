//==========================================================================
// *** CHECK CONNECTIVITY MATRIX A2S ***
//==========================================================================
clear; clc; xdel(winsid())

// REAMARK: works only if structural mesh is a plane // Z-Axis

// Load Structural mesh
fid = mopen("StructuralMesh.vertices");
Nv = mfscanf(fid, "%i\n");
for k = 1:Nv
   vertices(k,:) = mfscanf(fid,"%lf %lf %lf\n" );
end
mclose(fid);
fid = mopen("StructuralMesh.elements");
Ne = mfscanf(fid, "%i\n");
for k = 1:Ne
   elements(k,:) = mfscanf(fid,"%i %i %i\n" );
end
mclose(fid);
vertices = vertices';
elements = elements';

// Plot Structural Mesh
scf(1);
for i = 1:length(elements(1, :))
   id_ele = [ elements(:, i) ]; 
   id_ele = [ id_ele; id_ele(1) ];
   plot(vertices(1, id_ele), vertices(2, id_ele), 'k-' ); 
end

// Color array (200x)
colorv = [ 'r', 'g', 'b', 'c' ];
colorv = [ colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv, ...
          colorv, colorv, colorv, colorv, colorv ];
          
// Load Aerodynamic CGs and connectivity
fid = mopen("CG_a.txt");
k = 1;
while ( ~meof(fid) )
   CGa(k,:) = mfscanf(fid,"%lf %lf %lf\n" );
   k = k + 1;
end
mclose(fid);
fid = mopen("id_a2s.txt");
k = 1;
while ( ~meof(fid) )
   a2s(k,:) = mfscanf(fid,"%i %i\n" );
   k = k + 1;
end
mclose(fid);
a2s  = a2s(:,2)' + 1;
CGa = CGa';

// Plot CGs
for i = 1:Ne
   CCGai = [];
   for j = 1: length(a2s)
      if ( a2s(j) == i ) 
         CCGai = [ CCGai, CGa(:,j)];
      end
   end
   Interface(i).entries = CCGai;
end
for i = 1: Ne
   plot(Interface(i).entries(1,:), Interface(i).entries(2,:), colorv(i) + '.', 'markersize', 2) 
end

// Save figure
xs2eps(1, "../Fig/CheckA2S.eps", 1)
