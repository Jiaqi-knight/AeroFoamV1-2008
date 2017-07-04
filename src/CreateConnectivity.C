//==============================================================================
// *** CREATE CONNECTIVITY ***
//==============================================================================
Info << "========================================"           << nl;
Info << " Creating Mesh Connectivity..."                     << nl;
Info << "========================================"           << nl;

// Variables definition
label i, id_L, id_R, id_LL, id_RR, id_adv;
vector x_LL, x_RR;
scalar Sf = 0.0, dxInternal = 0.0, dxBoundary = 0.0;

// Adaptive correction of dx until id_LL/id_RR changes to increase routine 
// robustness for skewed elements ( works fine only for hexa cells )
scalar adaptiveConnectivitySearch = 0.0;
if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("adaptiveConnectivity") )
adaptiveConnectivitySearch = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("adaptiveConnectivity") ); 
if ( adaptiveConnectivitySearch == 1 )
{
    Info << "Using adaptive algorithm..." << nl;
}
else
{
    Info << "Using standard algorithm..." << nl;
}

// Create id_RR = (j + 2) and id_LL = (j - 1) lists ( see Jameson )
labelField id_LL_i( mesh.faces().size(), 0 );
labelField id_RR_i( mesh.faces().size(), 0 ); 

//------------------------------------------------------------------------------
// *** Loop on internal faces ***
//------------------------------------------------------------------------------
Info << "Loop on Internal Faces..." << nl;
id_adv = mesh.Sf().size()/10;
forAll( mesh.Sf(), i )
{

    // Input 
    id_L = mesh.faceOwner()[i];
    id_R = mesh.faceNeighbour()[i];     
    
    // *** ADAPTIVE VERSION *** ( works better with quad/hexa cells )
    if ( adaptiveConnectivitySearch == 1 )
    {
        // L) Incremental search ( until id_L != id_LL ) should increase code robustness
        dxInternal = 1.0;
        id_LL      = id_L;
        while( ( id_LL == id_L ) && ( dxInternal < 5.0 ) )
        {
            x_LL  = mesh.Cf()[i] - ( mesh.C()[id_R] - mesh.C()[id_L] )*dxInternal;
            id_LL = meshSearch(mesh).findNearestCell( x_LL, id_L );
            if ( id_LL < 0 ) id_LL = id_L;   
            dxInternal = dxInternal + 0.5;    
        }
        
        // R) Incremental search ( until id_R != id_RR ) should increase code robustness
        dxInternal = 0.25;
        id_RR      = id_R;
        while( ( id_RR == id_R ) && ( dxInternal < 5.0 ) )
        {
            x_RR  = mesh.Cf()[i] + ( mesh.C()[id_R] - mesh.C()[id_L] )*dxInternal;
            id_RR = meshSearch(mesh).findNearestCell( x_RR, id_R );
            if ( id_RR < 0 ) id_RR = id_R;   
            dxInternal = dxInternal + 0.25; 
        }
    }
    // *** NON ADAPTIVE VERSION *** ( works better for tri/tetra cells )
    else
    {  
        // Find id_RR = (j + 2) and id_LL = (j - 1) ( see Jameson )
        dxInternal = 1.5; 
        x_LL  = mesh.C()[id_L] - ( mesh.C()[id_R] - mesh.C()[id_L] )*dxInternal;
        x_RR  = mesh.C()[id_R] + ( mesh.C()[id_R] - mesh.C()[id_L] )*dxInternal;
        id_LL = meshSearch(mesh).findNearestCell( x_LL, id_L );
        id_RR = meshSearch(mesh).findNearestCell( x_RR, id_R );
        if ( id_LL < 0 ) id_LL = id_L;
        if ( id_RR < 0 ) id_RR = id_R;
    } 
        
    // Update connectivity array
    id_LL_i[i] = id_LL;
    id_RR_i[i] = id_RR;
        
    // Wait bar
    if ( i >= id_adv )
    {
        printf("%3.0f %% \n", i*100.0/mesh.Sf().size() );
        id_adv = id_adv + mesh.Sf().size()/10;
    }   
     
}
Info << "----------------------------------------" << nl;

//------------------------------------------------------------------------------
// *** Loop on boundary faces ***
//------------------------------------------------------------------------------
Info << "Loop on Boundary Faces..." << nl;
id_adv = ( mesh.faces().size() - mesh.Sf().size() - 1 )/10;
forAll( mesh.boundaryMesh(), iPatch )
{
    forAll( mesh.boundaryMesh()[iPatch].faceAreas(), ii )
    {
        
        // Input
        id_L = mesh.boundaryMesh()[iPatch].faceCells()[ii];           
        i    = mesh.boundaryMesh()[iPatch].start()+ii;
                 
        // *** ADAPTIVE VERSION *** ( works better with quad/hexa cells )
        if ( adaptiveConnectivitySearch == 1 )
        {
            // L) Incremental search ( until id_L != id_LL ) should increase code robustness
            dxBoundary = 0.5;
            id_LL      = id_L;
            Sf = mag( mesh.boundaryMesh()[iPatch].faceAreas()[ii] );
            while( ( id_LL == id_L ) && ( dxBoundary < 5.0 ) )
            {
                x_LL  = mesh.C()[id_L] - ( mesh.boundaryMesh()[iPatch].faceAreas()[ii]/Sf*mesh.V()[id_L]/Sf )*dxBoundary;
                id_LL = meshSearch(mesh).findNearestCell( x_LL, id_L );
                if ( id_LL < 0 ) id_LL = id_L;   
                dxBoundary = dxBoundary + 0.5;    
            }
        }     
        // *** NON ADAPTIVE VERSION *** ( works better for tri/tetra cells )
        else
        {      
            // Find id_LL = (j - 1) ( see Jameson ) 
            dxBoundary = 3.0;
            Sf = mag( mesh.boundaryMesh()[iPatch].faceAreas()[ii] );
            x_LL  = mesh.C()[id_L] - ( mesh.boundaryMesh()[iPatch].faceAreas()[ii]/Sf*mesh.V()[id_L]/Sf )*dxBoundary;
            id_LL = meshSearch(mesh).findNearestCell( x_LL, id_L );
            if ( id_LL < 0 ) id_LL = id_L;
        }    
        
        // Update connectivity array
        id_LL_i[i] = id_LL;
        id_RR_i[i] = -1;
                
        // Wait bar
        if ( ( i - mesh.Sf().size() ) >= id_adv )
        {
            printf("%3.0f %% \n", ( i - mesh.Sf().size() - 1 )*100.0/( mesh.faces().size() - mesh.Sf().size() - 1 ) );
            id_adv = id_adv + ( mesh.faces().size() - mesh.Sf().size() )/10;
        }   
        
    }
}      
Info << "----------------------------------------"           << nl << nl;  

//------------------------------------------------------------------------------
// *** Load connectivity matrix corrected with MDS ***
//------------------------------------------------------------------------------
// REMARK: This should increase accuracy in the domain and on the boundary
scalar load_conn = 0.0;
if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("loadConnectivity") )
load_conn = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("loadConnectivity") ); 
if ( ( load_conn == 1 ) || ( load_conn == 2 ) )
{ 

    Info << "========================================" << nl;
    Info << " Loading Mesh Connectivity..."            << nl;
    Info << "========================================" << nl;
    // Open file anc check for errors
    FILE *ff;
    ff = fopen("./Data/LR2LLRR.txt", "r");
    if ( ff == NULL ) Info << "ERROR: File LR2LLRR not found!" << nl;
    
    // Loop on mesh faces
    forAll( id_LL_i, i )
    {
        fscanf(ff, "%i %i %i %i\n", &id_L, &id_LL, &id_R, &id_RR );
        id_LL_i[i] = id_LL-1;
        id_RR_i[i] = id_RR-1;
    } 
    
    // Close file
    fclose(ff);
    Info << " Done!"                                   << nl; 
    Info << "----------------------------------------" << nl << nl; 
}

//------------------------------------------------------------------------------
// *** Write connectivity matrix to be corrected with MDS ***
//------------------------------------------------------------------------------
if ( ( load_conn == -1 ) || ( load_conn == 2 ) )
{
    // Write connectivity structures id_L->id_LL, id_R->id_RR on file
    Info << "========================================" << nl;
    Info << " Writing Mesh Connectivity..."            << nl;
    Info << "========================================" << nl;
    FILE *ff;
    ff = fopen("./Log/LR2LLRR.txt", "w");
    forAll( id_LL_i, i )
    {
        // Input 
        id_L  = mesh.faceOwner()[i];
        id_R  = mesh.faceNeighbour()[i]; 
        id_LL = id_LL_i[i];
        id_RR = id_RR_i[i];
        fprintf( ff, "%i %i %i %i\n", id_L+1, id_LL+1, id_R+1, id_RR+1 );
        
    }
    fclose(ff);
    Info << " Done!"                                   << nl;
    Info << "----------------------------------------" << nl << nl; 
}              

//==============================================================================
// *** TEST CONNECTIVITY ***
//==============================================================================

// Read from control file
// BUG ??? Some frames are repeated 
scalar test_conn = 0.0;
if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("testMesh") )
test_conn = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("testMesh") ); 

// Test
if ( test_conn > 0 )
{
    Info << "========================================" << nl;
    Info << " Testing Mesh Connectivity..."            << nl;
    Info << "========================================" << nl;
    // Loop on internal and boundary Faces
    Info << "i"  << " " << "id_LL" << " " << "id_L" << " " << "id_R" << " " << "id_RR" << nl;
    forAll( id_LL_i, i )
    {
        // Input 
        id_L  = mesh.faceOwner()[i];
        id_R  = mesh.faceNeighbour()[i]; 
        id_LL = id_LL_i[i];
        id_RR = id_RR_i[i];
        
        // Statistics
        Info << i << ") " << id_LL << " " << id_L << " " << id_R << " " << id_RR << nl;
        
        // Plot
        T = 0*T; 
        T[id_L]  = 100;
        if ( id_R >= 0 ) T[id_R] = 200;
        T[id_LL] = 300;
        if ( id_RR >= 0 ) T[id_RR] = 400;
                
        if ( ( id_RR_i[i] == -1 ) || ( test_conn == 1.0 ) ) // plot only boundary elements when testMesh = 2
        {
        // Update statistics
        scalar dTime    = readScalar( runTime.controlDict().lookup("writeInterval") );
        runTime.value() = runTime.value() + dTime - runTime.deltaT().value();
        runTime++;
        runTime.write();
        }
   
    }
    Info << "----------------------------------------" << nl << nl;  
    
    // Output    
    return(0);  
}
