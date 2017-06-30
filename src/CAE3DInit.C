//==============================================================================
// *** CREATE CAE TOOLBOX DATA STRUCTURES (3D WINGS) ***
//==============================================================================
Info << "========================================"             << nl;
Info << " Creating CAE3D data structures..."                   << nl;
Info << "========================================"             << nl;

// Initialize cellDisplacement and pointDisplacement fields
# include "DisplacementInit.C"

// Variables definition
//-------------------------------------------// [S]
label is = 0, W = 0;                         // i, id_L previously defined in CreateConnectivity.C
label Nv_s, Ne_s, Nmodes, NactiveMode;       // Structural model dimensions
scalar restart, solutionType, t0Forcing_s;   // Structural solution parameters
scalar qMax_s, tauMax_s, kMax_s, epsU_s;     // Modal forced movement parameters
scalar tPrint, dtPrint, err_qq_s, err_QQ_s;  // Convergence on generalized displacements and forces           
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
labelField id_A_s, id_B_s, id_C_s;           // Structural mesh elements
vectorField xx_s_v;                          // Structural mesh vertices
vectorField nn0_s_e, nn_s_e;                 // Structural mesh elements normal vectors ( undeformed and deformed )
vectorField uu_s_v, uup_s_v;                 // Structural vertices displacement and velocity
vectorField RR_s_v;                          // Sum of aerodynamic forces distributed on structural vertices
scalarField ff_s, mm_s, cc_s;                // Structural experimental modal data
scalarField qq_s, qqp_s, qqo_s;              // Generalized structural displacements
scalarField QQ0_s, QQ_s, QQo_s;              // Generalized aerodynamic forces

//-------------------------------------------// [A]
label ia = 0, N = 0, id_bodyPatch = 0;       // Dimensions
scalar Poo_ = 0.0, Too_ = 0.0, Uoo_ = 0.0;   // Reference aerodynamic variables
scalar qoo_ = 0.0, Moo_ = 0.0, gamma_ = 0.0; // Derived reference aerodynamic variables 
vector i_wind_, j_wind_, k_wind_;            // Wind frame of reference
scalar Sref_a, Lref_a;                       // Aerodynamic reference surface and length
vector Xref_a;                               // Reference point for the moment
scalar C_L, C_D;                             // Aerodynamic lift and drag coefficients (in wind axes)
scalar C_FX, C_FY, C_FZ, C_MX, C_MY, C_MZ;   // Aerodynamic force coefficients
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
labelField id_bodyFace;                      // Connectivity array from local to global face ids
labelField id_bodyCell;                      // Connectivity array of adjacent cells
labelField id_a2s;                           // A2S connectivity matrix
vectorField xx_a;                            // Coordinates of boundary faces centres
vectorField nn0_a, nn_a;                     // Boundary faces normal vectors ( undeformed and deformed )
scalarField SS_a;                            // Boundary faces areas
scalarField VVbn_a;                          // Normal transpiration velocity
scalarField CCp_a;                           // Pressure coefficient
vectorField RR_a;                            // Aerodynamic force evaluated on aerodynamic elements
vectorField uu_a, uup_a;                     // Body velocity projected on aerodynamic nodes

/******************************************************************************/
/* INPUT PARAMTERS                                                            */
/******************************************************************************/

 // (AGARD 445.6 wing aerodynamic data obtained by ESagard.sci)
 Sref_a = 0.3513;                        // Half wing surface
 Lref_a = 0.4676/2.0;                    // Half mean aerodynamic chord (m.a.c.)
 Xref_a = vector( 0.4938, 0.3547, 0.0 ); // 1/4c Mean aerodynamic chord (X, Y, Z)

 // Read aerodynamic reference surface Sref_a
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("Sref_a") )
 Sref_a = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("Sref_a") );
 
 // Read aerodynamic reference length Lref_a
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("Lref_a") )
 Lref_a = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("Lref_a") );
 
 // Read aerodynamic reference point Xref_a
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("Xref_a.x()") )
 Xref_a.x() = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("Xref_a.x()") );
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("Xref_a.y()") )
 Xref_a.y() = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("Xref_a.y()") );
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("Xref_a.z()") )
 Xref_a.z() = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("Xref_a.z()") );

/******************************************************************************/

 // (AGARD 445.6 wing structural data obtained by ESagard.sci)
 restart      = 0;        // restart structural computation and load modal state from file
 solutionType = 1;        // structural solution type: (0) for static modal solution (1) forced modal motion
 Nmodes       = 4;        // # of active modes                                     
 NactiveMode  = 1;        // forced mode
 kMax_s       = 2.0;      // maximum reduced frequency of interest
 epsU_s       = 0.017455; // tan(1 [deg]) TODO: Linearity check
 t0Forcing_s  = 0.0;      // dimensional time for the modal forced motion to start
 dtPrint      = 1e-5;     // dimensional time to save CAE3D statistics

 // Read restart
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("restart") )
 restart = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("restart") );
 
 // Read modal forcing start time t0Forcing_s
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("forcingStartTime") )
 t0Forcing_s = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("forcingStartTime") );
 
  // Read modal forcing start time t0Forcing_s
 tPrint = runTime.value();
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("printInterval") )
 dtPrint = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("printInterval") );
 
 // Read solutionType
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("solutionType") )
 solutionType = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("solutionType") );

 // Read Nmodes
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("Nmodes") )
 Nmodes = readLabel( mesh.schemesDict().subDict("CAE3DToolbox").lookup("Nmodes") );
 
 // Read NactiveMode
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("NactiveMode") )
 NactiveMode = readLabel( mesh.schemesDict().subDict("CAE3DToolbox").lookup("NactiveMode") );
 
 // Read maximum reduced frequency kMax_s
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("kMax") )
 kMax_s = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("kMax") );
 
 // Read maximum velocity ratio epsU_s
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("epsU") )
 epsU_s = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("epsU") );
 
 // Read maximum modal displacement qMax_s ( optional, otherwise it is computed )
 qMax_s = 0.0;
 if ( mesh.schemesDict().subDict("CAE3DToolbox").found("qMax") )
 qMax_s = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("qMax") );
 
/******************************************************************************/

//------------------------------------------------------------------------------
// Read structural mesh and experimental modal data
//------------------------------------------------------------------------------
// Read structural mesh data
CAE3DreadStructralMesh( &xx_s_v, &id_A_s, &id_B_s, &id_C_s );
Nv_s = xx_s_v.size();
Ne_s = id_A_s.size();

// Read structural experimental modal data
scalar UUMax_s = 0.0;
Matrix<scalar> UUx_s(Nv_s, Nmodes);
Matrix<scalar> UUy_s(Nv_s, Nmodes);  
Matrix<scalar> UUz_s(Nv_s, Nmodes);
CAE3DreadStructuralModes( Nv_s, Nmodes, &UUx_s, &UUy_s, &UUz_s, &ff_s, &mm_s, &cc_s );
qq_s  = scalarField(Nmodes, 0.0);
qqp_s = scalarField(Nmodes, 0.0);
qqo_s = scalarField(Nmodes, 0.0);
if ( restart == 1 ) CAE3DreadModalState(&qq_s, &qqp_s);
err_qq_s = 0.0;

// Allocate memory for generalized forces arrays
QQ0_s = scalarField(Nmodes, 0.0);
QQ_s  = scalarField(Nmodes, 0.0);
QQo_s = scalarField(Nmodes, 0.0);
err_QQ_s = 0.0;

// Allocate memory for load and movement arrays
RR_s_v  = vectorField(Nv_s, vector(0.0, 0.0, 0.0) );
nn0_s_e = vectorField(Ne_s, vector(0.0, 0.0, 0.0) );
nn_s_e  = vectorField(Ne_s, vector(0.0, 0.0, 0.0) );
uu_s_v  = vectorField(Nv_s, vector(0.0, 0.0, 0.0) );
uup_s_v = vectorField(Nv_s, vector(0.0, 0.0, 0.0) );

// Initialize undeformed structural mesh normal vector
CAE3DcomputeNn_s_e( &xx_s_v, &uu_s_v, &id_A_s, &id_B_s, &id_C_s, &nn0_s_e );

//------------------------------------------------------------------------------
// Initialize reference variables and useful arrays
//------------------------------------------------------------------------------
// Loop on boundary patches
forAll( mesh.boundaryMesh(), iPatch )
{
    // *** Initialize reference variables ***
    if( ( mesh.boundaryMesh().physicalTypes()[iPatch] == "supersonicInlet" ) || 
        ( mesh.boundaryMesh().physicalTypes()[iPatch] == "Riemann"         ) )
    {
        id_L   = mesh.boundaryMesh()[iPatch].faceCells()[0];    

        Poo_   = p[id_L];
        Too_   = T[id_L];
        Uoo_   = mag( U[id_L] );

        gamma_ = gamma.value();
        Moo_   = Uoo_/Foam::sqrt( gamma_*R.value()*Too_ );
        qoo_   = 0.5*gamma_*Poo_*Moo_*Moo_;

        i_wind_ = U[id_L]/mag( U[id_L] );
        j_wind_ = vector( 0.0, 1.0, 0.0 );
        k_wind_ = ( i_wind_ ^ j_wind_ );         
        k_wind_ = k_wind_/mag(k_wind_);
    }

    // *** Initialize connectivity Aerodynamics2Structure (id_a2s) array ***
    // Each aerodynamic element id contains the corrisponding structural subdomain id
    if( ( mesh.boundaryMesh().names()[iPatch] == "body" ) || 
        ( mesh.boundaryMesh().names()[iPatch] == "wing" ) )
    {
        // Memory allocation
        id_bodyPatch = iPatch;
        N = mesh.boundaryMesh()[iPatch].faceAreas().size();      
        id_a2s       = labelField( N, -1 ); // Connectivity array from aerodynamic to structural elements   
        id_bodyFace  = labelField( N, 0 );  // Connectivity array from local to global face ids
        id_bodyCell  = labelField( N, 0 );  // Connectivity array of adjacent cells
        xx_a   = vectorField( N );          // Coordinates of boundary faces centres
        nn0_a  = vectorField( N );          // Boundary faces normal vectors ( undeformed )
        nn_a   = vectorField( N );          // Boundary faces normal vectors ( deformed )
        SS_a   = scalarField( N );          // Boundary faces areas
        VVbn_a = scalarField( N );          // Normal transpiration velocity
        CCp_a  = scalarField( N );          // Pressure coefficient 
        RR_a   = vectorField( N );          // Aerodynamic force
        uu_a   = vectorField( N );          // Body displacement projected on aerodynamic nodes
        uup_a  = vectorField( N );          // Body velocity projected on aerodynamic nodes
         
        // Loop on boundary faces        
        forAll( mesh.boundaryMesh()[iPatch].faceAreas(), ii )
        {
            // Local to Global index and arrays initialization  
            id_L = mesh.boundaryMesh()[iPatch].faceCells()[ii];      
            ia   = mesh.boundaryMesh()[iPatch].whichFace(ii);
            is   = CAE3DbuildA2S( &mesh, ia, &xx_s_v, &id_A_s, &id_B_s, &id_C_s );
            id_a2s[ii]      = is;
            id_bodyFace[ii] = ia;
            id_bodyCell[ii] = id_L;
            xx_a[ii]   = mesh.Cf()[ia];
            nn0_a[ii]  = mesh.boundaryMesh()[iPatch].faceAreas()[ii];
            // REMARK: nn0, nn vectors are tilted and point outside the body
            nn0_a[ii]  = -nn0_a[ii]/mag( nn0_a[ii] );
            nn_a[ii]   = nn0_a[ii];
            SS_a[ii]   = mag( mesh.boundaryMesh()[iPatch].faceAreas()[ii] );
            VVbn_a[ii] = 0.0;
            CCp_a[ii]  = 0.0; 
            RR_a[ii]   = vector(0.0, 0.0, 0.0);  
            uu_a[ii]   = vector(0.0, 0.0, 0.0);  
            uup_a[ii]  = vector(0.0, 0.0, 0.0);
        }      
    }
} 

//------------------------------------------------------------------------------
// Aerodynamic to Structural interface (piecewise linear)
//------------------------------------------------------------------------------

// Build interface matrix (Tas^T is computed in place when needed)
Matrix<scalar> Tas(N, Nv_s);               
CAE3DbuildTas( &id_A_s, &id_B_s, &id_C_s, &id_a2s, &xx_s_v, &xx_a, &Tas );  

//------------------------------------------------------------------------------
// Compute forced modal motion parameters
//------------------------------------------------------------------------------
// TODO: *** Linearity check with epsU_s = tan(2 [deg]) ***
if ( qMax_s == 0 ) CAE3DforcedModalMotionParam( &UUx_s, &UUy_s, &UUz_s, NactiveMode, kMax_s, Uoo_, Lref_a, epsU_s, &UUMax_s, &tauMax_s, &qMax_s );
 
// Print NAEMO inf file (see Cavagna et al.)
CAE3DprintNAEMOInf( NactiveMode, qMax_s, Uoo_, Poo_/(R.value()*Too_), Lref_a, Sref_a, dtPrint ); 
 
//------------------------------------------------------------------------------
// Print statistics
//------------------------------------------------------------------------------
Info << "----------------------------------------"             << nl;
Info << " Aerodynamic reference variables"                     << nl;
Info << "----------------------------------------"             << nl;
Info << " Sref_a [m^2] = " << Sref_a                           << nl;
Info << " Lref_a   [m] = " << Lref_a                           << nl;
Info << " Xref_a   [m] = " << Xref_a                           << nl;
Info << " gamma    [-] = " << gamma_                           << nl;
Info << " Poo     [Pa] = " << Poo_                             << nl;
Info << " Too      [K] = " << Too_                             << nl;
Info << " Uoo    [m/s] = " << Uoo_                             << nl;
Info << " Moo      [-] = " << Moo_                             << nl;
Info << " qoo     [Pa] = " << qoo_                             << nl;
Info << " i_wind   [-] = " << i_wind_                          << nl;
Info << " j_wind   [-] = " << j_wind_                          << nl;
Info << " k_wind   [-] = " << k_wind_                          << nl;
Info << "----------------------------------------"             << nl;
Info << " Structural reference variables"                      << nl;
Info << "----------------------------------------"             << nl;
Info << " SolutionType = " << solutionType                     << nl;
Info << " NactiveMode  = " << NactiveMode                      << nl;
Info << " # of Modes   = " << Nmodes                           << nl;
Info << " UUMax    [m] = " << UUMax_s                          << nl;
Info << " kMax     [-] = " << kMax_s                           << nl;
Info << " u/Uoo    [-] = " << epsU_s                           << nl;
Info << " qMax     [-] = " << qMax_s                           << nl;
Info << " tauMax   [-] = " << tauMax_s                         << nl;
Info << " tMax     [s] = " << tauMax_s*Lref_a/Uoo_             << nl;
Info << " tStart   [s] = " << t0Forcing_s                      << nl;
Info << " dtPrint  [s] = " << dtPrint                          << nl;
Info << "----------------------------------------"             << nl;
Info << " Interface [ A ] <=> [ S ]"                           << nl;
Info << "----------------------------------------"             << nl;
Info << " # of [ A ] elements = " << id_a2s.size()             << nl;
Info << "      |   ^           "                               << nl;
Info << "      v   |           "                               << nl;
Info << " # of [ S ] elements = " << id_A_s.size()             << nl;
forAll(id_a2s, ii)
{
    if( id_a2s[ii] == -1 ) W = 1;
}
if ( W == 1 )
{
    Info << " Interface = Warning!"                            << nl;
}
else
{
    Info << " Interface = OK!"                                 << nl;
}
Info << "----------------------------------------"             << nl << nl;  

// Statistics file initialization
if ( runTime.value() == 0.0 )
{
    FILE* fCAE3D; 
    fCAE3D = fopen("./Log/AerodynamicLoads.txt", "w");
    fclose(fCAE3D);
    fCAE3D = fopen("./Log/qqs.txt", "w");
    fclose(fCAE3D);
    fCAE3D = fopen("./Log/QQa.txt", "w");
    fclose(fCAE3D); 
}

//------------------------------------------------------------------------------
// Test connectivity
//------------------------------------------------------------------------------
scalar testInterface = 0.0;
if ( mesh.schemesDict().subDict("CAE3DToolbox").found("testInterface") )
testInterface = readScalar( mesh.schemesDict().subDict("CAE3DToolbox").lookup("testInterface") ); 

if ( testInterface == 1 )
{
    // Open files
    FILE *fa, *fb;
    fa = fopen( "./Log/CG_a.txt", "w" );  
    fb = fopen( "./Log/id_a2s.txt", "w" );

    // Loop on boundary patches
    forAll( mesh.boundaryMesh(), iPatch )
    {
        // Check for boundary conditions
        if( ( mesh.boundaryMesh().names()[iPatch] == "body" ) || 
            ( mesh.boundaryMesh().names()[iPatch] == "wing" ) )
        {
            // Loop on boundary faces
            T = 0.0*T;
            forAll( mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                // Local to Global index and arrays initialization  
                id_L = mesh.boundaryMesh()[iPatch].faceCells()[ii];      
                ia   = mesh.boundaryMesh()[iPatch].whichFace(ii);
                is   = id_a2s[ii];

                // Print on file
                fprintf( fa, "%lf %lf %lf\n", mesh.Cf()[ia].x(), mesh.Cf()[ia].y(), mesh.Cf()[ia].z() );
                fprintf( fb, "%i %i\n", ia, is);

                // Modify T field
                T[id_L] = is*Foam::pow(-1, is);
            }    
        }
    }

    // Close files
    fclose(fa);
    fclose(fb);  

    // Create folder with T field modified to check interface
    scalar dTime    = readScalar( runTime.controlDict().lookup("writeInterval") );
    runTime.value() = runTime.value() + dTime - runTime.deltaT().value();
    runTime++;
    runTime.write();
    return(0);
}
