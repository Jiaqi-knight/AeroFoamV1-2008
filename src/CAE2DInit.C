//==============================================================================
// *** CREATE CAE TOOLBOX DATA STRUCTURES (2D CASE) ***
//==============================================================================
Info << "========================================"           << nl;
Info << " Creating CAE2D data structures..."                 << nl;
Info << "========================================"           << nl;

// Initialize cellDisplacement and pointDisplacement fields
# include "DisplacementInit.C"

// Variables definition
label id_bodyPatch = 0, N = 0;                     // i, id_L previously defined in CreateConnectivity.C   
scalar pi = 3.14159265358979, gamma_ = 1.4;        // Constants
scalar DeltaZ = 1.0, meshCheck = 1.0;              // Extruded mesh checks                          
scalar Poo_ = 0.0, Too_ = 0.0, Uoo_ = 0.0;         // Reference primitive variables 
scalar qoo_ = 0.0, Moo_ = 0.0;                     // Derived reference variables (and domain width Dz)
word boundaryName, referenceFrame;                 // Boundary name and reference frame type
vector i_wind_, j_wind_, k_wind_;                  // Wind frame of reference axes
vector i_hat_, j_hat_, k_hat_;                     // Body frame of reference axes
scalar Cref, Lref;                                 // Reference lengths
vector Xref, Xhinge;                               // Reference posistion for C_MX computation and aileron deflection
scalar h, h_p, a, a_p, d, d_p;                     // Aeroelastic airfoil d.o.f. ( see Bisplinghoff pag. 533 + mods ) 
scalar C_L, C_D, C_MX, C_H;                        // Aerodynamic coefficients
FILE* fCAE2D;                                      // Statistics file pointer
//-------------------------------------------------//
scalar solutionType;                               // SolutionType (1 for forced motion, 2 for free motion ODEs with 2 d.o.f.)
label Mov_h, Mov_a, Mov_d;                         // Movement type
scalar A0_h, A1_h, f_h, tau_h;                     // Plunging movement parameters
scalar A0_a, A1_a, f_a, tau_a;                     // Pitching movement parameters
scalar A0_d, A1_d, f_d, tau_d;                     // Aileron movement parameters
scalar m, d_AECG, I_CG, K_hh, K_aa;                // 2 d.o.f. aero-elastic model parameters
scalar I_dd, C_dd, K_dd;                           // Aileron aero-elastic model parameters (for buzz)
//-------------------------------------------------//
labelField id_bodyFace;                            // Connectivity array from local to global face ids
labelField id_bodyCell;                            // Connectivity array of adjacent cells
vectorField xx;                                    // Coordinates of boundary faces centres
vectorField nn0, nn;                               // Boundary faces normal vectors ( undeformed and deformed )
vectorField uu, uu_p;                              // Boundary faces displacement and velocity
scalarField SS;                                    // Boundary faces areas
scalarField VVbn;                                  // Normal transpiration velocity
scalarField CCp;                                   // Pressure coefficient

//------------------------------------------------------------------------------
// Read CAE2D input parameters
//------------------------------------------------------------------------------

// Reference frame (body vs wind)
referenceFrame = word("body");
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("referenceFrame") )
referenceFrame = word( mesh.schemesDict().subDict("CAE2DToolbox").lookup("referenceFrame") );
i_wind_ = vector(1.0, 0.0, 0.0);
j_wind_ = vector(0.0, 1.0, 0.0);
k_wind_ = vector(0.0, 0.0, 1.0);
i_hat_  = vector(1.0, 0.0, 0.0);
j_hat_  = vector(0.0, 1.0, 0.0);
k_hat_  = vector(0.0, 0.0, 1.0);

// Cref (chord length)
Cref = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Cref") )
Cref = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Cref") );

// Lref (usually half chord length)
Lref = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Lref") )
Lref = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Lref") );

// Xref (usually at a quarter chord position)
Xref = vector(0.25, 0.0, 0.0);
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Xref.x()") )
Xref.x() = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Xref.x()") );
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Xref.y()") )
Xref.y() = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Xref.y()") );
Xref.z() = mesh.Cf()[0].z(); // to be verified with meshCheck

// Xhinge ( for aileron deflection )
Xhinge = vector(0.7, 0.0, 0.0);
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Xhinge.x()") )
Xhinge.x() = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Xhinge.x()") );
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Xhinge.y()") )
Xhinge.y() = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Xhinge.y()") );
Xhinge.z() = mesh.Cf()[0].z(); // to be verified with meshCheck

// SolutionType (1 for forced movement, 2 for free motion 2 d.o.f.)
solutionType = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("solutionType") )
solutionType = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("solutionType") );    
    
//----------------------------------    
// *** Forced motion parameters ***  
//----------------------------------
// Plunging movement ( translation h of the EA ) parameters
Mov_h = 0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Mov_h") )
Mov_h = readLabel( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Mov_h") ); 

A0_h = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("A0_h") )
A0_h = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("A0_h") ); 

A1_h = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("A1_h") )
A1_h = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("A1_h") ); 

f_h = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("f_h") )
f_h = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("f_h") ); 

tau_h = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("tau_h") )
tau_h = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("tau_h") );  
     
// Pitching movement ( rotation a around the EA ) parameters
Mov_a = 0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Mov_a") )
Mov_a = readLabel( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Mov_a") ); 

A0_a = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("A0_a") )
A0_a = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("A0_a") ); 

A1_a = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("A1_a") )
A1_a = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("A1_a") ); 

f_a = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("f_a") )
f_a = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("f_a") ); 

tau_a = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("tau_a") )
tau_a = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("tau_a") );       

// Aileron movement ( rotation d around hinge position ) parameters
Mov_d = 0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("Mov_d") )
Mov_d = readLabel( mesh.schemesDict().subDict("CAE2DToolbox").lookup("Mov_d") ); 

A0_d = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("A0_d") )
A0_d = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("A0_d") ); 

A1_d = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("A1_d") )
A1_d = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("A1_d") ); 

f_d = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("f_d") )
f_d = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("f_d") ); 

tau_d = 0.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("tau_d") )
tau_d = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("tau_d") );         
    
//---------------------------------   
// *** Free motion parameters ***  
//---------------------------------   

// Model mass
m = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("m") )
m = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("m") ); 

// Distance between the elastic axis (EA) and the center of gravity (CG) 
// WARNING: only in input EA and not AE
d_AECG = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("d_EACG") )
d_AECG = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("d_EACG") ); 
     
// Inertia moment around the center of gravity (CG)
I_CG = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("I_CG") )
I_CG = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("I_CG") ); 
     
// Translational stiffness
K_hh = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("K_hh") )
K_hh = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("K_hh") ); 

// Rotational stiffness
K_aa = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("K_aa") )
K_aa = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("K_aa") ); 

// Aileron moment of inertia around hinge axis
I_dd = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("I_dd") )
I_dd = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("I_dd") ); 

// Aileron rotational damping coefficient
C_dd = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("C_dd") )
C_dd = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("C_dd") ); 

// Aileron rotational stiffness
K_dd = 1.0;
if ( mesh.schemesDict().subDict("CAE2DToolbox").found("K_dd") )
K_dd = readScalar( mesh.schemesDict().subDict("CAE2DToolbox").lookup("K_dd") ); 

//------------------------------------------------------------------------------
// Initialize reference variables and useful arrays
//------------------------------------------------------------------------------
// Compute aerodynamic panels length along Z (in OpenFoam 2D triangular meshes are extruded along Z axis)
DeltaZ    = max( mesh.allPoints().component(2) ) - min( mesh.allPoints().component(2) );
meshCheck = mag( DeltaZ - 2*mesh.Cf()[0].z() );

// Loop on boundary patches
forAll( mesh.boundaryMesh(), iPatch )
{
    // *** Initialize reference variables ***
    if( ( mesh.boundaryMesh().physicalTypes()[iPatch] == "supersonicInlet" ) || 
        ( mesh.boundaryMesh().physicalTypes()[iPatch] == "Riemann" ) )
    {
        id_L    = mesh.boundaryMesh()[iPatch].faceCells()[0];   
        gamma_  = gamma.value(); 
        Poo_    = p[id_L];
        Too_    = T[id_L];
        Uoo_    = mag( U[id_L] );
        Moo_    = Uoo_/Foam::sqrt( gamma_*R.value()*Too_ );
        qoo_    = 0.5*gamma_*Poo_*Moo_*Moo_; 
        i_wind_ = U[id_L]/mag( U[id_L] );
        k_wind_ = vector(0.0, 0.0, 1.0);
        j_wind_ = ( k_wind_ ^ i_wind_ );         
    }
    
    // *** Initialize useful arrays ***
    if( ( mesh.boundaryMesh().names()[iPatch] == "body" ) || ( mesh.boundaryMesh().names()[iPatch] == "wing" ) )
    {
        // Memory allocation
        id_bodyPatch = iPatch;
        boundaryName = word( mesh.boundaryMesh().names()[id_bodyPatch] );  
        N            = mesh.boundaryMesh()[iPatch].faceAreas().size();   
        id_bodyFace  = labelField( N, 0 ); // Connectivity array from local to global face ids
        id_bodyCell  = labelField( N, 0 ); // Connectivity array of adjacent cells
        xx   = vectorField( N );           // Coordinates of boundary faces centres
        nn0  = vectorField( N );           // Boundary faces normal vectors ( undeformed )
        nn   = vectorField( N );           // Boundary faces normal vectors ( deformed )
        uu   = vectorField( N );           // Boundary faces displacement
        uu_p = vectorField( N );           // Boundary faces velocity
        SS   = scalarField( N );           // Boundary faces areas
        VVbn = scalarField( N );           // Normal transpiration velocity
        CCp  = scalarField( N );           // Pressure coefficient
        
        // Loop on boundary faces
        forAll( mesh.boundaryMesh()[iPatch].faceAreas(), ii )
        {
            // Local to Global index and arrays initialization  
            id_L = mesh.boundaryMesh()[iPatch].faceCells()[ii];       
            i    = mesh.boundaryMesh()[iPatch].whichFace(ii);
            id_bodyFace[ii] = i;
            id_bodyCell[ii] = id_L;
            xx[ii]   = mesh.Cf()[i];
            nn0[ii]  = mesh.boundaryMesh()[iPatch].faceAreas()[ii];
            // REMARK: nn0, nn vectors point outside the body
            nn0[ii]  = -nn0[ii]/mag( nn0[ii] );
            nn[ii]   = nn0[ii];
            uu[ii]   = vector(0.0, 0.0, 0.0);
            uu_p[ii] = vector(0.0, 0.0, 0.0);
            SS[ii]   = mag( mesh.boundaryMesh()[iPatch].faceAreas()[ii] )/DeltaZ;
            VVbn[ii] = 0.0;
            CCp[ii]  = 0.0;   
        }            
    }
}      

// Choose frame of reference
if ( referenceFrame == "wind" )
{  
    i_hat_ = i_wind_;
    j_hat_ = j_wind_;
    k_hat_ = k_wind_;     
}
else // if ( referenceFrame == "body" ) DEFAULT
{
    i_hat_ = vector(1.0, 0.0, 0.0);
    j_hat_ = vector(0.0, 1.0, 0.0);
    k_hat_ = vector(0.0, 0.0, 1.0);
}

// Airfoil state initialization
h   = 0.0;
a   = 0.0;
h_p = 0.0;
a_p = 0.0;
d   = 0.0;
d_p = 0.0;

// Aerodynamic coefficients initialization
C_L  = 0.0;
C_D  = 0.0;
C_MX = 0.0;
C_H  = 0.0;

//------------------------------------------------------------------------------
// Statistics
//------------------------------------------------------------------------------

// Statistics file initialization
if ( runTime.value() == 0.0 )
{
    fCAE2D = fopen("./Log/AerodynamicLoads.txt", "w");
    fclose(fCAE2D);
    fCAE2D = fopen("./Log/Displacements.txt", "w");
    fclose(fCAE2D);
}

// Print statistics and checks
Info << " gamma    [-] = " << gamma_                            << nl;
Info << " Poo     [Pa] = " << Poo_                              << nl;
Info << " Too      [K] = " << Too_                              << nl;
Info << " Uoo    [m/s] = " << Uoo_                              << nl;
Info << " Moo      [-] = " << Moo_                              << nl;
Info << " qoo     [Pa] = " << qoo_                              << nl;
Info << " Cref     [m] = " << Cref                              << nl;
Info << " Lref     [m] = " << Lref                              << nl;
Info << " Xref     [m] = " << Xref                              << nl;
Info << " Xhinge   [m] = " << Xhinge                            << nl;
Info << " i_hat    [-] = " << i_hat_                            << nl;
Info << " j_hat    [-] = " << j_hat_                            << nl;
Info << " k_hat    [-] = " << k_hat_                            << nl;
Info << " refFrame     = " << referenceFrame                    << nl;
if (meshCheck < 1e-6 )  Info << " ExtrudedMesh = OK! "          << nl; 
if (meshCheck >= 1e-6 ) Info << " ExtrudedMesh = WARNING! "     << nl;  
Info << " # of BndFaces @ " << boundaryName << " = " << N       << nl;
Info << "----------------------------------------"        << nl << nl;  
