/*----------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
--------------------------------------------------------------------------------

Application:
    AeroFoam

Description:
    3D Euler equations cell-centered solver for aerodynamic (and aeroelastic) applications

Authors:    
    Giulio Romanelli, Elisa Serioli, (C) 2007-2008
    giulio.romanelli@gmail.com, elisa.serioli@gmail.com

Implemented Monotone Flux-splitting schemes:     
    - Roe (Approximate Riemann Solver)
    - AUSM (Jameson)
    - CUSP (Jameson)
    - HLL (Hartman-Hayman-Lax)
    - HLLC (Hartman-Hayman-Lax contact with Toro and Batten variants)
    - OS (Osher-Solomon)
    
Implemented High Resolution Flux-splitting schemes:     
    - LW (Lax-Wendroff)
    - JST (Jameson-Schmidt-Turkel)
    
Implemented Viscous Flux-splitting schemes:
    - NS (explicit viscous fluxes for laminar compressible flows) *** TODO *** Test and debugging    

Implemented EntropyFix schemes:
    - Roe: a) HH1 (Harten-Hayman #1)
           b) HH2 (Harten-Hayman #2)
           c) HH2b (Harten-Hayman #2 with local Mach number correction by V. Selmin)
    - HLL(C): a) D1 (Davies #1)
              b) D2 (Davies #2)
              c) ER (Enfieldt-Roe)
              d) LV (LeVeque)
              e) Toro (only for HLLC scheme)        

Implemented Flux-Limiters schemes:
    - VL (Lax-Wendroff with VanLeer limiter)  
    - MinMod (Lax-Wendroff with MinMod limiter) 
    - SB (Lax-Wendroff with SuperBee limiter) 
    - MC (Lax-Wendroff with MonotonizedCentral limiter) 

Implemented norm for residual evaluation:
    - L1  = sum( |x_i| )/Nv
    - Loo = max( |x_i| )  
    
Implemented Timestepping schemes:
    - EE
    - RK2
    - RK3
    - RK4
    - RK4LS (low storage Runge-Kutta scheme, perhaps only for steady-state)   

Implemented Local Timestepping scheme:     
    - Speed up convergence to steady-state with dt = dt(CFL_i) such that CFL_i ~= CFL_MAX
      REMARK: CFL_i is linearly increased to CFL_MAX during the first Nsteps time steps to increase robustness
      REMARK: If local timestepping is active simulation is ended when residuals are less than minResidual
      WARNING: Local timestepping strategy may cause instabilities if 2nd order is active
    
Implemented boundary conditions:
    - supersonicInlet
    - extrapolatedOutlet
    - slip, symmetryPlane, inviscidWall
    - Riemann
    - transpiration 
    - viscousAdiabaticWall, viscousIsothermalWall
    - extrapolatedP, inlet
    - extrapolatedT, adiabaticWall
    - extrapolatedU
    - extrapolatedPT
    - extrapolatedPU
    - extrapolatedTU
    - swirl
    - Empty   
    The ghost cell solution is extrapolated:
       - extrapolateBC = 0 uniform extrapolation of solution on ghost cells
       - extrapolateBC = 1 linear extraplation of solution on ghost cells
    REMARK: wall boundary conditions are now re-enforced prescribing the numerical fluxes vector   
          
Toolboxes and utilities:
    - CAE2D = Computational Aero Elasticity Toolbox for 2D 2/3 d.o.f. (h, a, d) airfoil sections  
              > Aerodynamic loads are computed and printed ( with linear extrapolation of cell pressure )
              > A forced motion is imposed with transpiratin BCs
              > ODEs of free motion for h and a are solved, while aileron deflection d is imposed 
              > ODE of free motion for d is solved, for example for buzz numerical simulation *** TODO *** Test inviscid and viscous buzz
                
    - CAE3D = Computational Aero Elasticity Toolbox for 3D wings (with NAEMO and MASSA postprocessing)
              > REMARK: arrays known on a structural or aerodynamic mesh are denoted with (s) or (a)
              > Trimmed solution (both aerodynamic and structural) is computed with an aeroelastic
                modal iterative method ( [diag(m_i*w0_i^2)]{q}=[U]'{F_s} )
              > Modal movement (q_i=q_i(t)) is imposed (a step smoothed with 1-cos function)
              > Structural mesh movement and velocity ({u_s}=[U]{q}, {up_s}=[U]{qp}) are computed
              > Aerodynamic mesh movement and velocity ({u_a} and {up_a}) are interpolated with
                an energy preserving interface (such that {u_a}=[T]{u_s} and {F_s}=[T]'{F_a})
                thanks to id_a2s connectivity array 
              > Transpiration boundary conditions are enforced (the completely non linear formulation
                is used so that Vn = {V}'{Dn} + {Vb}'{n0 + Dn} rather than Vn = {V}'{Dn} + {Vb}'{n0})
              > 3D Euler solver    
              > Aerodynamic loads ({F_a}) are computed ( with linear extrapolation of cell pressure ) 
                and interpolated on the structural mesh ({F_s}) 
              > Generalized forces ({Q_s}=[U]'{F_s}) are computed and stored for postprocessing with NAEMO    
              
    - showDisplacement = Deformed shape (both 2D and 3D) can be plotted in paraFoam with this utility, which 
                         moves the mesh and also computes the volumes swept by faces (for ALE formulation)
                         > with -amp option a scale factor is applied (it does not work correctly for rotations)
                         > WARNING!!! cellDisplacement should be initialized with fixedValue BCs on the boundaries
                           to be deformed, otherwise deformation won't work
      
Other:
    - testInput = [0, 1] to verify input parameters before execution  
    - testMesh  = [0, 1] to verify connectivity in searching RR and LL cells (plot in ParaFoam cell T field)
    - adaptiveConnectivity = [0, 1] to change the mesh search algorithm used to build the extended cells connectivity
                             REMARK: The Adaptive algorithm (1) works fine with structured meshes while the 
                             NonAdaptive algorithm (0) is to be preferred for tetrahedral meshes
    - loadConnectivity = [-1, 0, 1, 2] to write on file Log/LR2LLRR.txt the connectivity matrix id_LL/RR_i if (-1),
                         to load from file Data/LR2LLRR.txt the MDS corrected connectivity matrix id_LL/RR_i if (1),
                         to do both if (2) 
                         REMARK: Using a MDS corrected connectivity matrix sohuld increase accuracy (e.g. ONERA M6 wing)
    - testInterface = [0, 1] to verify interface between aerodynamic an structural meshes in CAE3D (plot in ParaFoam cell T field) 
    - AeroFoam4Speed = this parser comments out unused options to optimize efficiency, compiles executables
                       and restores original source code files. Useful for automatic use of Toolboxes.
    
ToDo:
 - HARD
    - Fully implicit with (modified) NR method (an interface for lduMatrix must be built)
    - Node Centered vs. Cell Centered (Voronoi polyhedral mesh from a FE tetrahedral mesh)
    - Mixture of reacting real gases in thermodynamic equilibrium
    - Add viscous fluxes (and turbulence models) 
    - Axisymmetric solver without swirl (quasi-2D)
    - Parallelization (with OpenFOAM built-in functions)
    - ALE formulation (with OpenFOAM built-in functions)
 - EASY
    - CAE3D Structural static deformation is computed with Code_Aster FEM code   
    - CAE3D Add Wind Axes definition in input
    - CAE3D Load structural model data only if necessary
    - CAE3D Add modal scale factor in input ( in the first lines of .mode# files )
    - CAE3D Add rigid motion parameters in input
    
Known bugs: 
 - LL and RR cells internal search algorithm fails for particularly badly shaped volumes 
   (it should be fixed loading connectivity from file or with adaptive search algorithm)

\*----------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // System include

#   include "fvCFD.H"
#   include "meshSearch.H"
#   include "pointMesh.H"
#   include "pointFields.H"
#   include "Matrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // AeroFoam include
                                   
//---------------------------------// *** Flux splitting schemes ***
#   include "RoeFlux.C"            // Roe
#   include "AUSMFlux.C"           // Advection Upwind Splitting Method 
#   include "CUSPFlux.C"           // Convective Upwind and Split Pressure
#   include "HLLFlux.C"            // Harten-Lax-vanLeer
#   include "HLLCFlux.C"           // Harten-Lax-vanLeer Corrected 
#   include "OSFlux.C"             // Osher-Solomon
// #   include "NewFlux.C"         // Template to add a new flux splitting scheme

//---------------------------------// *** High resolution flux splitting schemes ***
#   include "LWHRFlux.C"           // Roe/Lax-Wendroff HR
#   include "JSTHRFlux.C"          // Jameson-Schmidt-Tuckel HR
#   include "ViscousFlux.C"        // Viscous Fluxes (laminar)
// #   include "NewFlux.C"         // Template to add a new HR flux splitting scheme

//---------------------------------// *** AeroFoam routines *** 
#   include "BuildFluxes.C"        // Build numerical fluxes
#   include "UpdateSolution.C"     // Update solution

//---------------------------------// *** CAE Toolbox *** 
#   include "CAE2DTools.C"         // 2D aeroelasticity (airfoil)
#   include "CAE3DTools.C"         // 3D aeroelasticity (wing)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//============================================================================//
//                                                                            //
//                                *** MAIN ***                                //
//                                                                            //
//============================================================================//
int main(int argc, char *argv[])
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // System include

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // AeroFoam include
 
#   include "ReadThermodynamicProperties.C"
#   include "CreateFields.C"
#   include "CreateConnectivity.C"
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // ComputationalAeroElasticity Toolboxes

//IF{ toolBox CAE2D
//#   include "CAE2DInit.C"
//FI} toolBox CAE2D
//IF{ toolBox CAE3D
//#   include "CAE3DInit.C"
//FI} toolBox CAE3D

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
            
    // Variables definition
    label k = 0;                                               // Dimensions
    label nV = mesh.V().size(), nS = mesh.Sf().size();         // 
    scalar t0cpu = 0.0, t1cpu = 0.0, ETA = 0.0, ERT = 0.0;     // Time statistics
    scalar dtcpu = 0.0, dtsim = 0.0, t1sim = 0.0, Co = 0.0;    //
    scalar aa1 = 0.0, aa2 = 0.0, aa3 = 0.0, aa4 = 0.0;         // RK coefficient
    scalar err_rrho = 0.0, err_mm = 0.0, err_eet_tilde = 0.0;  //  
    FILE* ff;                                                  // Statistics  
    
    // Epsilon
    dimensionedScalar epsRho      = dimensionedScalar("epsRho"     ,      rrho.dimensions(), 1e-10);
    dimensionedScalar epsM        = dimensionedScalar("epsM"       ,        mm.dimensions(), 1e-10);    
    dimensionedScalar epsEt_tilde = dimensionedScalar("epsEt_tilde", eet_tilde.dimensions(), 1e-10); 
          
    // Read time-stepping scheme
    // timeScheme = [ EE, AB2, RKp ]
    word timeScheme = word("EE");
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("timeScheme") ) 
    timeScheme = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("timeScheme") ); 

    // Read CourantMAX parameters
    // CoMAX = [ 0.0, depending on the time integration scheme (e.g. CoMAX = 2.0 for RK2) ]
    scalar CoMAX = 0.0, CoMAX_ = 0.0;
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("CourantMax") )
    CoMAX = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("CourantMax") ); 
    scalar Nsteps = 100.0;
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("CourantSteps") )
    Nsteps = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("CourantSteps") ); 

    // Read monotone flux splitting scheme
    // fs = [ Roe, OS, AUSM, CUSP, HLL, HLLC, KNP ]
    word mf = word( "Roe" );
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("MonotoneFlux") )
    mf = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("MonotoneFlux") ); 
        
    // Read high resolution flux splitting scheme
    // hr = [ none, LW, JST ]
    word hr = word( "LW" );
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("HighResolutionFlux") )
    hr = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("HighResolutionFlux") ); 
    
    // Read flux limiter
    word fl = word( "VL" ); 
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("fluxLimiter") )
    fl = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("fluxLimiter") );  
    
    // Read entropy fix
    // ef = [ Roe:{ HH1, HH2, HH2b, none }, HLL/C{ D1, D2, ER, LV, Toro } ] 
    word ef = word("HH2");
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("entropyFix") )
    ef = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("entropyFix") ); 

    // Read extrapolateBC
    // eBC = [ 0.0, 1.0 ] 
    scalar eBC = 0.0;
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("extrapolateBC") )
    eBC = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("extrapolateBC") );      
           
    // Read residualNorm
    // resNorm = [ L1, Loo ] 
    word resNorm = word("L1");
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("residualNorm") )
    resNorm = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("residualNorm") ); 
    
    // Read minResidual
    scalar minResidual = -1;
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("minResidual") )
    minResidual = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("minResidual") ); 
    
    // Read flowType
    // flowType = [ Euler, NS, ( *** TODO: RANS *** ) ] 
    word ft = word("Euler");
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("flowType") )
    ft = word( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("flowType") ); 
        
    // Read testInput
    // testInput = [ 0.0:1.0 ]
    scalar test_input = 0.0;
    if ( mesh.schemesDict().subDict("AeroFoamSchemes").found("testInput") )
    test_input = readScalar( mesh.schemesDict().subDict("AeroFoamSchemes").lookup("testInput") ); 
        
    // Print mesh statistics
    Info << "========================================"           << nl;
    Info << " Mesh "                                             << nl;
    Info << "========================================"           << nl;
    Info << " # of Faces    =  " << nS                           << nl;
    Info << " # of Cells    =  " << nV                           << nl;
    Info << " # of BndPatch =  " << mesh.boundaryMesh().size()   << nl;
    forAll( mesh.boundaryMesh(), iPatch ) 
    {
       Info << " " << iPatch << ") " << mesh.boundaryMesh().physicalTypes()[iPatch] << " @ " << mesh.boundaryMesh().names()[iPatch] << nl;
    }
    Info << "----------------------------------------"           << nl << nl;
    
    // Print numerical method statistics
    Info << "========================================"           << nl;
    Info << " Numerical Scheme "                                 << nl;
    Info << "========================================"           << nl;
    Info << " flowType      =  " << ft                           << nl;
    Info << " timeScheme    =  " << timeScheme                   << nl;
    Info << " CourantMax    =  " << CoMAX                        << nl;
    Info << " CourantSteps  =  " << Nsteps                       << nl;
    Info << " MonotoneFlux  =  " << mf                           << nl;
    Info << " HighResFlux   =  " << hr                           << nl;
    Info << " fluxLimiter   =  " << fl                           << nl;
    Info << " entropyFix    =  " << ef                           << nl;
    Info << " extrapolateBC =  " << eBC                          << nl;
    Info << " residualNorm  =  " << resNorm                      << nl;
    Info << " minResidual   =  " << minResidual                  << nl;
    Info << "----------------------------------------"           << nl;  
        
    // Check for input errors 
    if ( test_input == 1 )
    {
        Info << " Does everything look fine? [y] to continue [n] otherwise ";
        if ( getc(stdin) == 'n' ) return(0);
    }
    
    // Initialize convergence history files (clear only if simulation restarts from t = 0)
    if ( runTime.value() == 0 )
    {
        ff = fopen("./Log/residuals.txt", "w");
        fclose(ff);
    }
    
    //==========================================================================
    // *** TIME LOOP [index (k)] ***
    //==========================================================================
    Info << "\n Starting time loop...\n";    
    for (runTime++; !runTime.end(); runTime++)
    {
    
        // Time    
        t1sim  = runTime.value();
        dtsim  = min( localDt );
        
        //======================================================================
        // Explicit Euler temporal integration scheme
        //======================================================================
        //IF{ timeScheme EE
        if ( timeScheme == "EE" )
        {        
            // a) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho, &mm, &eet_tilde, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // b) Update solution
            update_solution( 1.0, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );

        }
        //FI} timeScheme EE
        //======================================================================
        // (Explicit) 2 substeps Runge-Kutta temporal integration scheme
        //======================================================================
        //IF{ timeScheme RK2
        else if ( timeScheme == "RK2" )
        {
            // Substep 1
            // 1.a) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_o, &mm_o, &eet_tilde_o, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
            
            // 1.b) Update solution
            update_solution( 0.5, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
 
            // Substep 2
            // 2.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( 1.0, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 2.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 2.c) Update solution
            update_solution( 0.5, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
        }
        //FI} timeScheme RK2
        //======================================================================
        // (Explicit) 3 substeps Runge-Kutta temporal integration scheme
        //======================================================================
        //IF{ timeScheme RK3
        else if ( timeScheme == "RK3" )
        {
            // Substep 1
            // 1.a) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_o, &mm_o, &eet_tilde_o, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 1.b) Update solution
            update_solution( 1.0/6.0, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
 
            // Substep 2
            // 2.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( 0.5, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 2.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 2.c) Update solution
            update_solution( 4.0/6.0, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
            // Substep 3
            // 3.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( 2.0, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 3.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_o, &mm_o, &eet_tilde_o, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 3.c) Update solution             
            update_solution( -1.0, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );       
            
            // 3.d) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 3.e) Update solution
            update_solution( 1.0/6.0, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
        }
        //FI} timeScheme RK3
        //======================================================================
        // (Explicit) 4 substeps Runge-Kutta temporal integration scheme
        //======================================================================
        //IF{ timeScheme RK4
        else if ( timeScheme == "RK4" )
        {
            // RK Coefficients
            aa1 = 1.0/6; aa2 = 1.0/3; aa3 = 1.0/3; aa4 = 1.0/6;  
            
            // Substep 1
            // 1.a) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_o, &mm_o, &eet_tilde_o, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 1.b) Update solution
            update_solution( aa1, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
            // Substep 2
            // 2.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( 0.5, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 2.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 2.c) Update solution
            update_solution( aa2, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
            // Substep 3
            // 3.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( 0.5, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 3.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 3.c) Update solution
            update_solution( aa3, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
            // Substep 4
            // 4.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( 1.0, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 4.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // 4.c) Update solution
            update_solution( aa4, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
            
        }
        //FI} timeScheme RK4
        //======================================================================
        // (Explicit) 4 substeps Runge-Kutta temporal integration scheme with 
        // modified coefficient to minimize dissipation ( from Hirsch )
        //======================================================================
        //IF{ timeScheme RK4LS
        else if ( timeScheme == "RK4LS" )
        {
        
            // RK coefficients from slide-show
            aa1 = 1.0/4;  aa2 = 1.0/3;  aa3 = 1.0/2;  aa4 = 1.0; 
            
            // Substep 1
            // 1.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( aa1, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
                     
            // 1.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // Substep 2
            // 2.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( aa2, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 2.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // Substep 3
            // 3.a) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( aa3, &localDt, &mesh, &rrho_oo, &mm_oo, &eet_tilde_oo, &FFrho, &FFm, &FFet_tilde );
            
            // 3.b) Build fluxes
            build_fluxes( &localDt, ft, &R, &Cv, &Pr, &mesh, &id_LL_i, &id_RR_i,
                          &rrho_oo, &mm_oo, &eet_tilde_oo, 
                          &p, &U, &T, &gradUx, &gradUy, &gradUz, &gradT,
                          &FFrho, &FFm, &FFet_tilde, 
                          mf, hr, fl, ef, eBC );
                          
            // Substep 4
            // 4) Update solution
            rrho_oo      = rrho_o; 
            mm_oo        = mm_o;
            eet_tilde_oo = eet_tilde_o;
            update_solution( aa4, &localDt, &mesh, &rrho, &mm, &eet_tilde, &FFrho, &FFm, &FFet_tilde );
        
        }
        //FI} timeScheme RK4LS
               
        //Write results and convergence history on file (L1 norm vs Loo norm)
        //IF{ residualNorm Loo 
        if ( resNorm == "Loo" )
        {
            err_rrho      = max( mag( rrho      - rrho_o      )/( mag( rrho_o      ) + epsRho      ) ).value();
            err_mm        = max( mag( mm        - mm_o        )/( mag( mm_o        ) + epsM        ) ).value();
            err_eet_tilde = max( mag( eet_tilde - eet_tilde_o )/( mag( eet_tilde_o ) + epsEt_tilde ) ).value();
        }
        //FI} residualNorm Loo
        //IF{ residualNorm L1 
        else if ( resNorm == "L1" )
        {
            err_rrho      = sum( mag( rrho      - rrho_o      )/( mag( rrho_o      ) + epsRho      ) ).value()/nV;
            err_mm        = sum( mag( mm        - mm_o        )/( mag( mm_o        ) + epsM        ) ).value()/nV;
            err_eet_tilde = sum( mag( eet_tilde - eet_tilde_o )/( mag( eet_tilde_o ) + epsEt_tilde ) ).value()/nV;
        }
        //FI} residualNorm L1 
        //IF{ residualNorm L2 
        else if ( resNorm == "L2" )
        {
            err_rrho      = Foam::sqrt( sum( magSqr( rrho      - rrho_o      ) ).value() )/Foam::sqrt( sum( ( magSqr( rrho_o      ) + epsRho      ) ).value() )/nV;
            err_mm        = Foam::sqrt( sum( magSqr( mm        - mm_o        ) ).value() )/Foam::sqrt( sum( ( magSqr( mm_o        ) + epsM        ) ).value() )/nV;
            err_eet_tilde = Foam::sqrt( sum( magSqr( eet_tilde - eet_tilde_o ) ).value() )/Foam::sqrt( sum( ( magSqr( eet_tilde_o ) + epsEt_tilde ) ).value() )/nV;
        }
        //FI} residualNorm L2 
        //IF{ writeResidual 1
        ff = fopen("./Log/residuals.txt", "a");
        fprintf(ff, "%e %e %e\n", err_rrho, err_mm, err_eet_tilde);
        fclose(ff);
        //FI} writeResidual 1
             
        // Update solution at time level (k-1)->oo and (k)->o 
        rrho_oo      = rrho_o;
        mm_oo        = mm_o;
        eet_tilde_oo = eet_tilde_o;
        rrho_o       = rrho;
        mm_o         = mm;
        eet_tilde_o  = eet_tilde;
        U            = mm/rrho;
        p            = ( gamma - 1 )*( eet_tilde - 0.5*rrho*magSqr(U) );
        T            = p/( rrho*R );
        
        // Courant Number   
        localCo = 0.0;
        forAll( mesh.Sf(), i )
        {
            label id_o = mesh.faceOwner()[i];
            label id_n = mesh.faceNeighbour()[i];
            localCo[ id_o ] = localCo[ id_o ] + ( ( mag( U[id_o] ) + Foam::sqrt( gamma.value()*R.value()*T[id_o] ) )*mesh.magSf()[i]/mesh.V()[id_o] )*localDt[id_o];
            localCo[ id_n ] = localCo[ id_n ] + ( ( mag( U[id_n] ) + Foam::sqrt( gamma.value()*R.value()*T[id_n] ) )*mesh.magSf()[i]/mesh.V()[id_n] )*localDt[id_n];
        }
        Co = max( localCo ); 
        
        // TimeStep Adaptivity
        CoMAX_ = ( k + 1 )*CoMAX/Nsteps*( k < Nsteps ) + CoMAX*( k >= Nsteps );
        if ( CoMAX > 0 )
        {
            forAll( mesh.V(), i )
            {
                localDt[i] = localDt[i]*CoMAX_/( localCo[i] + 1e-5 );
            }
        }
        else
        {
            localDt = runTime.deltaT().value();
        }
                       
        // Print statistics
        runTime.write();
        t1cpu  = runTime.elapsedCpuTime(); 
        dtcpu  = ( t1cpu - t0cpu );                      
        ERT    = runTime.elapsedClockTime();
        ETA    = ( runTime.endTime().value() - t1sim )/( dtsim + 1e-10 )*ERT/( k + 1 );
        Info                                               << nl;                                          
        Info << "========================================" << nl;
        Info << " Iteration # " << k << " @ " << argv[2]   << nl;
        Info << "========================================" << nl;
        Info << "TimeStep      [s] = " << dtsim            << nl;
        Info << "CourantNumber [-] = " << Co               << nl;
        Info << "RunTime       [s] = " << t1sim            << nl;  
        Info << "IterationTime [s] = " << dtcpu            << nl;
        Info << "ExecutionTime [h] = " << ERT/3600         << nl;
        Info << "RemainingTime [h] = " << ETA/3600         << nl;
        Info << "----------------------------------------" << nl;
        Info << "ContinuityEq  [-] = " << err_rrho         << nl;
        Info << "MomentumEq    [-] = " << err_mm           << nl;
        Info << "EnergyEq      [-] = " << err_eet_tilde    << nl;
        Info << "----------------------------------------" << nl;
        t0cpu = t1cpu;  
        k     = k + 1;
        
        // Check residual and stop simulation if Residual_i < minResidual for all i
        if ( ( k > 100 ) & ( err_rrho < minResidual ) & ( err_mm < minResidual ) & ( err_eet_tilde < minResidual ) )
        {
            Info << "Convergence to minResidual reached!" << nl;
            Info << "AeroFoam terminated succesfully!" << endl;
            system("date");
            return(0);
        }
        
        //======================================================================
        // 2D aero-elastic airfoil ( see Bisplinghoff pag. 533 )
        //======================================================================
//IF{ toolBox CAE2D
//#       include "CAE2DTimeStep.C"    
//FI} toolBox CAE2D
        
        //======================================================================
        // 3D wing flutter analysis ( see Cavagna et al. )
        //====================================================================== 
//IF{ toolBox CAE3D       
//#       include "CAE3DTimeStep.C" 
//FI} toolBox CAE3D

    }    

    // End
    Info << "AeroFoam terminated succesfully!" << endl;
    system("date");
    return(0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

