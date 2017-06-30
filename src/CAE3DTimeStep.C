//==============================================================================
// 3D aero-elastic wing ( see Cavagna et al. )
//==============================================================================
{
//..............................................................................AERODYNAMIC LOADS   
   
    // Compute aerodynamic loads 
    CAE2DcomputeCp( &p, Poo_, qoo_, &mesh, &nn0_a, &id_bodyCell, &id_bodyFace, &id_LL_i, &CCp_a );
    CAE3DcomputeLoads( &xx_a, &SS_a, &nn_a, &CCp_a, i_wind_, j_wind_, k_wind_, Sref_a, Lref_a, Xref_a, 
                       &C_L, &C_D, &C_FX, &C_FY, &C_FZ, &C_MX, &C_MY, &C_MZ );
                       
    // Print statistics
    Info << "C_Lift        [-] = " << C_L                       << nl;
    Info << "C_Drag        [-] = " << C_D                       << nl;
    Info << "C_ForceX      [-] = " << C_FX                      << nl;
    Info << "C_ForceY      [-] = " << C_FY                      << nl;
    Info << "C_ForceZ      [-] = " << C_FZ                      << nl;
    Info << "C_MomentumX   [-] = " << C_MX                      << nl;
    Info << "C_MomentumY   [-] = " << C_MY                      << nl;
    Info << "C_MomentumZ   [-] = " << C_MZ                      << nl;
    Info << "----------------------------------------"          << nl; 
    if( t1sim >= tPrint ) 
    {
        CAE3DprintLoads( t1sim, C_FX, C_FY, C_FZ, C_MX, C_MY, C_MZ );
        CAE3DprintCp( &xx_a, &CCp_a  );
    }
       
//..............................................................................GENERALIZED AERODYNAMIC FORCES          
///*                  
    // Forces interface     
    CAE3DcomputeR_a( qoo_, &SS_a, &nn_a, &CCp_a, &RR_a );    
    CAE3DcomputeR_s_v( &Tas, &RR_a, &RR_s_v );  

    // Compute generalized aerodynamic forces 
    CAE3DcomputeQ_s( &UUx_s, &UUy_s, &UUz_s, &RR_s_v, &QQ_s );  
    if ( t1sim <= ( t0Forcing_s + dtsim ) ) QQ0_s = QQ_s; // After the first forcing iteration save QQ0_s = QQ_s as a reference 

    // Check convergence on modal displacements and generalized aerodynamic forces
    if ( resNorm == "Loo" )
    {
        err_qq_s = max( qq_s - QQo_s )/max( qqo_s + 1e-10 )/Nmodes;
        err_QQ_s = max( QQ_s - QQo_s )/max( QQo_s + 1e-10 )/Nmodes;
    }
    else 
    {
        err_qq_s = sum( mag( qq_s - qqo_s ) )/( sum( mag( qqo_s ) ) + 1e-10 )/Nmodes;
        err_QQ_s = sum( mag( QQ_s - QQo_s ) )/( sum( mag( QQo_s ) ) + 1e-10 )/Nmodes;
    }
    qqo_s = qq_s;
    QQo_s = QQ_s;  
    
    // Print statistics
    Info << "tau           [-] = " << ( t1sim - t0Forcing_s )*Uoo_/Lref_a                         << nl;
    forAll(QQ_s, ii) Info << "qs_" << ii+1 << "          [m] = " << qq_s[ii]                      << nl;
    Info << "Residual_qs   [-] = " << err_qq_s                                                    << nl;  
    forAll(QQ_s, ii) Info << "Qa_" << ii+1 << "        [m^2] = " << ( QQ_s[ii] - QQ0_s[ii] )/qoo_ << nl;
    Info << "Residual_Qa   [-] = " << err_QQ_s                                                    << nl;   
    Info << "----------------------------------------"                                            << nl; 
    if( t1sim >= tPrint ) 
    {
        CAE3DprintNAEMOWrk( NactiveMode, qoo_, &qq_s, &QQ0_s, &QQ_s, k );
        CAE3DprintGenDisplacement( t1sim - t0Forcing_s, Uoo_, Lref_a, &qq_s );
        CAE3DprintGenForces( t1sim - t0Forcing_s, Uoo_, Lref_a, qoo_, &QQ0_s, &QQ_s );
    }       
//*/
//..............................................................................FORCED RIGID MOTION
/*
    // Forced rigid motion (useful to check transpiration velocity BCs validity)
    scalar pi = 3.14159265358979;
    label Mov_h, Mov_a;
    vector h_hat, a_hat;
    scalar A0_h, A1_h, f_h, tau_h, h = 0.0, hp = 0.0;
    scalar A0_a, A1_a, f_a, tau_a, a = 0.0, ap = 0.0;
      
    // Forced motion input parameters
    // Translation h of the EA
    Mov_h = 0;
    h_hat = vector(0.0, 0.0, -1.0);
    A0_h  = 0.0;
    A1_h  = 0.0;
    f_h   = 0.0;
    tau_h = 0.0;
    // Rotation alpha around the EA
    Mov_a = 1;
    a_hat = vector(0.0, 1.0, 0.0);
    A0_a  = 6.0*pi/180;
    A1_a  = 0.0*pi/180;
    f_a   = 0.0;
    tau_a = 0.0;          
     
    // Forced motion
    CAE3DforcedRigidMotion( Mov_h, A0_h, A1_h, f_h, tau_h, t1sim, &h, &hp );
    CAE3DforcedRigidMotion( Mov_a, A0_a, A1_a, f_a, tau_a, t1sim, &a, &ap );
    
    // Update normals and body velocity    
    CAE3DcomputeRigidNn_a( h, a, hp, ap, h_hat, a_hat, &xx_a, Xref_a, &nn0_a, &nn_a );
    CAE3DcomputeRigidUup_a( h, a, hp, ap, h_hat, a_hat, &xx_a, Xref_a, &uup_a );   
    if( t1sim >= tPrint ) CAE3DprintRigidMode( h, a, hp, ap, h_hat, a_hat, &xx_s_v, Xref_a );
*/    
//..............................................................................FORCED MODAL MOTION
///*
    // Only aerodynamics and structure still
    if ( solutionType == 0 )
    {
        qq_s  = 0.0;
        qqp_s = 0.0;
    }
    // Static modal solution or iterative aeroelastic solver (e.g. for trim solution)
    else if ( solutionType == 1 )
    {
        CAE3DstaticModalSolution( &mm_s, &cc_s, &ff_s, &QQ_s, &qq_s, &qqp_s );
    }
    // Dynamic modal solution non linear coupling (e.g. for postflutter direct coupled nonlinear simulation)
    // *** Fully implicit solver mandatory to obtain significant results in a reasonable amount of time ***
    else if ( solutionType == 2 )
    {
        CAE3DdynamicModalSolution( &mm_s, &cc_s, &ff_s, &QQ_s, dtsim, &qq_s, &qqp_s );
    }
    // Forced modal motion with a blended step to build linearized ROM [Ham(k, M, Re)] with NAEMO 
    else if ( solutionType == 3 )
    {
        CAE3DforcedModalMotion( NactiveMode, qMax_s, tauMax_s, Uoo_, Lref_a, t1sim, t0Forcing_s, &qq_s, &qqp_s );
    }
    // Forced rigid pitch mode motion (yet another way to impose a forced rigid motion)
    // REMARK: - The rigid motion mode is a combination of pitching (a) and 
    //           plunging (h) and can be written with function CAEprintRigidMode
    //         - The rigid mode file StructuralModel.mode0 should be renamed as 
    //           StructuralModel.mode<Nmodes+1> ( Nmodes and NactiveMode should be 
    //           then modified accordingly )
    //         - This option is also used to check the aerodynamic matrix static gain [Ham(0)]
    else if ( solutionType == 4 )
    {
        qq_s  = 0.0;
        qqp_s = 0.0;
        qq_s[NactiveMode-1] = qMax_s;
    }
    // Save solution for restart
    CAE3DprintModalState( &qq_s, &qqp_s ); 
    
    // Structural displacement and velocity
    CAE3DcomputeUu_s_v( &UUx_s, &UUy_s, &UUz_s, &qq_s, &uu_s_v );
    CAE3DcomputeUu_s_v( &UUx_s, &UUy_s, &UUz_s, &qqp_s, &uup_s_v );
        
    // Displacement interface
    CAE3DcomputeNn_s_e( &xx_s_v, &uu_s_v, &id_A_s, &id_B_s, &id_C_s, &nn_s_e );
    CAE3DcomputeNn_a( &id_a2s, &id_A_s, &id_B_s, &id_C_s, &xx_a, &xx_s_v, &nn0_s_e, &nn_s_e, &nn0_a, &nn_a );
    CAE3DcomputeUu_a( &Tas, &uup_s_v, &uup_a );  
    if( t1sim >= tPrint ) 
    {  
        CAE3DprintDeformedShape( &xx_s_v, &uu_s_v, &id_A_s, &id_B_s, &id_C_s );
        CAE3DcomputeUu_a( &Tas, &uu_s_v, &uu_a ); 
        //------------------------------------------------------------------------------ 
        // REMARK: Always check that at staringTime cellDisplacement BCs are fixedValue
        //------------------------------------------------------------------------------ 
        CAE3DupdateCellDisplacement( id_bodyPatch, &cellDisplacement, &uu_a );
    }    
//*/
//..............................................................................TRANSPIRATION VELOCITY BCs
///*
    // Compute and set transpiration BCs
    CAE3DcomputeTranspirationVelocity( &uup_a, &nn0_a, &nn_a, &U, &id_bodyCell, &VVbn_a );
    CAE3DsetTranspirationVelocityBCs( &mesh, id_bodyPatch, &U, &VVbn_a, &nn_a );
    if( t1sim >= tPrint ) tPrint = tPrint + dtPrint;
//*/        
//..............................................................................END    
}
