//==============================================================================
// 2D aero-elastic airfoil ( see Bisplinghoff pag. 533 )
//==============================================================================
{     
    // Compute aerodynamic loads  
    CAE2DcomputeCp( &p, Poo_,  qoo_, &mesh, &nn0, &id_bodyCell, &id_bodyFace, &id_LL_i, &CCp );   
    CAE2DcomputeLoads( &xx, &SS, &nn, &CCp, i_hat_, j_hat_, k_hat_, Cref, Lref, Xref, Xhinge, &C_L, &C_D, &C_MX, &C_H );
          
    // Print statistics
    Info << "C_Lift        [-] = " << C_L                       << nl;
    Info << "C_Drag        [-] = " << C_D                       << nl;
    Info << "C_MomentumX   [-] = " << C_MX                      << nl;
    Info << "C_Hinge       [-] = " << C_H                       << nl;
    Info << "----------------------------------------"          << nl;
    CAE2DprintLoads( t1sim, C_L, C_D, C_MX, C_H );
    CAE2DprintCp( &xx, &CCp, Cref );
      
//..............................................................................FORCED MOTION  
if ( solutionType == 1 )
{     
    // Forced motion
    CAE2DforcedMotion( Mov_h, A0_h, A1_h, f_h, tau_h, t1sim, &h, &h_p );
    CAE2DforcedMotion( Mov_a, A0_a, A1_a, f_a, tau_a, t1sim, &a, &a_p );
    CAE2DforcedMotion( Mov_d, A0_d, A1_d, f_d, tau_d, t1sim, &d, &d_p );
}    
//..............................................................................FREE MOTION (h, a)
else if ( solutionType == 2 )  
{          
    // Solve ODEs of 2 d.o.f aeroelastic airfoil ( ~= Bisplinghoff pag. 533, BACT model ) 
    // TODO: Add d d.o.f. residualized as input (but with hinge moment equation to be checked)
    // REMARK: a > 0 if airfoil is pitched up
    CAE2DsolveODE_ha( m, d_AECG, I_CG, K_hh, K_aa, qoo_, Cref, Lref, C_L, C_MX, dtsim, &h, &a, &h_p, &a_p );
} 
//..............................................................................BUZZ (d)
else if ( solutionType == 3 )  
{
    // REMARK: d > 0 if aileron is pitched down increasing airfoil camber
    CAE2DsolveODE_d( I_dd, C_dd, K_dd, qoo_, Cref, Lref, C_H, dtsim, &d, &d_p );
}
//..............................................................................TRANSPIRATION BCs
  
    // Compute transpiration velocity BCs   
    CAE2DcomputeNormals( h, a, d, &xx, Xref, Xhinge, &nn0, &nn );
    CAE2DcomputeDisplacements( h, a, d, &xx, Xref, Xhinge, &uu );
    CAE2DcomputeVelocities( h_p, a_p, d_p, &xx, Xref, Xhinge, &uu_p ); 
    CAE2DcomputeTranspirationVelocity( &uu_p, &nn0, &nn, &U, &id_bodyCell, &VVbn );
    CAE2DsetTranspirationVelocityBCs( &mesh, id_bodyPatch, &U, &VVbn, &nn );
      
    // Print statistics
    Info << "h(t)          [m] = " << h                         << nl;
    Info << "h_p(t)      [m/s] = " << h_p                       << nl;
    Info << "a(t)        [deg] = " << a*180.0/pi                << nl;
    Info << "a_p(t)    [deg/s] = " << a_p*180.0/pi              << nl;
    Info << "d(t)        [deg] = " << d*180.0/pi                << nl;
    Info << "d_p(t)    [deg/s] = " << d_p*180.0/pi              << nl;
    Info << "----------------------------------------"          << nl;  
    CAE2DupdateCellDisplacement( id_bodyPatch, &cellDisplacement, &uu );   
    CAE2DprintMovement( t1sim, h, a, d, h_p, a_p, d_p ); 
                
}
