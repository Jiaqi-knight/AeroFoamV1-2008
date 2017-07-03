//=========================================================================//
//                                                                         //
//                           *** COMPUTE_CP ***                            //
//                                                                         //
//=========================================================================//
int CAE2DcomputeCp( volScalarField *P, scalar Poo, scalar qoo, 
                    fvMesh *mesh, vectorField * nn, labelField *id_bodyCell, labelField *id_bodyFace, labelField *id_LL_i,
                    scalarField *CCp )
{
    // Variables definition
    label i, id_L, id_LL;
    scalar p_L, p_LL, p_R, dx;
    
    // Loop on boundary faces
    forAll( (*id_bodyCell), ii )
    {
        id_L  = (*id_bodyCell)[ii];
        // Linear extrapolation (when possible) of pressure field on the body surface
        i     = (*id_bodyFace)[ii];
        id_LL = (*id_LL_i)[i];
        dx    = mag( ( (*mesh).C()[id_L] - (*mesh).C()[id_LL] ) & (*nn)[ii] )/
                mag( ( (*mesh).Cf()[i]   - (*mesh).C()[id_LL] ) & (*nn)[ii] );
        p_L   = (*P)[id_L];
        p_LL  = (*P)[id_LL];
        p_R   = p_L;
        if ( dx > 0.5 ) p_R = p_LL + ( p_L - p_LL )/dx;
        (*CCp)[ii] = ( p_R - Poo )/( qoo );
    }

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                              *** PRINT_CP ***                              //
//                                                                            //
//============================================================================//
int CAE2DprintCp( vectorField *xx, scalarField *CCp, scalar Cref )
{
    // Variables definition
    FILE *f1;

    // Print Cp results on file ( compatible with gnuplot )
    f1 = fopen("./Log/CCp_up.txt", "w");
    forAll( (*CCp), ii )
    {  
        if ( (*xx)[ii].y() >= 0 ) fprintf(f1, "%lf %lf\n", (*xx)[ii].x()/Cref, -(*CCp)[ii] );
    }
    fclose(f1);   
    f1 = fopen("./Log/CCp_lo.txt", "w");
    forAll( (*CCp), ii )
    {  
        if ( (*xx)[ii].y() <= 0 ) fprintf(f1, "%lf %lf\n", (*xx)[ii].x()/Cref, -(*CCp)[ii] );
    }
    fclose(f1);  
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** COMPUTE_LOADS ***                            //
//                                                                            //
//============================================================================//
int CAE2DcomputeLoads( vectorField *xx, scalarField *SS, vectorField *nn, scalarField *CCp, 
                       vector i_hat, vector j_hat, vector k_hat, scalar Cref, scalar Lref, vector Xref, vector Xhinge, 
                       scalar *C_L, scalar *C_D, scalar *C_MX, scalar *C_H )
{
    
    // Compute Lift, Drag and Moment coefficients
    (*C_L)  = 0.0; 
    (*C_D)  = 0.0;
    (*C_MX) = 0.0;
    (*C_H)  = 0.0;
    forAll( (*SS), ii )
    {
        // REMARK: - nn0, nn vectors point outside the body
        //         - C_MX positive if pitches the airfoil nose up   
        //         - C_H positive if aileron is deflected down increasing camber  
        (*C_L)  = (*C_L)  - (*SS)[ii]*( (*nn)[ii] & j_hat )*(*CCp)[ii];
        (*C_D)  = (*C_D)  - (*SS)[ii]*( (*nn)[ii] & i_hat )*(*CCp)[ii]; 
        (*C_MX) = (*C_MX) + (*SS)[ii]*( ( ( (*xx)[ii] - Xref ) ^ (*nn)[ii]*(*CCp)[ii] ) & k_hat ); 
        if ( (*xx)[ii].x() >= Xhinge.x() )
        {
            (*C_H)  = (*C_H)  + (*SS)[ii]*( ( ( (*xx)[ii] - Xhinge ) ^ (*nn)[ii]*(*CCp)[ii] ) & k_hat ); 
        }         
    }
    (*C_L)  = (*C_L)/Cref;
    (*C_D)  = (*C_D)/Cref;
    (*C_MX) = (*C_MX)/( Cref*Lref );  
    (*C_H)  = (*C_H)/( Cref*Lref );  

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                             *** PRINT_LOADS ***                            //
//                                                                            //
//============================================================================//
int CAE2DprintLoads( scalar t, scalar C_L, scalar C_D, scalar C_MX, scalar C_H )
{
    // Variables definition
    FILE *f1;

    // Print loads on file ( compatible with gnuplot )
    // C_L C_D C_MX
    f1 = fopen("./Log/AerodynamicLoads.txt", "a");
    fprintf(f1, "%g %g %g %g %g\n", t, C_L, C_D, C_MX, C_H );
    fclose(f1); 
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** FORCED_MOTION ***                            //
//                                                                            //
//============================================================================//
int CAE2DforcedMotion( label flag, scalar A0, scalar A1, scalar f, scalar tau, scalar t, 
                       scalar *y, scalar *y_p )
{
    // Variables definition
    scalar pi = 3.14159265358979;
    
    // Oscillatory motion
    if ( flag == 1 )  
    {
        (*y)   = A0 + A1*Foam::sin( 2*pi*f*( t + tau ) );
        (*y_p) =    f*A1*Foam::cos( 2*pi*f*( t + tau ) ); 
    }
    else if ( flag == 2 )  
    {
        (*y)   = A0 + A1*( t + tau );
        (*y_p) = A1; 
    }
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                          *** COMPUTE_NORMALS ***                           //
//                                                                            //
//============================================================================//
int CAE2DcomputeNormals( scalar h, scalar a, scalar d, vectorField *xx, vector Xref, vector Xhinge, vectorField *nn0, vectorField *nn )
{
    // Variables definition
    vector ni;

    // Compute deformed normal vectors
    forAll( (*nn), ii )
    {
        // Pitching movement
        (*nn)[ii].x() =  Foam::cos(a)*(*nn0)[ii].x() + Foam::sin(a)*(*nn0)[ii].y();
        (*nn)[ii].y() = -Foam::sin(a)*(*nn0)[ii].x() + Foam::cos(a)*(*nn0)[ii].y();
        (*nn)[ii].z() = (*nn0)[ii].z(); 
        
        // Aileron deflection correction (in 2D rotation is additive)
        if ( (*xx)[ii].x() >= Xhinge.x() )
        {
            ni = (*nn)[ii];
            (*nn)[ii].x() =  Foam::cos(d)*ni.x() + Foam::sin(d)*ni.y();
            (*nn)[ii].y() = -Foam::sin(d)*ni.x() + Foam::cos(d)*ni.y();
            (*nn)[ii].z() = ni.z(); 
        }
    }
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                       *** COMPUTE_DISPLACEMENTS ***                        //
//                                                                            //
//============================================================================//
int CAE2DcomputeDisplacements( scalar h, scalar a, scalar d, vectorField *xx, vector Xref, vector Xhinge, vectorField *uu )
{
    // Variables definition
    vector ui;

    // Compute displacements
    forAll( (*uu), ii )
    {
        // Pitching and plunging displacement 
        (*uu)[ii].x() =  ( (*xx)[ii].x() - Xref.x() )*Foam::cos( a ) 
                      +  ( (*xx)[ii].y() - Xref.y() )*Foam::sin( a ) 
                      +  Xref.x() - (*xx)[ii].x();
        (*uu)[ii].y() = -( (*xx)[ii].x() - Xref.x() )*Foam::sin( a ) 
                      +  ( (*xx)[ii].y() - Xref.y() )*Foam::cos( a ) 
                      +  Xref.y() - (*xx)[ii].y() + h;  
        (*uu)[ii].z() = 0.0;
        
        // Aileron deflection correction (in 2D rotation is additive) 
        if ( (*xx)[ii].x() >= Xhinge.x() )
        {
            ui = (*uu)[ii];
            (*uu)[ii].x() = ui.x() + ( (*xx)[ii].x() - Xhinge.x() )*Foam::cos( d ) 
                                   + ( (*xx)[ii].y() - Xhinge.y() )*Foam::sin( d ) 
                                   + Xhinge.x() - (*xx)[ii].x();
            (*uu)[ii].y() = ui.y() - ( (*xx)[ii].x() - Xhinge.x() )*Foam::sin( d ) 
                                   + ( (*xx)[ii].y() - Xhinge.y() )*Foam::cos( d ) 
                                   + Xhinge.y() - (*xx)[ii].y();  
            (*uu)[ii].z() = ui.z();                   
        }        
    }

    // Return
    return(0);
}


//============================================================================//
//                                                                            //
//                          *** COMPUTE_VELOCITIES ***                        //
//                                                                            //
//============================================================================//
int CAE2DcomputeVelocities( scalar h_p, scalar a_p, scalar d_p, vectorField *xx, vector Xref, vector Xhinge, vectorField *uu_p )
{
    // Variables definition
    vector i_hat, j_hat, k_hat;
    i_hat = vector(1.0, 0.0, 0.0);
    j_hat = vector(0.0, 1.0, 0.0);
    k_hat = vector(0.0, 0.0, 1.0);

    // Compute velocities
    forAll( (*uu_p), ii )
    {
        // Plunging and pitching velocities
        (*uu_p)[ii] =  h_p*j_hat + ( ( -a_p*k_hat ) ^ ( (*xx)[ii] - Xref ) );
       
        // Aileron deflection correction
        if ( (*xx)[ii].x() >= Xhinge.x() )
        {
            (*uu_p)[ii] = (*uu_p)[ii] + ( ( -d_p*k_hat ) ^ ( (*xx)[ii] - Xhinge ) );
        }
    }

    // Return
    return(0);
}
    
//============================================================================//
//                                                                            //
//                 *** COMPUTE_TRANSPIRATION_VELOCITY ***                     //
//                                                                            //
//============================================================================//
int CAE2DcomputeTranspirationVelocity( vectorField *uu_p, vectorField *nn0, vectorField *nn, volVectorField *U, labelField *id_bodyCell,
                                       scalarField *VVbn )
{
    // Variables definition
    label id_L;
    
    // Loop on control points (boundary faces centres)
    forAll( (*VVbn), ii )
    {
        id_L          = (*id_bodyCell)[ii];      
        // REMARK: nn0, nn vectors point outside the body
        // Warning: in OpenFOAM nn points outside the computational domain
        (*VVbn)[ii]   = -( (*U)[id_L] & ( (*nn)[ii] - (*nn0)[ii] ) ) + ( (*uu_p)[ii] & (*nn)[ii] ); 
    }
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                 *** SET_TRANSPIRATION_VELOCITY_BCS ***                     //
//                                                                            //
//============================================================================//
int CAE2DsetTranspirationVelocityBCs( fvMesh *mesh, label id_bodyPatch, volVectorField *U, scalarField *VVbn, vectorField *nn )
{
    // *** BCs ( Boundary type must not be empty ) ***
    if ( (*U).boundaryField()[id_bodyPatch].size() > 0 )
    {
        // Loop on id_bodyPatch boundary patch faces
        forAll( (*mesh).boundaryMesh()[id_bodyPatch].faceAreas(), ii )
        {
            // REMARK: nn0, nn vectors point outside the body
            // Warning: in OpenFOAM nn points outside the computational domain
            // 1) Classical transpiration boundary condition
            (*U).boundaryFieldRef()[id_bodyPatch][ii].x() = -(*VVbn)[ii];
            // 2) Modified transpiration boundary condition (with deformed normal vector)
            //(*U).boundaryFieldRef()[id_bodyPatch][ii] = (*VVbn)[ii]*(*nn)[ii];
        } 
    }
     
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** SOLVE_ODE_ha ***                             //
//                                                                            //
//============================================================================//
int CAE2DsolveODE_ha( scalar m, scalar d_AECG, scalar I_CG, scalar K_hh, scalar K_aa, // LHS
                      scalar qoo, scalar Cref, scalar Lref, scalar C_L, scalar C_MAE, // RHS
                      scalar dt, scalar *h, scalar *a, scalar *h_p, scalar *a_p  )    // Variables of state
{
    // Variables definition
    scalar M_hh, M_ha, M_aa, detM;
    scalar h_o, a_o, h_p_o, a_p_o;
      
    // Build mass matrix
    M_hh = m;
    M_ha = -m*d_AECG;
    M_aa = I_CG + m*d_AECG*d_AECG;
    detM = M_hh*M_aa - M_ha*M_ha;   
        
    // Solve ODEs for 2 d.o.f elastic airfoil dynamics ( see Bisplinghoff pag. 533 ) (EE)
    h_o   = (*h);
    a_o   = (*a);
    h_p_o = (*h_p);
    a_p_o = (*a_p);
    (*h)   = (*h)   + dt*h_p_o;
    (*a)   = (*a)   + dt*a_p_o;
    (*h_p) = (*h_p) + dt*1.0/detM*(  M_aa*K_hh*h_o - M_ha*K_aa*a_o ) + dt*1.0/detM*qoo*Cref*(  M_aa*C_L - M_ha*Lref*C_MAE );
    (*a_p) = (*a_p) + dt*1.0/detM*( -M_ha*K_hh*h_o + M_hh*K_aa*a_o ) + dt*1.0/detM*qoo*Cref*( -M_ha*C_L + M_ha*Lref*C_MAE );

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** SOLVE_ODE_d ***                             //
//                                                                            //
//============================================================================//
int CAE2DsolveODE_d( scalar I_dd, scalar C_dd, scalar K_dd, scalar qoo, scalar Cref, scalar Lref, scalar C_H, scalar dt, scalar *d, scalar *d_p ) 
{
    // Variables definition
    scalar d_o, d_p_o; 
        
    // Solve ODE for free moving elastic aileron
    d_o    = (*d);
    d_p_o  = (*d_p);
    (*d)   = (*d)   + dt*d_p_o;
    (*d_p) = (*d_p) + dt*( qoo*Cref*Lref*C_H/I_dd - K_dd/I_dd*d_o - C_dd/I_dd*d_p_o );

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** PRINT_MOVEMENT ***                           //
//                                                                            //
//============================================================================//
int CAE2DprintMovement( scalar t, scalar h, scalar alpha, scalar delta, scalar h_p, scalar alpha_p, scalar delta_p )
{
    // Variables definition
    scalar pi = 3.14159265358979;
    FILE *f1;

    // Print loads on file ( compatible with geplot.out )
    // h, a, d, h_p, a_p, d_p
    f1 = fopen("./Log/Displacements.txt", "a");
    fprintf(f1, "%g %g %g %g %g %g %g\n", t, h, alpha*180.0/pi, delta*180.0/pi, h_p, alpha_p*180.0/pi, delta_p*180.0/pi );
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                       *** UPDATECELLDISPLACEMENT ***                       //
//                                                                            //
//============================================================================//
int CAE2DupdateCellDisplacement( label id_bodyPatch, volVectorField *cellDisplacement, vectorField *uu )
{
    // Update cellDisplacement for later plot deformed mesh with showDisplacement and paraFoam
    forAll( (*uu), ii )
    {
        (*cellDisplacement).boundaryFieldRef()[id_bodyPatch][ii] = (*uu)[ii];  
    }

    // Return
    return(0);
}
