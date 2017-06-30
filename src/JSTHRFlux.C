//============================================================================//
//                                                                            //
//                              *** JSTHR_FLUX ***                            //
//                                                                            //
//============================================================================//
int jsthr_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, 
                               scalar rho_LL, scalar rho_RR, vector mm_LL, vector mm_RR, scalar et_tilde_LL, scalar et_tilde_RR, 
                               vector n_i, scalar *Frho, vector *Fm, scalar *Fet_tilde ) 
{  
    
    // Variables definition
    scalar ux_L, uy_L, uz_L, ux_R, uy_R, uz_R, u_L, u_R;
    scalar e_tilde_L, e_tilde_R, magU_L, magU_R, p_L, p_R, p_LL, p_RR, c_L, c_R;
    scalar nx, ny, nz, norm, eps = 1e-10;
    scalar lambda1_L, lambda2_L, lambda3_L;
    scalar lambda1_R, lambda2_R, lambda3_R;
    scalar lambda_L, lambda_R, lambda;
    scalar nu_L, nu_R, k2, k4, e2, e4;
    scalar Fjstrho, Fjstet_tilde;
    vector Fjstm;

    // Change frame of reference
    ux_L = mm_L.x()/rho_L;
    uy_L = mm_L.y()/rho_L;
    uz_L = mm_L.z()/rho_L;
    ux_R = mm_R.x()/rho_R;
    uy_R = mm_R.y()/rho_R;
    uz_R = mm_R.z()/rho_R;
    nx = n_i.x();
    ny = n_i.y();
    nz = n_i.z();  
    if ( mag(nx) < eps ) nx = 0.0;
    if ( mag(ny) < eps ) ny = 0.0;
    if ( mag(nz) < eps ) nz = 0.0;
    norm = Foam::sqrt( sqr(nx) + sqr(ny) + sqr(nz) );
    nx   = nx/norm;
    ny   = ny/norm;
    nz   = nz/norm;    
    u_L  = nx*ux_L + ny*uy_L + nz*uz_L;
    u_R  = nx*ux_R + ny*uy_R + nz*uz_R;
    magU_L = Foam::sqrt( sqr(ux_L) + sqr(uy_L) + sqr(uz_L) );
    magU_R = Foam::sqrt( sqr(ux_R) + sqr(uy_R) + sqr(uz_R) );
    e_tilde_L = et_tilde_L - 0.5*rho_L*sqr(magU_L);
    e_tilde_R = et_tilde_R - 0.5*rho_R*sqr(magU_R);
    p_L = ( gamma_ - 1 )*e_tilde_L;
    p_R = ( gamma_ - 1 )*e_tilde_R;  
    c_L = Foam::sqrt( gamma_*p_L/rho_L );
    c_R = Foam::sqrt( gamma_*p_R/rho_R );
    
    // Eigenvalues
    lambda1_L = u_L - c_L;
    lambda2_L = u_L;
    lambda3_L = u_L + c_L;
    lambda1_R = u_R - c_R;
    lambda2_R = u_R;
    lambda3_R = u_R + c_R;
    lambda_L = lambda1_L;
    if ( lambda2_L > lambda_L ) lambda_L = lambda2_L;
    if ( lambda3_L > lambda_L ) lambda_L = lambda3_L;
    lambda_L = mag( lambda_L );
    lambda_R = lambda1_R;
    if ( lambda2_R > lambda_R ) lambda_R = lambda2_R;
    if ( lambda3_R > lambda_R ) lambda_R = lambda3_R;
    lambda_R = mag( lambda_R );
    
    // Initialization
    Fjstrho      = 0.0;
    Fjstm.x()    = 0.0;
    Fjstm.y()    = 0.0;
    Fjstm.z()    = 0.0; 
    Fjstet_tilde = 0.0;
    
    // Mean value
    lambda = 0.5*( lambda_L + lambda_R );
    //lambda = 0.5*( u_L + u_R ) + 0.5*( c_L + c_R );  
    
    // Parameters evaluation
    p_LL = ( gamma_ - 1 )*( et_tilde_LL - 0.5*magSqr(mm_LL)/rho_LL );
    p_RR = ( gamma_ - 1 )*( et_tilde_RR - 0.5*magSqr(mm_RR)/rho_RR ); 
    nu_L = mag( p_R - 2*p_L + p_LL )/( p_R + 2*p_L + p_LL ); 
    nu_R = mag( p_RR - 2*p_R + p_L )/( p_RR + 2*p_R + p_L );
    k2 = 1.0/4.0;
    k4 = 1.0/256.0;
    e2 = k2*nu_L; if ( k2*nu_R > e2 ) e2 = k2*nu_R; 
    e4 = 0.0; if ( ( k4 - e2 ) > e4 ) e4 = k4 - e2; 
    
    // Update JST fluxes
    Fjstrho      = -lambda*( e2*( rho_R      - rho_L      ) - e4*( rho_RR      - 3*rho_R      + 3*rho_L      - rho_LL      ) );
    Fjstm.x()    = -lambda*( e2*( mm_R.x()   - mm_L.x()   ) - e4*( mm_RR.x()   - 3*mm_R.x()   + 3*mm_L.x()   - mm_LL.x()   ) );
    Fjstm.y()    = -lambda*( e2*( mm_R.y()   - mm_L.y()   ) - e4*( mm_RR.y()   - 3*mm_R.y()   + 3*mm_L.y()   - mm_LL.y()   ) );
    Fjstm.z()    = -lambda*( e2*( mm_R.z()   - mm_L.z()   ) - e4*( mm_RR.z()   - 3*mm_R.z()   + 3*mm_L.z()   - mm_LL.z()   ) );
    Fjstet_tilde = -lambda*( e2*( et_tilde_R - et_tilde_L ) - e4*( et_tilde_RR - 3*et_tilde_R + 3*et_tilde_L - et_tilde_LL ) );    
  
    // Update global fluxes
    *Frho      = *Frho      + Fjstrho;
    (*Fm).x()  = (*Fm).x()  + Fjstm.x();
    (*Fm).y()  = (*Fm).y()  + Fjstm.y();
    (*Fm).z()  = (*Fm).z()  + Fjstm.z();
    *Fet_tilde = *Fet_tilde + Fjstet_tilde;
    
    return(0);
    
}    
