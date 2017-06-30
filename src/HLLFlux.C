//============================================================================//
//                                                                            //
//                              *** HLL_FLUX ***                              //
//                                                                            //
//============================================================================//
int hll_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, vector n_i, 
              scalar *Frho, vector *Fm, scalar *Fet_tilde, word method )
{  
    
    // Variables definition
    scalar nx, ny, nz, norm, eps;                                  // Reference frame   
    scalar ux_L, uy_L, uz_L, ux_R, uy_R, uz_R;                     // 
    scalar  u_L,  v_L,  w_L,  u_R,  v_R,  w_R, magU_L, magU_R;     // L and R variables
    scalar  p_L, p_R, c_L, c_R, e_tilde_L, e_tilde_R, ht_L, ht_R;  //    
    scalar lambda1_L, lambda2_L, lambda3_L, S_L;                   // S_L and S_R variables ( depending on method )
    scalar lambda1_R, lambda2_R, lambda3_R, S_R;                   //
    scalar Rrho, u_hat, v_hat, w_hat, c_hat, ht_hat;               // Roe average ( Einfeldt method )
    scalar lambda1_hat, lambda2_hat, lambda3_hat;                  //
    scalar Fhllrho, Fhllet_tilde;                                  // Fluxes
    vector Fhllm;                                                  //
       
    // Change frame of reference 
    eps = 1e-10;
    nx  = n_i.x();
    ny  = n_i.y();
    nz  = n_i.z();  
    if ( mag(nx) < eps ) nx = 0.0;
    if ( mag(ny) < eps ) ny = 0.0;
    if ( mag(nz) < eps ) nz = 0.0;
    norm = Foam::sqrt( sqr(nx) + sqr(ny) + sqr(nz) );
    nx   = nx/norm;
    ny   = ny/norm;
    nz   = nz/norm;
    
    // L and R variables
    ux_L = mm_L.x()/rho_L;
    uy_L = mm_L.y()/rho_L;
    uz_L = mm_L.z()/rho_L;
    ux_R = mm_R.x()/rho_R;
    uy_R = mm_R.y()/rho_R;
    uz_R = mm_R.z()/rho_R;    
    u_L = nx*ux_L + ny*uy_L + nz*uz_L;
    u_R = nx*ux_R + ny*uy_R + nz*uz_R;
    if ( mag(nz) < 1.0 )
    {
        v_L = -ny/Foam::sqrt( 1 - sqr(nz) )*ux_L + nx/Foam::sqrt( 1 - sqr(nz) )*uy_L;
        v_R = -ny/Foam::sqrt( 1 - sqr(nz) )*ux_R + nx/Foam::sqrt( 1 - sqr(nz) )*uy_R;
        w_L = -nx*nz/Foam::sqrt( 1 - sqr(nz) )*ux_L - ny*nz/Foam::sqrt( 1 - sqr(nz) )*uy_L + Foam::sqrt( 1 - sqr(nz) )*uz_L;
        w_R = -nx*nz/Foam::sqrt( 1 - sqr(nz) )*ux_R - ny*nz/Foam::sqrt( 1 - sqr(nz) )*uy_R + Foam::sqrt( 1 - sqr(nz) )*uz_R;
    }
    else
    {
        v_L = ux_L*nz;
        v_R = ux_R*nz;
        w_L = uy_L;
        w_R = uy_R;
    }
    magU_L = Foam::sqrt( sqr(ux_L) + sqr(uy_L) + sqr(uz_L) );
    magU_R = Foam::sqrt( sqr(ux_R) + sqr(uy_R) + sqr(uz_R) );
    e_tilde_L = et_tilde_L - 0.5*rho_L*sqr(magU_L);
    e_tilde_R = et_tilde_R - 0.5*rho_R*sqr(magU_R);
    p_L = ( gamma_ - 1 )*e_tilde_L;
    p_R = ( gamma_ - 1 )*e_tilde_R;
    ht_L = ( et_tilde_L + p_L )/rho_L;
    ht_R = ( et_tilde_R + p_R )/rho_R;
    c_L = Foam::sqrt( gamma_*p_L/rho_L );
    c_R = Foam::sqrt( gamma_*p_R/rho_R );
    
    // Eigenvalues
    lambda1_L = u_L - c_L;
    lambda2_L = u_L;
    lambda3_L = u_L + c_L;
    lambda1_R = u_R - c_R;
    lambda2_R = u_R;
    lambda3_R = u_R + c_R; 

    // S_L and S_R variables ( Davies # 1 )
    //IF{ entropyFix D1
    if ( method == "D1" )
    {
        S_L = lambda1_L;
        S_R = lambda3_R;
    }
    //FI} entropyFix D1
    // S_L and S_R variables ( Davies # 2 )
    //IF{ entropyFix D2
    else if ( method == "D2" )
    {
        S_L = lambda1_L;
        if ( lambda1_R < S_L ) S_L = lambda1_R;
        S_R = lambda3_R;
        if ( lambda3_L > S_R ) S_R = lambda3_L;
    }
    //FI} entropyFix D2
    // S_L and S_R variables ( Einfeldt-Roe )
    //IF{ entropyFix ER
    else if ( method == "ER" )
    {
        // Roe average state
        Rrho   = Foam::sqrt( rho_R/rho_L );
        u_hat  = ( u_L  + u_R*Rrho  )/( 1 + Rrho );
        v_hat  = ( v_L  + v_R*Rrho  )/( 1 + Rrho );
        w_hat  = ( w_L  + w_R*Rrho  )/( 1 + Rrho );
        ht_hat = ( ht_L + ht_R*Rrho )/( 1 + Rrho );
        c_hat  = Foam::sqrt( ( gamma_ - 1 )*( ht_hat - 0.5*( Foam::pow( u_hat, 2.0 ) + Foam::pow( v_hat, 2.0 ) + Foam::pow( w_hat, 2.0 ) ) ) );
        S_L    = lambda1_L;
        S_R    = lambda3_R;
        lambda1_hat = u_hat - c_hat;
        lambda2_hat = u_hat;
        lambda3_hat = u_hat + c_hat;        
        if ( lambda1_hat < S_L ) S_L = lambda1_hat;
        if ( lambda3_hat > S_R ) S_R = lambda3_hat; 
    }
    //FI} entropyFix ER
    // DEFAULT S_L and S_R variables ( LeVeque ) 
    //IF{ entropyFix LV
    else // if ( method == "LV" )
    {
        S_L = lambda1_L;
        if ( lambda1_R < S_L ) S_L = lambda1_R;
        if ( lambda2_L < S_L ) S_L = lambda2_L;
        if ( lambda2_R < S_L ) S_L = lambda2_R;
        if ( lambda3_L < S_L ) S_L = lambda3_L;
        if ( lambda3_R < S_L ) S_L = lambda3_R;
        S_R = lambda1_L;
        if ( lambda1_R > S_R ) S_R = lambda1_R;
        if ( lambda2_L > S_R ) S_R = lambda2_L;
        if ( lambda2_R > S_R ) S_R = lambda2_R;
        if ( lambda3_L > S_R ) S_R = lambda3_L;
        if ( lambda3_R > S_R ) S_R = lambda3_R;
    }
    //FI} entropyFix LV
    
    // HLL Fluxes
    if ( S_L > 0 )
    {
        Fhllrho      = rho_L*u_L;
        Fhllm.x()    = rho_L*u_L*u_L + p_L;
        Fhllm.y()    = rho_L*u_L*v_L;
        Fhllm.z()    = rho_L*u_L*w_L;
        Fhllet_tilde = u_L*( et_tilde_L + p_L );
    }
    else if ( S_R < 0 )
    {
        Fhllrho      = rho_R*u_R;
        Fhllm.x()    = rho_R*u_R*u_R + p_R;
        Fhllm.y()    = rho_R*u_R*v_R;
        Fhllm.z()    = rho_R*u_R*w_R;
        Fhllet_tilde = u_R*( et_tilde_R + p_R );
    }
    else
    {
        Fhllrho      = ( S_R*rho_L*u_L - S_L*rho_R*u_R + S_L*S_R*( rho_R - rho_L ) )/( S_R - S_L + eps );
        Fhllm.x()    = ( S_R*( rho_L*u_L*u_L + p_L ) - S_L*( rho_R*u_R*u_R + p_R ) + S_L*S_R*( rho_R*u_R - rho_L*u_R ) )/( S_R - S_L + eps );
        Fhllm.y()    = ( S_R*( rho_L*u_L*v_L ) - S_L*( rho_R*u_R*v_R ) + S_L*S_R*( rho_R*v_R - rho_L*v_R ) )/( S_R - S_L + eps );
        Fhllm.z()    = ( S_R*( rho_L*u_L*w_L ) - S_L*( rho_R*u_R*w_R ) + S_L*S_R*( rho_R*w_R - rho_L*w_R ) )/( S_R - S_L + eps );
        Fhllet_tilde = ( S_R*( u_L*( et_tilde_L + p_L ) ) - S_L*( u_R*( et_tilde_R + p_R ) ) + S_L*S_R*( et_tilde_R - et_tilde_L ) )/( S_R - S_L + eps );
    }

    // Back in the Cartesian reference frame 
    *Frho      = Fhllrho;
    if ( mag(nz) < 1.0 )
    {
    (*Fm).x()  = Fhllm.x()*nx - Fhllm.y()*ny/( Foam::sqrt( 1 - sqr(nz) ) ) - Fhllm.z()*nx*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).y()  = Fhllm.x()*ny + Fhllm.y()*nx/( Foam::sqrt( 1 - sqr(nz) ) ) - Fhllm.z()*ny*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).z()  = Fhllm.x()*nz + Fhllm.z()*( Foam::sqrt( 1 - sqr(nz) ) );
    }
    else
    {
    (*Fm).x()  = 0.0;
    (*Fm).y()  = 0.0;
    (*Fm).z()  = Fhllm.x()*nz;
    }
    *Fet_tilde = Fhllet_tilde;
    return(0);
    
}    
