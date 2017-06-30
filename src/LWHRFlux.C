//============================================================================//
//                                                                            //
//                              *** LWDU ***                                  //
//                                                                            //
//============================================================================//
scalar lwDU( scalar dt, scalar dx, scalar lambda, scalar alpha_, scalar alpha_L_, scalar alpha_R_, word fluxLim )
{
    // Variables definition
    scalar Co, a, b, theta, min1, min2, max1, max2, psi, DU, eps = 1e-10;

    // Lax-Wendroff increment with VanLeer limiter 
    Co  = lambda*dt/( dx + eps );
    a   = alpha_;
    b   = 0.5*( alpha_L_ + alpha_R_ ) + 0.5*Foam::sign( Co )*( alpha_L_ - alpha_R_ );
    //IF{ fluxLim VL
    if ( fluxLim == "VL" ) // VanLeer
    {
        psi = ( a*mag(b) + mag(a)*b )/( mag(a) + mag(b) + eps );
    }
    //FI} fluxLim VL
    //IF{ fluxLim MinMod
    else if ( fluxLim == "MinMod" ) // MinMod
    {
        theta = ( mag(a) < mag(b) );
        psi = ( 0.5*Foam::sign(a) + 0.5*Foam::sign(b) )*( mag(a)*theta + mag(b)*( 1 - theta ) );
    }
    //FI} fluxLim MinMod
    //IF{ fluxLim SB
    else if ( fluxLim == "SB" ) // Superbee ( to be optimized )
    {
        min1 = 2*mag(b);
        if ( mag(a) < min1 ) min1 = mag(a);
        min2 = 2*mag(a);
        if ( mag(b) < min2 ) min2 = mag(b);
        theta = ( min1 > min2 );
        psi = ( 0.5*Foam::sign(a) + 0.5*Foam::sign(b) )*( min1*theta + min2*( 1 - theta ) );
    }
    //FI} fluxLim SB
    //IF{ fluxLim MC
    else //if ( fluxLim == "MC" ) // Monotonized-Central ( to be optimized )
    {
        min1 = 0.5*( a + b );
        if ( 2*a < min1 ) min1 = 2*a;
        if ( 2*b < min1 ) min1 = 2*b;
        max2 = 0.5*( a + b );
        if ( 2*a > max2 ) max2 = 2*a;
        if ( 2*b > max2 ) max2 = 2*b;
        max1 = min1;
        if ( max1 < 0 ) max1 = 0;
        min2 = max2;
        if ( min2 > 0 ) min2 = 0;
        psi = max1 + min2;
    }
    //FI} fluxLim MC
    DU  = lambda*psi*0.5*( Foam::sign(Co) - Co );
    
    // Output
    return(DU);
}

//============================================================================//
//                                                                            //
//                            *** LWHR_FLUX ***                               //
//                                                                            //
//============================================================================//
int lwhr_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, 
                              scalar rho_LL, scalar rho_RR, vector mm_LL, vector mm_RR, scalar et_tilde_LL, scalar et_tilde_RR, 
                              vector n_i, scalar dt, scalar dx,
                              scalar *Frho, vector *Fm, scalar *Fet_tilde, word fluxLim ) 
{  
    
    // Variables definition
    scalar ux_L,  uy_L,  uz_L,  ux_R,  uy_R,  uz_R;
    scalar ux_LL, uy_LL, uz_LL, ux_RR, uy_RR, uz_RR;
    scalar  u_L,  v_L,   w_L,   u_R,   v_R,   w_R;
    scalar u_LL,  v_LL,  w_LL,  u_RR,  v_RR,  w_RR;
    scalar e_tilde_L, e_tilde_R, ht_L, ht_R;
    scalar magU_L, magU_R, p_L, p_R, c_L, c_R;
    scalar nx, ny, nz, norm, eps;
    
    scalar u_hat, v_hat, w_hat, c_hat, magU_hat, ht_hat;
    
    scalar lambda1, lambda2, lambda3;
    
    scalar R11, R12, R13, R14, R15;
    scalar R21, R22, R23, R24, R25; 
    scalar R31, R32, R33, R34, R35; 
    scalar R41, R42, R43, R44, R45; 
    scalar R51, R52, R53, R54, R55;
    
    scalar L11, L12, L13, L14, L15; 
    scalar L21, L22, L23, L24, L25; 
    scalar L31, L32, L33, L34, L35; 
    scalar L41, L42, L43, L44, L45; 
    scalar L51, L52, L53, L54, L55;
    
    scalar DU1, DU2, DU3, DU4, DU5;

    scalar Flwrho, Flwet_tilde;
    vector Flwm;
    
    scalar alpha_rho, alpha_et_tilde, alpha_L_rho, alpha_L_et_tilde, alpha_R_rho, alpha_R_et_tilde;
    vector alpha_m, alpha_L_m, alpha_R_m;
    
    // Change frame of reference  
    eps = 1e-10;
    ux_L = mm_L.x()/rho_L;
    uy_L = mm_L.y()/rho_L;
    uz_L = mm_L.z()/rho_L;
    ux_R = mm_R.x()/rho_R;
    uy_R = mm_R.y()/rho_R;
    uz_R = mm_R.z()/rho_R;
    ux_LL = mm_LL.x()/rho_LL;
    uy_LL = mm_LL.y()/rho_LL;
    uz_LL = mm_LL.z()/rho_LL;
    ux_RR = mm_RR.x()/rho_RR;
    uy_RR = mm_RR.y()/rho_RR;
    uz_RR = mm_RR.z()/rho_RR;
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
    u_LL = nx*ux_LL + ny*uy_LL + nz*uz_LL;
    u_RR = nx*ux_RR + ny*uy_RR + nz*uz_RR; 
    if ( mag(nz) < 1.0 )
    {
        v_L  = -ny/Foam::sqrt( 1 - sqr(nz) )*ux_L + nx/Foam::sqrt( 1 - sqr(nz) )*uy_L;
        v_R  = -ny/Foam::sqrt( 1 - sqr(nz) )*ux_R + nx/Foam::sqrt( 1 - sqr(nz) )*uy_R;
        v_LL = -ny/Foam::sqrt( 1 - sqr(nz) )*ux_LL + nx/Foam::sqrt( 1 - sqr(nz) )*uy_LL;
        v_RR = -ny/Foam::sqrt( 1 - sqr(nz) )*ux_RR + nx/Foam::sqrt( 1 - sqr(nz) )*uy_RR;
        w_L  = -nx*nz/Foam::sqrt( 1 - sqr(nz) )*ux_L - ny*nz/Foam::sqrt( 1 - sqr(nz) )*uy_LL + Foam::sqrt( 1 - sqr(nz) )*uz_L;
        w_R  = -nx*nz/Foam::sqrt( 1 - sqr(nz) )*ux_R - ny*nz/Foam::sqrt( 1 - sqr(nz) )*uy_RR + Foam::sqrt( 1 - sqr(nz) )*uz_R;
        w_LL = -nx*nz/Foam::sqrt( 1 - sqr(nz) )*ux_LL - ny*nz/Foam::sqrt( 1 - sqr(nz) )*uy_LL + Foam::sqrt( 1 - sqr(nz) )*uz_LL;
        w_RR = -nx*nz/Foam::sqrt( 1 - sqr(nz) )*ux_RR - ny*nz/Foam::sqrt( 1 - sqr(nz) )*uy_RR + Foam::sqrt( 1 - sqr(nz) )*uz_RR;
    }
    else
    {
        v_L  = nz*ux_L;
        v_R  = nz*ux_R;
        v_LL = nz*ux_LL;
        v_RR = nz*ux_RR;
        w_L  = uy_L;
        w_R  = uy_R;
        w_LL = uy_LL;
        w_RR = uy_RR;
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
    
    // Roe's average   
    u_hat    = ( u_L*Foam::sqrt(rho_L) + u_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    v_hat    = ( v_L*Foam::sqrt(rho_L) + v_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    w_hat    = ( w_L*Foam::sqrt(rho_L) + w_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    ht_hat   = ( ht_L*Foam::sqrt(rho_L) + ht_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    magU_hat = Foam::sqrt( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) );
    c_hat    = Foam::sqrt( ( gamma_ - 1 )*( ht_hat - 0.5*Foam::sqr( magU_hat ) ) );
    
    // Right eigenvectors
    R11 = 1.0;
    R12 = 1.0;
    R13 = 1.0;
    R14 = 1.0;
    R15 = 1.0;
    R21 = u_hat - c_hat;
    R22 = u_hat;
    R23 = u_hat;
    R24 = u_hat;
    R25 = u_hat + c_hat;
    R31 = v_hat;
    R32 = v_hat;
    R33 = v_hat - c_hat;
    R34 = v_hat;
    R35 = v_hat;
    R41 = w_hat;
    R42 = w_hat;
    R43 = w_hat;
    R44 = w_hat - c_hat;
    R45 = w_hat;
    R51 = 0.5*sqr(magU_hat) + sqr(c_hat)/( gamma_ - 1 ) - u_hat*c_hat;
    R52 = 0.5*sqr(magU_hat);
    R53 = 0.5*sqr(magU_hat) - v_hat*c_hat;
    R54 = 0.5*sqr(magU_hat) - w_hat*c_hat;
    R55 = 0.5*sqr(magU_hat) + sqr(c_hat)/( gamma_ - 1 ) + u_hat*c_hat;
    
    // Left eigenvectors
    L11 = 0.5*( 0.5*( gamma_ - 1 )*sqr(magU_hat) + u_hat*c_hat );
    L12 = -0.5*( c_hat + ( gamma_ - 1 )*u_hat );
    L13 = -0.5*( gamma_ - 1 )*v_hat;
    L14 = -0.5*( gamma_ - 1 )*w_hat;
    L15 = 0.5*( gamma_ - 1 );
    L21 = sqr(c_hat) - c_hat*( v_hat + w_hat ) - 0.5*( gamma_ - 1 )*sqr(magU_hat);
    L22 = ( gamma_ - 1 )*u_hat;
    L23 = c_hat + ( gamma_ - 1 )*v_hat;
    L24 = c_hat + ( gamma_ - 1 )*w_hat;
    L25 = 1 - gamma_;
    L31 = v_hat*c_hat;
    L32 = 0;
    L33 = -c_hat;
    L34 = 0;
    L35 = 0;
    L41 = w_hat*c_hat;
    L42 = 0;
    L43 = 0;
    L44 = -c_hat;
    L45 = 0;
    L51 = 0.5*( 0.5*( gamma_ - 1 )*sqr(magU_hat) - u_hat*c_hat );
    L52 = 0.5*( c_hat - ( gamma_ - 1 )*u_hat );
    L53 = -0.5*( gamma_ - 1 )*v_hat;
    L54 = -0.5*( gamma_ - 1 )*w_hat;
    L55 = 0.5*( gamma_ - 1 );
       
    //==========================================================================    
    // *** HIGH ORDER LW FLUXES *** 
    //==========================================================================   
    
    // (j + 1/2) cell
    DU1 = ( rho_R      - rho_L      )/sqr(c_hat);
    DU2 = ( rho_R*u_R  - rho_L*u_L  )/sqr(c_hat);
    DU3 = ( rho_R*v_R  - rho_L*v_L  )/sqr(c_hat);
    DU4 = ( rho_R*w_R  - rho_L*w_L  )/sqr(c_hat);
    DU5 = ( et_tilde_R - et_tilde_L )/sqr(c_hat);
    alpha_rho      = L11*DU1 + L12*DU2 + L13*DU3 + L14*DU4 + L15*DU5;
    alpha_m.x()    = L21*DU1 + L22*DU2 + L23*DU3 + L24*DU4 + L25*DU5;
    alpha_m.y()    = L31*DU1 + L32*DU2 + L33*DU3 + L34*DU4 + L35*DU5;
    alpha_m.z()    = L41*DU1 + L42*DU2 + L43*DU3 + L44*DU4 + L45*DU5;
    alpha_et_tilde = L51*DU1 + L52*DU2 + L53*DU3 + L54*DU4 + L55*DU5;
    
    // (j + 3/2) cell
    DU1 = ( rho_RR      - rho_R      )/sqr(c_hat);
    DU2 = ( rho_RR*u_RR - rho_R*u_R  )/sqr(c_hat);
    DU3 = ( rho_RR*v_RR - rho_R*v_R  )/sqr(c_hat);
    DU4 = ( rho_RR*w_RR - rho_R*w_R  )/sqr(c_hat);
    DU5 = ( et_tilde_RR - et_tilde_R )/sqr(c_hat);
    alpha_R_rho      = L11*DU1 + L12*DU2 + L13*DU3 + L14*DU4 + L15*DU5;
    alpha_R_m.x()    = L21*DU1 + L22*DU2 + L23*DU3 + L24*DU4 + L25*DU5;
    alpha_R_m.y()    = L31*DU1 + L32*DU2 + L33*DU3 + L34*DU4 + L35*DU5;
    alpha_R_m.z()    = L41*DU1 + L42*DU2 + L43*DU3 + L44*DU4 + L45*DU5;
    alpha_R_et_tilde = L51*DU1 + L52*DU2 + L53*DU3 + L54*DU4 + L55*DU5;
    
    // (j - 1/2) cell
    DU1 = ( rho_L      - rho_LL      )/sqr(c_hat);
    DU2 = ( rho_L*u_L  - rho_LL*u_LL )/sqr(c_hat);
    DU3 = ( rho_L*v_L  - rho_LL*v_LL )/sqr(c_hat);
    DU4 = ( rho_L*w_L  - rho_LL*w_LL )/sqr(c_hat);
    DU5 = ( et_tilde_L - et_tilde_LL )/sqr(c_hat);
    alpha_L_rho      = L11*DU1 + L12*DU2 + L13*DU3 + L14*DU4 + L15*DU5;
    alpha_L_m.x()    = L21*DU1 + L22*DU2 + L23*DU3 + L24*DU4 + L25*DU5;
    alpha_L_m.y()    = L31*DU1 + L32*DU2 + L33*DU3 + L34*DU4 + L35*DU5;
    alpha_L_m.z()    = L41*DU1 + L42*DU2 + L43*DU3 + L44*DU4 + L45*DU5;
    alpha_L_et_tilde = L51*DU1 + L52*DU2 + L53*DU3 + L54*DU4 + L55*DU5;
    
    // Eigenvalues 
    lambda1 = u_hat - c_hat;
    lambda2 = u_hat;
    lambda3 = u_hat + c_hat;
    
    // Lax-Wendroff increment with limiter
    DU1 = lwDU( dt, dx, lambda1, alpha_rho, alpha_L_rho, alpha_R_rho, fluxLim );
    DU2 = lwDU( dt, dx, lambda2, alpha_m.x(), alpha_L_m.x(), alpha_R_m.x(), fluxLim );
    DU3 = lwDU( dt, dx, lambda2, alpha_m.y(), alpha_L_m.y(), alpha_R_m.y(), fluxLim );
    DU4 = lwDU( dt, dx, lambda2, alpha_m.z(), alpha_L_m.z(), alpha_R_m.z(), fluxLim );
    DU5 = lwDU( dt, dx, lambda3, alpha_et_tilde, alpha_L_et_tilde, alpha_R_et_tilde, fluxLim );
    
    // Update fluxes
    Flwrho      = ( R11*DU1 + R12*DU2 + R13*DU3 + R14*DU4 + R15*DU5 );
    Flwm.x()    = ( R21*DU1 + R22*DU2 + R23*DU3 + R24*DU4 + R25*DU5 );
    Flwm.y()    = ( R31*DU1 + R32*DU2 + R33*DU3 + R34*DU4 + R35*DU5 );
    Flwm.z()    = ( R41*DU1 + R42*DU2 + R43*DU3 + R44*DU4 + R45*DU5 );
    Flwet_tilde = ( R51*DU1 + R52*DU2 + R53*DU3 + R54*DU4 + R55*DU5 );
    
    // Back in the Cartesian frame of reference   
    *Frho      = *Frho      + Flwrho;
    if ( mag(nz) < 1.0 )
    {
    (*Fm).x()  = (*Fm).x()  + Flwm.x()*nx - Flwm.y()*ny/( Foam::sqrt( 1 - sqr(nz) ) ) - Flwm.z()*nx*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).y()  = (*Fm).y()  + Flwm.x()*ny + Flwm.y()*nx/( Foam::sqrt( 1 - sqr(nz) ) ) - Flwm.z()*ny*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).z()  = (*Fm).z()  + Flwm.x()*nz + Flwm.z()*( Foam::sqrt( 1 - sqr(nz) ) );
    }
    else
    {
    (*Fm).x()  = (*Fm).x();
    (*Fm).y()  = (*Fm).y();
    (*Fm).z()  = (*Fm).z()  + Flwm.x()*nz;
    }
    *Fet_tilde = *Fet_tilde + Flwet_tilde;
    
    return(0);
    
}    
