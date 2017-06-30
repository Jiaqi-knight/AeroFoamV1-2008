//============================================================================//
//                                                                            //
//                            *** ROE_FLUX ***                                //
//                                                                            //
//============================================================================//
int roe_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, vector n_i, 
              scalar *Frho, vector *Fm, scalar *Fet_tilde, word entropyFix )
{  
    
    // Variables definition
    scalar ux_L, uy_L, uz_L, ux_R, uy_R, uz_R;
    scalar  u_L,  v_L,  w_L,  u_R,  v_R,  w_R;
    scalar e_tilde_L, e_tilde_R, ht_L, ht_R;
    scalar magU_L, magU_R, p_L, p_R, c_L, c_R;
    scalar nx, ny, nz, norm;
    
    scalar u_hat, v_hat, w_hat, c_hat, magU_hat, ht_hat;
    
    scalar lambda1_L,  lambda2_L,  lambda3_L;
    scalar lambda1_R,  lambda2_R,  lambda3_R;
    scalar lambda1_EF, lambda2_EF, lambda3_EF;
    scalar delta, eps;
    
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

    scalar Fsyrho, Fsyet_tilde, Fuprho, Fupet_tilde;
    vector Fsym, Fupm;
    
    eps = 1e-10;
    
    //==========================================================================   
    // *** SYMMETRIC FLUXES *** 
    //==========================================================================    
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
        v_L = nz*ux_L;
        v_R = nz*ux_R;
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
        
    // Assembly symmetric fluxes
    Fsyrho      = 0.5*( rho_L*u_L + rho_R*u_R );
    Fsym.x()    = 0.5*( rho_L*u_L*( ux_L ) + p_L*nx + rho_R*u_R*( ux_R ) + p_R*nx );
    Fsym.y()    = 0.5*( rho_L*u_L*( uy_L ) + p_L*ny + rho_R*u_R*( uy_R ) + p_R*ny );
    Fsym.z()    = 0.5*( rho_L*u_L*( uz_L ) + p_L*nz + rho_R*u_R*( uz_R ) + p_R*nz );
    Fsyet_tilde = 0.5*( u_L*( et_tilde_L + p_L ) + u_R*( et_tilde_R + p_R ) );
    
    //==========================================================================    
    // *** UPWIND FLUXES *** 
    //==========================================================================    
    lambda1_L = u_L - c_L;
    lambda2_L = u_L;
    lambda3_L = u_L + c_L;
    lambda1_R = u_R - c_R;
    lambda2_R = u_R;
    lambda3_R = u_R + c_R;
    
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
    
    // Jump variables
    DU1 = ( rho_R      - rho_L      )/sqr(c_hat);
    DU2 = ( rho_R*u_R  - rho_L*u_L  )/sqr(c_hat);
    DU3 = ( rho_R*v_R  - rho_L*v_L  )/sqr(c_hat);
    DU4 = ( rho_R*w_R  - rho_L*w_L  )/sqr(c_hat);
    DU5 = ( et_tilde_R - et_tilde_L )/sqr(c_hat);

    lambda1_EF = ( u_hat - c_hat );
    lambda2_EF = ( u_hat );
    lambda3_EF = ( u_hat + c_hat );
    //IF{ entropyFix HH1
    if ( entropyFix == "HH1" ) 
    {
        // Harten and Hyman 1 entropy fix (HH1)
        delta = 0.1;
        if ( ( lambda1_EF - lambda1_L ) > delta ) delta = lambda1_EF - lambda1_L;
        if ( ( lambda1_R - lambda1_EF ) > delta ) delta = lambda1_R - lambda1_EF;
        lambda1_EF = mag( lambda1_EF );
        if ( lambda1_EF < delta ) lambda1_EF = delta;
        
        delta = 0.1;
        if ( ( lambda2_EF - lambda2_L ) > delta ) delta = lambda2_EF - lambda2_L;
        if ( ( lambda2_R - lambda2_EF ) > delta ) delta = lambda2_R - lambda2_EF;
        lambda2_EF = mag( lambda2_EF );
        if ( lambda2_EF < delta ) lambda2_EF = delta;
        
        delta = 0.1;
        if ( ( lambda3_EF - lambda3_L ) > delta ) delta = lambda3_EF - lambda3_L;
        if ( ( lambda3_R - lambda3_EF ) > delta ) delta = lambda3_R - lambda3_EF;
        lambda3_EF = mag( lambda3_EF );
        if ( lambda3_EF < delta ) lambda3_EF = delta;
    }
    //FI} entropyFix HH1
    //IF{ entropyFix HH2
    else if ( entropyFix == "HH2" )
    {  
        // Harten and Hyman 2 entropy fix (HH2)
        delta = 0.1;
        //delta = 0.1*c_hat;
        if ( ( lambda1_EF - lambda1_L ) > delta ) delta = lambda1_EF - lambda1_L;
        if ( ( lambda1_R - lambda1_EF ) > delta ) delta = lambda1_R - lambda1_EF;
        lambda1_EF = mag( lambda1_EF );
        if ( lambda1_EF < delta ) lambda1_EF = 0.5*( sqr(lambda1_EF)/delta + delta );
        
        delta = 0.1;
        //delta = 0.1*c_hat;
        if ( ( lambda2_EF - lambda2_L ) > delta ) delta = lambda2_EF - lambda2_L;
        if ( ( lambda2_R - lambda2_EF ) > delta ) delta = lambda2_R - lambda2_EF;
        lambda2_EF = mag( lambda2_EF );
        if ( lambda2_EF < delta ) lambda2_EF = 0.5*( sqr(lambda2_EF)/delta + delta );
        
        delta = 0.1;
        //delta = 0.1*c_hat;
        if ( ( lambda3_EF - lambda3_L ) > delta ) delta = lambda3_EF - lambda3_L;
        if ( ( lambda3_R - lambda3_EF ) > delta ) delta = lambda3_R - lambda3_EF;
        lambda3_EF = mag( lambda3_EF );
        if ( lambda3_EF < delta ) lambda3_EF = 0.5*( sqr(lambda3_EF)/delta + delta );
    }   
    //FI} entropyFix HH2
    //IF{ entropyFix HH2b 
    else // if ( entropyFix == "HH2b" ) DEFAULT
    {    
         // Harten and Hyman 2b entropy fix (HH2b), correction on the local Mach number
        delta = 0.1*( u_hat/c_hat + 1 );
        lambda1_EF = mag( lambda1_EF );
        lambda2_EF = mag( lambda2_EF );
        lambda3_EF = mag( lambda3_EF );
        lambda1_EF = ( 0.5 + 0.5*Foam::sign( lambda1_EF - delta ) )*lambda1_EF + 
                     ( 0.5 - 0.5*Foam::sign( lambda1_EF - delta ) )*( sqr( lambda1_EF ) + sqr( delta ) )/( 2*delta + eps );
        lambda2_EF = ( 0.5 + 0.5*Foam::sign( lambda2_EF - delta ) )*lambda2_EF + 
                     ( 0.5 - 0.5*Foam::sign( lambda2_EF - delta ) )*( sqr( lambda2_EF ) + sqr( delta ) )/( 2*delta + eps );
        lambda3_EF = ( 0.5 + 0.5*Foam::sign( lambda3_EF - delta ) )*lambda3_EF + 
                     ( 0.5 - 0.5*Foam::sign( lambda3_EF - delta ) )*( sqr( lambda3_EF ) + sqr( delta ) )/( 2*delta + eps );       
    }
    //FI} entropyFix HH2b 
        
    // Assembly upwind fluxes
    Fuprho      = lambda1_EF*( L11*DU1 + L12*DU2 + L13*DU3 + L14*DU4 + L15*DU5 );
    Fupm.x()    = lambda2_EF*( L21*DU1 + L22*DU2 + L23*DU3 + L24*DU4 + L25*DU5 );
    Fupm.y()    = lambda2_EF*( L31*DU1 + L32*DU2 + L33*DU3 + L34*DU4 + L35*DU5 );
    Fupm.z()    = lambda2_EF*( L41*DU1 + L42*DU2 + L43*DU3 + L44*DU4 + L45*DU5 );
    Fupet_tilde = lambda3_EF*( L51*DU1 + L52*DU2 + L53*DU3 + L54*DU4 + L55*DU5 );
    
    DU1 = Fuprho;
    DU2 = Fupm.x();
    DU3 = Fupm.y();
    DU4 = Fupm.z();
    DU5 = Fupet_tilde;
    
    Fuprho      = -0.5*( R11*DU1 + R12*DU2 + R13*DU3 + R14*DU4 + R15*DU5 );
    Fupm.x()    = -0.5*( R21*DU1 + R22*DU2 + R23*DU3 + R24*DU4 + R25*DU5 );
    Fupm.y()    = -0.5*( R31*DU1 + R32*DU2 + R33*DU3 + R34*DU4 + R35*DU5 );
    Fupm.z()    = -0.5*( R41*DU1 + R42*DU2 + R43*DU3 + R44*DU4 + R45*DU5 );
    Fupet_tilde = -0.5*( R51*DU1 + R52*DU2 + R53*DU3 + R54*DU4 + R55*DU5 );
    
    //==========================================================================    
    // *** ASSEMBLY FLUXES *** 
    //==========================================================================         
    *Frho      = Fsyrho      + Fuprho;
    if ( mag(nz) < 1.0 )
    {
    (*Fm).x()  = Fsym.x()    + Fupm.x()*nx - Fupm.y()*ny/( Foam::sqrt( 1 - sqr(nz) ) ) - Fupm.z()*nx*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).y()  = Fsym.y()    + Fupm.x()*ny + Fupm.y()*nx/( Foam::sqrt( 1 - sqr(nz) ) ) - Fupm.z()*ny*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).z()  = Fsym.z()    + Fupm.x()*nz + Fupm.z()*( Foam::sqrt( 1 - sqr(nz) ) );
    }
    else
    {
    (*Fm).x()  = Fsym.x();
    (*Fm).y()  = Fsym.y();
    (*Fm).z()  = Fsym.z()    + Fupm.x()*nz;
    }
    *Fet_tilde = Fsyet_tilde + Fupet_tilde;
    
    return(0);
    
}    
