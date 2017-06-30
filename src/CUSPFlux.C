//============================================================================//
//                                                                            //
//                             *** CUSP_FLUX ***                              //
//                                                                            //
//============================================================================//
int cusp_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, vector n_i, 
               scalar *Frho, vector *Fm, scalar *Fet_tilde )
{  

    // Variables definition
    scalar eps;
    
    // Reference frame
    scalar nx, ny, nz, norm;
    scalar ux_L, uy_L, uz_L, ux_R, uy_R, uz_R;
    
    // L and R variables
    scalar  u_L,  v_L,  w_L,  u_R,  v_R,  w_R;
    scalar magU_L, magU_R, p_L, p_R, c_L, c_R;
    scalar e_tilde_L, e_tilde_R, ht_L, ht_R;
    
    // 1/2 variables
    scalar u12, c12, M12;
    scalar F1, F2;
    scalar aa0, aa2, aa4;
    
    // aa0 manopolina per la dissipazione
    aa0 = 0.25;
    
    // Fluxes
    scalar Fsyrho, Fsyet_tilde, Fcusprho, Fcuspet_tilde;
    vector Fsym, Fcuspm;   
    
    // Initialization
    Fsyrho      = 0.0;
    Fsym.x()    = 0.0;
    Fsym.y()    = 0.0;
    Fsym.z()    = 0.0;
    Fsyet_tilde = 0.0;
    Fcusprho      = 0.0;
    Fcuspm.x()    = 0.0;
    Fcuspm.y()    = 0.0;
    Fcuspm.z()    = 0.0;
    Fcuspet_tilde = 0.0;
       
    // Change frame of reference 
    eps = 1e-10;
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
    c_L = Foam::sqrt( gamma_*p_L/rho_L );
    c_R = Foam::sqrt( gamma_*p_R/rho_R );
    ht_L = ( p_L + et_tilde_L )/rho_L;
    ht_R = ( p_R + et_tilde_R )/rho_R;
          
    // Assembly symmetric fluxes
    Fsyrho      = 0.5*( rho_L*u_L + rho_R*u_R );
    Fsym.x()    = 0.5*( rho_L*u_L*( ux_L ) + p_L*nx + rho_R*u_R*( ux_R ) + p_R*nx );
    Fsym.y()    = 0.5*( rho_L*u_L*( uy_L ) + p_L*ny + rho_R*u_R*( uy_R ) + p_R*ny );
    Fsym.z()    = 0.5*( rho_L*u_L*( uz_L ) + p_L*nz + rho_R*u_R*( uz_R ) + p_R*nz );
    Fsyet_tilde = 0.5*( u_L*( et_tilde_L + p_L ) + u_R*( et_tilde_R + p_R ) );
    
    // Build dissipative contribution
    u12 = 0.5*( u_L + u_R );
    c12 = 0.5*( c_L + c_R );
    M12 = u12/c12;
    if ( mag(M12) <= 1 )
    {
       aa2 = 1.5 - 2*aa0;
       aa4 = aa0 - 0.5;
       F1 = aa0 + aa2*Foam::pow(M12, 2) + aa4*Foam::pow(M12, 4);
       F2 = 0.5*M12*( 3 - Foam::sqr(M12) );
    }
    else 
    {
       F1 = mag(M12);
       F2 = Foam::sign(M12);
    }
    Fcusprho      = -( F1*c12*( rho_R      - rho_L )                        );
    Fcuspm.x()    = -( F1*c12*( rho_R*u_R  - rho_L*u_L ) + F2*( p_R - p_L ) );
    Fcuspm.y()    = -( F1*c12*( rho_R*v_R  - rho_L*v_L )                    );
    Fcuspm.z()    = -( F1*c12*( rho_R*w_R  - rho_L*w_L )                    );
    Fcuspet_tilde = -( F1*c12*( et_tilde_R - et_tilde_L )                   );
    
     // Back in the Cartesian reference frame 
    *Frho      = Fsyrho      + Fcusprho;
    if ( mag(nz) < 1.0 )
    {
    (*Fm).x()  = Fsym.x()    + Fcuspm.x()*nx - Fcuspm.y()*ny/( Foam::sqrt( 1 - sqr(nz) ) ) - Fcuspm.z()*nx*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).y()  = Fsym.y()    + Fcuspm.x()*ny + Fcuspm.y()*nx/( Foam::sqrt( 1 - sqr(nz) ) ) - Fcuspm.z()*ny*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).z()  = Fsym.z()    + Fcuspm.x()*nz + Fcuspm.z()*( Foam::sqrt( 1 - sqr(nz) ) );
    }
    else
    {
    (*Fm).x()  = Fsym.x();
    (*Fm).y()  = Fsym.y();
    (*Fm).z()  = Fsym.z()    + Fcuspm.x()*nz;
    }
    *Fet_tilde = Fsyet_tilde + Fcuspet_tilde;
          
    return(0);

}
