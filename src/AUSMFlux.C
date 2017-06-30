//============================================================================//
//                                                                            //
//                             *** AUSM_FLUX ***                              //
//                                                                            //
//============================================================================//
int ausm_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, vector n_i, 
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
    scalar e_tilde_L, e_tilde_R, ht_L, ht_R, M_L, M_R;
    
    // +(p) and -(m) variables
    scalar M_Lp, M_Rm, M12;
    scalar p_Lp, p_Rm;
      
    // Fluxes
    scalar Fausmrho, Fausmet_tilde;
    vector Fausmm;   
    
    // Initialization
    Fausmrho      = 0.0; 
    Fausmet_tilde = 0.0;
    Fausmm.x()    = 0.0;
    Fausmm.y()    = 0.0;
    Fausmm.z()    = 0.0;  
     
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
    M_L = u_L/c_L;
    M_R = u_R/c_R;  
    
    if ( ( mag(M_L) < 1 ) || ( mag(M_R) < 1 ) )
    {    
    
       // Build flux
       M_Lp =  0.25*Foam::sqr( M_L + 1 );
       M_Rm = -0.25*Foam::sqr( M_R - 1 );
       M12  = M_Lp + M_Rm;
       // (3)
       //p_Lp = 0.25*p_L*Foam::sqr( M_L + 1 )*( 2 - M_L );
       //p_Rm = 0.25*p_R*Foam::sqr( M_R - 1 )*( 2 + M_R );
       // (1)
       p_Lp = 0.5*p_L*( 1 + M_L );
       p_Rm = 0.5*p_R*( 1 - M_R );
       
       Fausmrho      = 0.5*M12*( rho_R*c_R      + rho_L*c_L     )  - 0.5*mag(M12)*( rho_R*c_R      - rho_L*c_L     );
       Fausmm.x()    = 0.5*M12*( rho_R*c_R*u_R  + rho_L*c_L*u_L )  - 0.5*mag(M12)*( rho_R*c_R*u_R  - rho_L*c_L*u_L ) + ( p_Lp + p_Rm );
       Fausmm.y()    = 0.5*M12*( rho_R*c_R*v_R  + rho_L*c_L*v_L )  - 0.5*mag(M12)*( rho_R*c_R*v_R  - rho_L*c_L*v_L );
       Fausmm.z()    = 0.5*M12*( rho_R*c_R*w_R  + rho_L*c_L*w_L )  - 0.5*mag(M12)*( rho_R*c_R*w_R  - rho_L*c_L*w_L );
       Fausmet_tilde = 0.5*M12*( rho_R*c_R*ht_R + rho_L*c_L*ht_L ) - 0.5*mag(M12)*( rho_R*c_R*ht_R - rho_L*c_L*ht_L );
       
    }
    else
    {
       
       if ( M_L >= 1 )
       {
           Fausmrho      = rho_L*u_L;
           Fausmm.x()    = rho_L*u_L*u_L + p_L;
           Fausmm.y()    = rho_L*u_L*v_L;
           Fausmm.z()    = rho_L*u_L*w_L;
           Fausmet_tilde = rho_L*u_L*ht_L;
       }
       else if ( M_R <= -1 )
       {
           Fausmrho      = rho_R*u_R;
           Fausmm.x()    = rho_R*u_R*u_R + p_R;
           Fausmm.y()    = rho_R*u_R*v_R;
           Fausmm.z()    = rho_R*u_R*w_R;
           Fausmet_tilde = rho_R*u_R*ht_R;
       }
       
    }   
    // Back in the Cartesian reference frame 
    *Frho      = Fausmrho;
    if ( mag(nz) < 1.0 )
    {
    (*Fm).x()  = Fausmm.x()*nx - Fausmm.y()*ny/( Foam::sqrt( 1 - sqr(nz) ) ) - Fausmm.z()*nx*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).y()  = Fausmm.x()*ny + Fausmm.y()*nx/( Foam::sqrt( 1 - sqr(nz) ) ) - Fausmm.z()*ny*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).z()  = Fausmm.x()*nz + Fausmm.z()*( Foam::sqrt( 1 - sqr(nz) ) );
    }
    else
    {
    (*Fm).x()  = 0.0;
    (*Fm).y()  = 0.0;
    (*Fm).z()  = Fausmm.x()*nz;
    }
    *Fet_tilde = Fausmet_tilde;
     
    return(0);

}
