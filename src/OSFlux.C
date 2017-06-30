//============================================================================//
//                                                                            //
//                             *** EVAL_FP ***                                //
//                                                                            //
//============================================================================//
int eval_Fp( scalar sign, scalar gamma_, scalar rho, scalar u, scalar v, scalar w, scalar p, scalar *Frho, scalar *Fmx, scalar *Fmy, scalar *Fmz, scalar *Fet_tilde )
{
  
    (*Frho)      = (*Frho)      + sign*( rho*u );
    (*Fmx)       = (*Fmx)       + sign*( rho*u*u + p );
    (*Fmy)       = (*Fmy)       + sign*( rho*u*v );
    (*Fmz)       = (*Fmz)       + sign*( rho*u*w );
    (*Fet_tilde) = (*Fet_tilde) + sign*( u*( p/( gamma_ - 1 ) + 0.5*rho*( sqr(u) + sqr(v) + sqr(w) ) + p ) ); 
    return(0);

}

//============================================================================//
//                                                                            //
//                             *** OS_FLUX ***                                //
//                                                                            //
//============================================================================//
int os_flux( scalar gamma_, scalar rho_L, scalar rho_R, vector mm_L, vector mm_R, scalar et_tilde_L, scalar et_tilde_R, vector n_i, 
              scalar *Frho, vector *Fm, scalar *Fet_tilde )
{  
    
    // Variables definition
    // Reference frame
    scalar nx, ny, nz, norm;
    scalar ux_L, uy_L, uz_L, ux_R, uy_R, uz_R;
    
    // L and R variables
    scalar  u_L,  v_L,  w_L,  u_R,  v_R,  w_R;
    scalar magU_L, magU_R, p_L, p_R, c_L, c_R;
    scalar e_tilde_L, e_tilde_R, s_L, s_R;
    
    // 1 and 4 variables
    scalar  u_1,  v_1,  w_1,  u_4,  v_4,  w_4;
    scalar  p_1,  p_4,  c_1,  c_4,  z_L,  z_R,  rho_1,  rho_4;
    
    // L* and R* variables
    scalar  u_Ls,  v_Ls,  w_Ls,  u_Rs,  v_Rs,  w_Rs;
    scalar  p_Ls,  p_Rs,  c_Ls,  c_Rs,  rho_Ls,  rho_Rs;
    
    // Building Tree
    scalar alpha, beta, gamma, delta, eps, omega;
    
    // Fluxes
    scalar Fosrho, Foset_tilde;
    vector Fosm;   
       
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
    c_L = Foam::sqrt( gamma_*p_L/rho_L );
    c_R = Foam::sqrt( gamma_*p_R/rho_R );
    s_L = p_L/Foam::pow(rho_L, gamma_); 
    s_R = p_R/Foam::pow(rho_R, gamma_); 
    
    // 1 and 4 variables 
    z_L   = 0.5*( gamma_ - 1 )*u_L + c_L;
    z_R   = 0.5*( gamma_ - 1 )*u_R - c_R;
    alpha = Foam::pow( s_R/s_L, 0.5/gamma_ );
    c_1   = ( z_L - z_R )/( 1 + alpha );
    c_4   = alpha*c_1;
    u_1   = 2/( gamma_ - 1 )*( z_L - c_1 );
    u_4   = u_1;
    v_1   = v_L;
    v_4   = v_R;
    w_1   = w_L;
    w_4   = w_R;
    rho_1 = Foam::pow( ( c_1/c_L ), 2.0/( gamma_ - 1 ) )*rho_L;
    rho_4 = Foam::pow( alpha, -2.0 )*rho_1;
    p_1   = rho_1*c_1*c_1/gamma_;
    p_4   = rho_4*c_4*c_4/gamma_;
    
    // L* and R* variables
    c_Ls   = 2/( gamma_ + 1 )*z_L; 
    c_Rs   = -2/( gamma_ + 1 )*z_R; 
    rho_Ls = Foam::pow( ( c_Ls/c_L ), 2.0/( gamma_ - 1 ) )*rho_L;
    rho_Rs = Foam::pow( ( c_Rs/c_R ), 2.0/( gamma_ - 1 ) )*rho_R;
    u_Ls   = c_Ls; 
    u_Rs   = -c_Rs;
    v_Ls   = v_L;
    v_Rs   = v_R;
    w_Ls   = w_L;
    w_Rs   = w_R;
    p_Ls   = rho_Ls*c_Ls*c_Ls/gamma_;
    p_Rs   = rho_Rs*c_Rs*c_Rs/gamma_;  
        
    // Initialization ( necessary since eval_Fp updates and does not overwrite elements )
    Fosrho      = 0.0;
    Fosm.x()    = 0.0;
    Fosm.y()    = 0.0;
    Fosm.z()    = 0.0;
    Foset_tilde = 0.0;
    
    // Build Osher-Solomon fluxes ( see tree of pag 228 of "Mathematical and Computational Methods for Compressible Flows" )   
    alpha = u_L - c_L;
    beta  = u_1 - c_1;
    gamma = u_1;
    delta = u_4 + c_4;
    omega = u_R + c_R; 
    if ( alpha >= 0 )
    {
        if ( beta >= 0 )
        {
            if ( omega >= 0 )
            {
                eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );   
                //Info << "1\n";
            }
            else
            {
                eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                eval_Fp( -1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                //Info << "2\n";
            }
        }
        else
        {   
            if ( gamma >= 0 )
            {
                if ( omega >= 0 )
                {
                    eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( 1.0, gamma_, rho_1, u_1, v_1, w_1, p_1, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( -1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    //Info << "3\n";
                }
                else
                {
                    eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( 1.0, gamma_, rho_1, u_1, v_1, w_1, p_1, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( -1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                    eval_Fp( -1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                    //Info << "4\n";
                }
            }
            else
            {
                if ( delta >= 0 )
                {
                    if ( omega >= 0 )
                    {
                        eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( 1.0, gamma_, rho_4, u_4, v_4, w_4, p_4, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( -1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        //Info << "5\n";
                    }
                    else
                    {
                        eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( 1.0, gamma_, rho_4, u_4, v_4, w_4, p_4, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( -1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        eval_Fp( -1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        //Info << "6\n";
                    }
                }
                else
                {
                    if ( omega >= 0 )
                    {
                        eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( 1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        eval_Fp( -1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        //Info << "7\n";
                    }
                    else
                    {
                        eval_Fp( 1.0, gamma_, rho_L, u_L, v_L, w_L, p_L, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        eval_Fp( -1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        //Info << "8\n";
                    }
                }
            }
        }
    }
    else
    {
        if ( beta >= 0 )
        {
            if ( omega >= 0 )
            {
                eval_Fp( 1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                //Info << "9\n";  
            }
            else
            {
                eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                eval_Fp( 1.0, gamma_, rho_Ls, u_Ls, v_Ls, w_Ls, p_Ls, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                eval_Fp( -1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                //Info << "10\n";
            }
        }
        else
        {   
            if ( gamma >= 0 )
            {
                if ( omega >= 0 )
                {
                    eval_Fp( 1.0, gamma_, rho_1, u_1, v_1, w_1, p_1, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );   
                    //Info << "11\n";
                }
                else
                {
                    eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( 1.0, gamma_, rho_1, u_1, v_1, w_1, p_1, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    eval_Fp( -1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                    //Info << "12\n";
                }
            }
            else
            {
                if ( delta >= 0 )
                {
                    if ( omega >= 0 )
                    {
                        eval_Fp( 1.0, gamma_, rho_4, u_4, v_4, w_4, p_4, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        //Info << "13\n"; 
                    }
                    else
                    {
                        eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( 1.0, gamma_, rho_4, u_4, v_4, w_4, p_4, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        eval_Fp( -1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde ); 
                        //Info << "14\n";
                    }
                }
                else
                {
                    if ( omega >= 0 )
                    {
                        eval_Fp( 1.0, gamma_, rho_Rs, u_Rs, v_Rs, w_Rs, p_Rs, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );  
                        //Info << "15\n";
                    }
                    else
                    {
                        eval_Fp( 1.0, gamma_, rho_R, u_R, v_R, w_R, p_R, &Fosrho, &Fosm.x(), &Fosm.y(), &Fosm.z(), &Foset_tilde );
                        //Info << "16\n";
                    }
                }
            }
        }    
    }
    
    // Back in the Cartesian reference frame 
    *Frho      = Fosrho;
    if ( mag(nz) < 1.0 )
    {
    (*Fm).x()  = Fosm.x()*nx - Fosm.y()*ny/( Foam::sqrt( 1 - sqr(nz) ) ) - Fosm.z()*nx*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).y()  = Fosm.x()*ny + Fosm.y()*nx/( Foam::sqrt( 1 - sqr(nz) ) ) - Fosm.z()*ny*nz/( Foam::sqrt( 1 - sqr(nz) ) );
    (*Fm).z()  = Fosm.x()*nz + Fosm.z()*( Foam::sqrt( 1 - sqr(nz) ) );
    }
    else
    {
    (*Fm).x()  = 0.0;
    (*Fm).y()  = 0.0;
    (*Fm).z()  = Fosm.x()*nz;
    }
    *Fet_tilde = Foset_tilde;
    
    return(0);
    
}    
