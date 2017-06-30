//============================================================================//
//                                                                            //
//                            *** SUTHERLAND ***                              //
//                                                                            //
//============================================================================//
scalar Sutherland( scalar T )
{
    // Variables definition
    scalar mu, mu_0, T_0, S;
    
    // Standard Air
    mu_0 = 1.78e-5; // [ Pa*s ] 
    T_0  = 288.15;  // [ K ] 
    S    = 110.00;  // [ K ]     

    // Power law
    mu = mu_0*( S + T_0 )/( S + T )*Foam::pow( T/T_0, 1.5 );
    
    // Return
    return( mu );
}

//============================================================================//
//                                                                            //
//                           *** VISCOUS_FLUX ***                             //
//                                                                            //
//============================================================================//
int viscous_flux( scalar gamma_, scalar R_, scalar Pr_,                                                       // Constants
                  vector gradUx, vector gradUy, vector gradUz, vector gradT, vector U, scalar T, vector n_ij, // Input variables
                  scalar *Frho_, vector *Fm_, scalar *Fet_tilde_ )                                            // Output variables
{
    // Variables definition
    scalar mu, lambda, Cp_, kappa;
    scalar divU;
    vector StrainRateU;
    vector Gm_;
    scalar Get_tilde_;
    
    // Thermodynamics
    mu     = Sutherland( T );            // Sutherland power law
    lambda = -2.0/3.0*mu;                // Stokes hypothesis
    Cp_    = R_*gamma_/( gamma_ - 1.0 ); // Polytropic Ideal Gas
    kappa  = Cp_*mu/Pr_;                 // Constant Prandtl number flows 
       
    // Strain rate tensor and divergence scalar fields
    StrainRateU.x() = gradUx.x()*n_ij.x() + gradUx.y()*n_ij.y() + gradUx.z()*n_ij.z() + 
                      gradUx.x()*n_ij.x() + gradUy.x()*n_ij.y() + gradUz.x()*n_ij.z();
    StrainRateU.y() = gradUy.x()*n_ij.x() + gradUy.y()*n_ij.y() + gradUy.z()*n_ij.z() + 
                      gradUx.y()*n_ij.x() + gradUy.y()*n_ij.y() + gradUz.y()*n_ij.z();  
    StrainRateU.z() = gradUz.x()*n_ij.x() + gradUz.y()*n_ij.y() + gradUz.z()*n_ij.z() + 
                      gradUx.z()*n_ij.x() + gradUy.z()*n_ij.y() + gradUz.z()*n_ij.z(); 
    divU = gradUx.x() + gradUy.y() + gradUz.z();                  
    
    // Build viscous fluxes *** CHECK FOR ERRORS ***
    Gm_        = mu*StrainRateU + lambda*divU*n_ij;
    Get_tilde_ = ( Gm_ & U ) + kappa*( gradT & n_ij ); 

    // Update
    //(*Frho_)      = (*Frho_)     - 0.0;
    (*Fm_)        = (*Fm_)        - Gm_;
    (*Fet_tilde_) = (*Fet_tilde_) - Get_tilde_;
    
    // Return
    return(0);
}
