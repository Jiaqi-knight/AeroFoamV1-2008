//============================================================================//
//                                                                            //
//                             *** BUILD_FLUXES ***                           //
//                                                                            //
//============================================================================//
int build_fluxes( scalarField *localDt, 
                  word flowType, dimensionedScalar *R, dimensionedScalar *Cv, dimensionedScalar *Pr, 
                  fvMesh *mesh, labelField *id_LL_i, labelField *id_RR_i,
                  volScalarField     *rrho,  volVectorField     *mm,  volScalarField     *eet_tilde, 
                  volScalarField     *p,     volVectorField     *U,   volScalarField     *T, 
                  volVectorField     *gradUx, volVectorField    *gradUy, volVectorField  *gradUz, volVectorField *gradT,
                  surfaceScalarField *FFrho, surfaceVectorField *FFm, surfaceScalarField *FFet_tilde, 
                  word MonoFlux, word HiReFlux, word fluxLim, word entropyFix, scalar extrapolateBC )
{

    // Variables definition
    label  id_L, id_R, id_LL, id_RR, i;                      // Mesh faces, cells connectivity
    vector CG_L, CG_R, CG_LL, CG_RR, Cf, n_i_;               //  
    scalar rho_L_, et_tilde_L_, rho_R_, et_tilde_R_;         // Temporary conservative variables L, R
    vector m_L_, m_R_;                                       //
    scalar rho_LL_, et_tilde_LL_, rho_RR_, et_tilde_RR_;     // Temporary conservative variables LL, RR for high resolution
    vector m_LL_, m_RR_;                                     //
    scalar rho_B_, et_tilde_B_;                              // Extrapolated conservative variables at boundary
    vector m_B_;                                             //
    scalar dx_, slope_rho, slope_m, slope_et_tilde;          // Linear extrapolation of jump variables R, RR at boundary
    scalar phi_rho, phi_m, phi_et_tilde;                     // 
    scalar patchS, faceS;                                    // Boundary sub-patch id ( e.g. for swirl BCs, may be generalized )
    vector patchCG;                                          // 
    scalar Frho_, Fet_tilde_;                                // Temporary fluxes 
    vector Fm_;                                              //
    scalar p_L_, p_R_, T_L_, T_R_, c_R_, Un_R_;              // Temporary primitive variables
    vector U_L_, U_R_;                                       //
    scalar gamma_ = ( 1 + (*R).value()/(*Cv).value() );      // Thermodynamic properties  
    scalar R_ = (*R).value(), Cv_ = (*Cv).value();           //
    vector NablaUx_, NablaUy_, NablaUz_, NablaT_, U_;        // Viscous fluxes 
    scalar T_, Pr_ = (*Pr).value();                          // 
    scalar epsDist, epsGrad, epsPosi;                        // Epsilon
                    
    // Gradient of the primitive variables to compute viscous fluxes
    //IF{ flowType NS RANS
    if ( flowType == "NS" )       
    { 
        (*gradUx) = fvc::grad( (*mm).component(0)/(*rrho) );
        (*gradUy) = fvc::grad( (*mm).component(1)/(*rrho) );    
        (*gradUz) = fvc::grad( (*mm).component(2)/(*rrho) );   
        (*gradT)  = fvc::grad( ( (*eet_tilde) - 0.5*magSqr( (*mm) )/(*rrho) )/( (*rrho)*(*Cv) ) );
    } 
    //FI} flowType NS RANS
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        
    //==========================================================================       
    // A) Loop on internal faces [index (i)]
    //==========================================================================
    forAll( (*mesh).Sf(), i )
    {
        // Input 
        id_L        = (*mesh).faceOwner()[i];
        id_R        = (*mesh).faceNeighbour()[i];   
        rho_L_      = (*rrho)[id_L];
        rho_R_      = (*rrho)[id_R];
        m_L_        = (*mm)[id_L];
        m_R_        = (*mm)[id_R];
        et_tilde_L_ = (*eet_tilde)[id_L];
        et_tilde_R_ = (*eet_tilde)[id_R];
        n_i_        = (*mesh).Sf()[i]/(*mesh).magSf()[i];  
        
        // Find id_RR = (j + 2) and id_LL = (j - 1) ( see Jameson ) 
        id_LL        = (*id_LL_i)[i];
        id_RR        = (*id_RR_i)[i];
        rho_LL_      = (*rrho)[id_LL];
        rho_RR_      = (*rrho)[id_RR];
        m_LL_        = (*mm)[id_LL];
        m_RR_        = (*mm)[id_RR];
        et_tilde_LL_ = (*eet_tilde)[id_LL];
        et_tilde_RR_ = (*eet_tilde)[id_RR];
        
        // Cells CGs
        CG_L  = (*mesh).C()[id_L];
        CG_R  = (*mesh).C()[id_R];
        CG_LL = (*mesh).C()[id_LL];
        CG_RR = (*mesh).C()[id_RR];
        Cf    = (*mesh).Cf()[i];
           
        //---------------------------------------------------------------------- 
        // *** First order fluxes ***
        //----------------------------------------------------------------------
        //IF{ MonotoneFlux Roe
        if ( MonoFlux == "Roe" )
        {
            // Call roe_flux
            roe_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_, entropyFix );
        }
        //FI} MonotoneFlux Roe
        //IF{ MonotoneFlux AUSM
        else if ( MonoFlux == "AUSM" )
        {
            // Call ausm_flux
            ausm_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_ ); 
        }
        //FI} MonotoneFlux AUSM
        //IF{ MonotoneFlux CUSP
        else if ( MonoFlux == "CUSP" )
        {
            // Call cusp_flux
            cusp_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_ ); 
        }
        //FI} MonotoneFlux CUSP
        //IF{ MonotoneFlux HLL
        else if ( MonoFlux == "HLL" )
        {
            // Call hll_flux
            hll_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_, entropyFix ); 
        } 
        //FI} MonotoneFlux HLL
        //IF{ MonotoneFlux HLLC
        else if ( MonoFlux == "HLLC" )
        {
            // Call hllc_flux
            hllc_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_, entropyFix); 
        }
        //FI} MonotoneFlux HLLC
        //IF{ MonotoneFlux OS
        else if ( MonoFlux == "OS" )
        {
            // Call os_flux
            os_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_ );
        }
        //FI} MonotoneFlux OS
        
        //----------------------------------------------------------------------
        // *** High resolution fluxes ***
        //---------------------------------------------------------------------- 
        //IF{ HighResolutionFlux LW
        if ( HiReFlux == "LW" )
        {
            // Call lwhr_flux
            lwhr_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, 
                               rho_LL_, rho_RR_, m_LL_, m_RR_, et_tilde_LL_, et_tilde_RR_,
                               n_i_, 0.5*( (*localDt)[id_L] + (*localDt)[id_R] ), mag( CG_R - CG_L ),
                               &Frho_, &Fm_, &Fet_tilde_, fluxLim ); 
        }
        //FI} HighResolutionFlux LW
        //IF{ HighResolutionFlux JST
        else if ( HiReFlux == "JST" )
        {
            // Call jsthr_flux
            jsthr_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, 
                                rho_LL_, rho_RR_, m_LL_, m_RR_, et_tilde_LL_, et_tilde_RR_,
                                n_i_, &Frho_, &Fm_, &Fet_tilde_ );         
        }
        //FI} HighResolutionFlux JST
                
        //----------------------------------------------------------------------
        // *** Viscous fluxes for NS ***
        //----------------------------------------------------------------------
        //IF{ flowType NS RANS
        if ( flowType == "NS" )
        {
            NablaUx_ = 0.5*( (*gradUx)[id_L] + (*gradUx)[id_R] );
            NablaUy_ = 0.5*( (*gradUy)[id_L] + (*gradUy)[id_R] );
            NablaUz_ = 0.5*( (*gradUz)[id_L] + (*gradUz)[id_R] );
            NablaT_  = 0.5*( (*gradT)[id_L]  + (*gradT)[id_R]  );
            U_L_     = m_L_/rho_L_;
            U_R_     = m_R_/rho_R_;
            U_       = 0.5*( U_L_ + U_R_ );
            T_L_     = ( et_tilde_L_ - 0.5*magSqr(m_L_)/rho_L_ )/( rho_L_*Cv_ );
            T_R_     = ( et_tilde_R_ - 0.5*magSqr(m_R_)/rho_R_ )/( rho_R_*Cv_ );
            T_       = 0.5*( T_L_ + T_R_ );
            viscous_flux( gamma_, R_, Pr_, NablaUx_, NablaUy_, NablaUz_, NablaT_, U_, T_, 
                                  n_i_, &Frho_, &Fm_, &Fet_tilde_ );   
        }
        //FI} flowType NS RANS
        
        // Output            
        (*FFrho)[i]      = Frho_;
        (*FFm)[i]        = Fm_;
        (*FFet_tilde)[i] = Fet_tilde_;

    } 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    
    //==========================================================================
    // B) Loop on boundary faces and BCs [local index(ii), global index(i)]
    //==========================================================================
    forAll( (*mesh).boundaryMesh(), iPatch )
    {
    
        // Patch BCtype
        word BCtype = (*mesh).boundaryMesh().physicalTypes()[iPatch];
        
        // For swirl/other boundary conditions compute the patch CG 
        // (weighted with the face areas to increase accuracy)
        patchS = 0.0;
        patchCG = vector(0.0, 0.0, 0.0);
        if ( BCtype == "swirl" )  
        {
            forAll( (*mesh).boundaryMesh()[iPatch].faceAreas(), ii )
            {
                i = (*mesh).boundaryMesh()[iPatch].start()+ii;
                faceS   = mag((*mesh).boundaryMesh()[iPatch].faceAreas()[ii]);
                patchS  = patchS  + faceS;
                patchCG = patchCG + faceS*(*mesh).Cf()[i];
            }
            patchCG = patchCG/patchS;
        }
    
        //======================================================================
        // Loop on patch BCtype faces 
        //======================================================================
        forAll( (*mesh).boundaryMesh()[iPatch].faceAreas(), ii )
        {
        
            // *** BCs ( Boundary type must not be empty ) ***
            if ( (*FFrho).boundaryField()[iPatch].size() > 0 )
            {
        
                // Input
                id_L        = (*mesh).boundaryMesh()[iPatch].faceCells()[ii];    
                id_R        = id_L;          
                i           = (*mesh).boundaryMesh()[iPatch].start()+ii;
                rho_L_      = (*rrho)[id_L];  
                m_L_        = (*mm)[id_L];
                et_tilde_L_ = (*eet_tilde)[id_L];
                U_L_        = m_L_/rho_L_;
                T_L_        = ( et_tilde_L_ - 0.5*rho_L_*magSqr(U_L_) )/( rho_L_*Cv_ );
                p_L_        = rho_L_*R_*T_L_;
                n_i_        = (*mesh).boundaryMesh()[iPatch].faceAreas()[ii];
                n_i_        = n_i_/mag(n_i_);
                
                // Find id_LL = (j - 1) ( see Jameson ) 
                id_LL        = (*id_LL_i)[i];
                id_RR        = id_R;
                rho_LL_      = (*rrho)[id_LL];
                m_LL_        = (*mm)[id_LL];
                et_tilde_LL_ = (*eet_tilde)[id_LL];
                
                // Cells CGs
                CG_L  = (*mesh).C()[id_L];
                CG_R  = 2*(*mesh).Cf()[i] - (*mesh).C()[id_L];
                CG_LL = (*mesh).C()[id_LL];
                CG_RR = 2*(*mesh).Cf()[i] - (*mesh).C()[id_LL];
                Cf    = (*mesh).Cf()[i];
                          
                //--------------------------------------------------------------              
                // *** Ghost cells ***
                //--------------------------------------------------------------
                // Default: > ExtrapolatedOutlet ( p, T, U zeroGradient )
                //          > Constant or linear extrapolation
                // REMARKS: - compatibility with High Resolution procedures?
                //          - distance should be projected along boundary face normal. 
                //            This may lead to division by zero and a check must be
                //            added to increase code robustness ( linear extrapolation  
                //            is not accurate when dx_ << 1 )
                //**************************************************************
                epsDist = 0.5;
                epsGrad = 0.25;
                epsPosi = 0.05;  
                //**************************************************************
                //dx_ = mag( ( CG_L - CG_LL ) & n_i_ )/mag( ( Cf   - CG_L ) & n_i_ );  
                dx_ = mag( CG_L - CG_LL )/mag( Cf   - CG_L );               
                if ( ( dx_ < epsDist ) || ( extrapolateBC == 0.0 )  )
                {
                    // Constant extrapolation ( *** TODO: Check how accuracy is affected *** )
                    rho_B_      = rho_L_;
                    m_B_        = m_L_;
                    et_tilde_B_ = et_tilde_L_;
                }
                else
                {
                    // Linear limited extrapolation ( *** TODO: Check how accuracy is affected *** ) 
                    slope_rho      = mag( rho_L_      - rho_LL_      )/dx_/( mag( rho_LL_      ) + 1e-12 );
                    slope_m        = mag( m_L_        - m_LL_        )/dx_/( mag( m_LL_        ) + 1e-12 );
                    slope_et_tilde = mag( et_tilde_L_ - et_tilde_LL_ )/dx_/( mag( et_tilde_LL_ ) + 1e-12 );   
                    phi_rho      = 1.0;
                    phi_m        = 1.0;
                    phi_et_tilde = 1.0; 
                    if ( slope_rho      > epsGrad ) phi_rho      = epsGrad/slope_rho;   
                    if ( slope_m        > epsGrad ) phi_m        = epsGrad/slope_m;  
                    if ( slope_et_tilde > epsGrad ) phi_et_tilde = epsGrad/slope_et_tilde; 
                    rho_B_      = rho_L_      +      phi_rho*( rho_L_      - rho_LL_      )/dx_;     
                    m_B_        = m_L_        +        phi_m*( m_L_        - m_LL_        )/dx_; 
                    et_tilde_B_ = et_tilde_L_ + phi_et_tilde*( et_tilde_L_ - et_tilde_LL_ )/dx_; 
                }

                // Check positivity ( add mag(*) to increase code robustness )
                if ( ( rho_B_ <= epsPosi*rho_L_ ) || ( et_tilde_B_ <= epsPosi*et_tilde_L_ ) ) 
                {
                   rho_B_      = rho_L_;
                   m_B_        = m_L_;
                   et_tilde_B_ = et_tilde_L_;
                }
                
                // Default initialization: constant extrapolation
                rho_R_      = rho_L_;
                m_R_        = m_L_;
                et_tilde_R_ = et_tilde_L_;
                rho_RR_      = rho_R_;
                m_RR_        = m_R_;
                et_tilde_RR_ = et_tilde_R_; 
                                        
                //--------------------------------------------------------------                     
                // *** BCs list ***
                //--------------------------------------------------------------
                // 1) SupersonicInlet ( p, T, U fixedValue )
                if ( BCtype == "supersonicInlet" )  
                {
                    // 1st ghost cell
                    p_R_        = (*p).boundaryField()[iPatch][ii];
                    T_R_        = (*T).boundaryField()[iPatch][ii];
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ ); 
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;                       
                }
                // 2) ExtrapolatedOutlet ( p, T, U zeroGradient )
                else if ( BCtype == "extrapolatedOutlet" )  
                {
                    // 1st ghost cell
                    rho_R_      = rho_L_;
                    m_R_        = m_L_;    
                    et_tilde_R_ = et_tilde_L_;
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;    
                    et_tilde_RR_ = et_tilde_R_;   
                }
                // 3) Slip, symmetryPlane, inviscidWall ( p, T, v, w zeroGradient, u fixedValue = 0 )
                else if ( ( BCtype == "slip" ) || ( BCtype == "symmetryPlane" ) || ( BCtype == "inviscidWall" ) )  
                {
                    // 1st ghost cell
                    rho_R_      = rho_L_;
                    m_R_        = m_L_ - 2*( m_L_ & n_i_ )*n_i_;    
                    et_tilde_R_ = et_tilde_L_ - 0.5*magSqr( m_L_ )/rho_L_ + 0.5*magSqr( m_R_ )/rho_R_;
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;
                }
                // 4) Transpiration ( p, T, v, w zeroGradient, u fixedValue = Vbn )
                // Vbn is expressed in the face local frame of reference that is U.x() = u, U.y() = v, U.z() = w 
                else if ( BCtype == "transpiration" )  
                {    
                    // 1st ghost cell 
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    rho_R_      = rho_L_;
                    m_R_        = m_L_ - 2*( m_L_ & n_i_ )*n_i_ + 2*rho_R_*U_R_.x()*n_i_;   
                    et_tilde_R_ = et_tilde_L_ - 0.5*magSqr( m_L_ )/rho_L_ + 0.5*magSqr( m_R_ )/rho_R_;
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;  
                }  
                // 5) ViscousAdiabaticWall ( p, T zeroGradient, U fixedValue = 0 )
                else if ( BCtype == "viscousAdiabaticWall" )  
                {
                    // 1st ghost cell
                    rho_R_      = rho_L_;
                    m_R_        = -m_L_;    
                    et_tilde_R_ = et_tilde_L_ - 0.5*magSqr( m_L_ )/rho_L_ + 0.5*magSqr( m_R_ )/rho_R_;
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;
                }  
                // 6) ViscousIsothermalWall ( p zeroGradient, T, U fixedValue )
                else if ( BCtype == "viscousIsothermalWall" )  
                {
                    // 1st ghost cell
                    T_R_        = (*T).boundaryField()[iPatch][ii];
                    rho_R_      = rho_L_;
                    m_R_        = -m_L_;     
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*magSqr( m_R_ )/rho_R_;
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;
                }  
                // 0) Riemann boundary conditions (SI + I + EO + O)
                else if ( BCtype == "Riemann" )
                {
                    // Default: Extrapolation
                    // 1st ghost cell
                    rho_R_      = rho_L_;
                    m_R_        = m_L_;    
                    et_tilde_R_ = et_tilde_L_;
                    U_R_        = m_R_/rho_R_;
                    T_R_        = ( et_tilde_R_ - 0.5*magSqr(m_R_)/rho_R_ )/( rho_R_*Cv_ );
                    p_R_        = rho_R_*R_*T_R_;
                    c_R_        = Foam::sqrt( gamma_*R_*T_R_ );
                    Un_R_       = U_R_ & n_i_;
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;  
                    
                    // *** Inlet ***
                    if ( Un_R_ < 0.0 )
                    {
                        // Supersonic Inlet -> p, T, u fixedValue   
                        if ( mag( Un_R_ ) >= c_R_ )  
                        {
                            // 1st ghost cell
                            p_R_        = (*p).boundaryField()[iPatch][ii];
                            T_R_        = (*T).boundaryField()[iPatch][ii];
                            U_R_        = (*U).boundaryField()[iPatch][ii];
                            rho_R_      = p_R_/( R_*T_R_ );
                            m_R_        = rho_R_*U_R_;                        
                            et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );
                            // 2nd ghost cell
                            rho_RR_      = rho_R_;
                            m_RR_        = m_R_;
                            et_tilde_RR_ = et_tilde_R_;    
                        }
                        // Subsonic Inlet -> T, u fixedValue, p extrapolated
                        else
                        {
                            // 1st ghost cell
                            T_R_        = (*T).boundaryField()[iPatch][ii];
                            U_R_        = (*U).boundaryField()[iPatch][ii];
                            rho_R_      = p_R_/( R_*T_R_ );
                            m_R_        = rho_R_*U_R_;                        
                            et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );
                            // 2nd ghost cell
                            rho_RR_      = rho_R_;
                            m_RR_        = m_R_;
                            et_tilde_RR_ = et_tilde_R_;   
                        }
                    }
                    // *** Outlet ***
                    else
                    {
                        // Supersonic Outlet -> p, T, u extrapolated
                        if ( mag( Un_R_ ) >= c_R_ )  
                        {
                            /*
                            // 1st ghost cell
                            rho_R_      = rho_L_;
                            m_R_        = m_L_;    
                            et_tilde_R_ = et_tilde_L_;
                            // 2nd ghost cell
                            rho_RR_      = rho_R_;
                            m_RR_        = m_R_;
                            et_tilde_RR_ = et_tilde_R_; 
                            */
                        }
                        // Subsonic Outlet -> p fixedValue, T, u extrapolated
                        else
                        {
                            // 1st ghost cell
                            p_R_        = (*p).boundaryField()[iPatch][ii];
                            rho_R_      = p_R_/( R_*T_R_ );
                            m_R_        = rho_R_*U_R_;                        
                            et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );
                            // 2nd ghost cell
                            rho_RR_      = rho_R_;
                            m_RR_        = m_R_;
                            et_tilde_RR_ = et_tilde_R_; 
                        }
                    }                      
                }                
                // 7) ExtrapolatedP, Inlet ( p zeroGradient, T, u fixedValue ) *** CONSTANT EXTRAPOLATION ***
                else if ( ( BCtype == "extrapolatedP" ) || ( BCtype == "inlet" ) )  
                {    
                    // 1st ghost cell
                    p_R_        = p_L_;
                    T_R_        = (*T).boundaryField()[iPatch][ii];
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );   
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;  
                }
                // 8) ExtrapolatedT ( T zeroGradient, p, U fixedValue ) *** CONSTANT EXTRAPOLATION ***
                else if ( BCtype == "extrapolatedT" )  
                {    
                    // 1st ghost cell    
                    p_R_        = (*p).boundaryField()[iPatch][ii];
                    T_R_        = T_L_;
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );    
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;  
                } 
                // 9) ExtrapolatedU ( U zeroGradient, p, T fixedValue ) *** CONSTANT EXTRAPOLATION ***
                else if ( BCtype == "extrapolatedT" )  
                {   
                    // 1st ghost cell     
                    p_R_        = (*p).boundaryField()[iPatch][ii];
                    T_R_        = (*T).boundaryField()[iPatch][ii];
                    U_R_        = U_L_;
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );   
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;   
                } 
                // 10) ExtrapolatedPT, AdiabaticWall ( p, T zeroGradient, U fixedValue ) *** CONSTANT EXTRAPOLATION ***
                else if ( ( BCtype == "extrapolatedPT" ) || ( BCtype == "adiabaticWall" ) )  
                {   
                    // 1st ghost cell     
                    p_R_        = p_L_;
                    T_R_        = T_L_;
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );   
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;   
                }   
                // 11) ExtrapolatedPU ( p, U zeroGradient, T fixedValue ) *** CONSTANT EXTRAPOLATION ***
                else if ( ( BCtype == "extrapolatedPU" ) )  
                {   
                    // 1st ghost cell 
                    p_R_        = p_L_;
                    T_R_        = (*T).boundaryField()[iPatch][ii];
                    U_R_        = U_L_;
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;      
                }  
                // 11) ExtrapolatedTU ( T, U zeroGradient, p fixedValue ) *** CONSTANT EXTRAPOLATION ***
                else if ( ( BCtype == "extrapolatedTU" ) )  
                {    
                    // 1st ghost cell
                    p_R_        = (*p).boundaryField()[iPatch][ii];
                    T_R_        = T_L_;
                    U_R_        = U_L_;
                    rho_R_      = p_R_/( R_*T_R_ );
                    m_R_        = rho_R_*U_R_;                        
                    et_tilde_R_ = rho_R_*Cv_*T_R_ + 0.5*rho_R_*magSqr( U_R_ );
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;      
                }    
                // 12) swirl ( p, T, U zeroGradient and U = U + Vt = U + W x (P_i - CG) ) *** CONSTANT EXTRAPOLATION ***
                else if ( ( BCtype == "swirl" ) )  
                {    
                    // 1st ghost cell
                    rho_R_      = rho_L_;
                    m_R_        = m_L_;    
                    et_tilde_R_ = et_tilde_L_;
                    // Rotation vector omega is temporarily stored in U_R_ vector
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    m_R_        = m_R_ + rho_R_*( U_R_ ^ ( (*mesh).Cf()[i] - patchCG ) );
                    et_tilde_R_ = et_tilde_L_ - 0.5*magSqr( m_L_ )/rho_L_ + 0.5*magSqr( m_R_ )/rho_R_;  
                    // 2nd ghost cell
                    rho_RR_      = rho_R_;
                    m_RR_        = m_R_;
                    et_tilde_RR_ = et_tilde_R_;     
                }    
                // *** TODO: Add a swirl BC where p, U, T and omega are prescribed ***
                 
                //--------------------------------------------------------------
                // *** First order fluxes *** 
                //--------------------------------------------------------------
                //IF{ MonotoneFlux Roe                                                              
                if ( MonoFlux == "Roe" )
                {
                    // Call roe_flux
                    roe_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_, entropyFix );
                }
                //FI} MonotoneFlux Roe    
                //IF{ MonotoneFlux AUSM
                else if ( MonoFlux == "AUSM" )
                {
                    // Call ausm_flux
                    ausm_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_ );  
                }
                //FI} MonotoneFlux AUSM
                //IF{ MonotoneFlux CUSP
                else if ( MonoFlux == "CUSP" )
                {
                    // Call cusp_flux
                    cusp_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_ ); 
                }
                //FI} MonotoneFlux CUSP
                //IF{ MonotoneFlux HLL
                else if ( MonoFlux == "HLL" )
                {
                    // Call hll_flux
                    hll_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_, entropyFix ); 
                }
                //FI} MonotoneFlux HLL
                //IF{ MonotoneFlux HLLC
                else if ( MonoFlux == "HLLC" )
                {
                    // Call hllc_flux
                    hllc_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_, entropyFix ); 
                }
                //FI} MonotoneFlux HLLC
                //IF{ MonotoneFlux OS 
                else if ( MonoFlux == "OS" )
                {
                    // Call os_flux
                    os_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, n_i_, &Frho_, &Fm_, &Fet_tilde_ );
                }
                //FI} MonotoneFlux OS 
                
                //--------------------------------------------------------------
                // *** High resolution fluxes ***
                //--------------------------------------------------------------
                //IF{ HighResolutionFlux LW
                if ( HiReFlux == "LW" )
                {
                    // Call lwhr_flux
                    lwhr_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, 
                                       rho_LL_, rho_RR_, m_LL_, m_RR_, et_tilde_LL_, et_tilde_RR_,
                                       n_i_, (*localDt)[id_L], mag( CG_R - CG_L ),
                                       &Frho_, &Fm_, &Fet_tilde_, fluxLim ); 
                }
                //FI} HighResolutionFlux LW
                //IF{ HighResolutionFlux JST
                else if ( HiReFlux == "JST" )
                {
                    // Call jsthr_flux
                    jsthr_flux( gamma_, rho_L_, rho_R_, m_L_, m_R_, et_tilde_L_, et_tilde_R_, 
                                        rho_LL_, rho_RR_, m_LL_, m_RR_, et_tilde_LL_, et_tilde_RR_,
                                        n_i_, &Frho_, &Fm_, &Fet_tilde_ );         
                }
                //FI} HighResolutionFlux JST
                
                //--------------------------------------------------------------
                // *** Re-enforce BCs for 3, 5, 6 and 4 ***
                //--------------------------------------------------------------
                if ( ( BCtype == "slip" ) || ( BCtype == "symmetryPlane" ) || ( BCtype == "inviscidWall" ) || // 3)  
                     ( BCtype == "viscousAdiabaticWall" ) ||                                                  // 5)
                     ( BCtype == "viscousIsothermalWall" ) )                                                  // 6)
                {    
                    // Only P is needed
                    rho_R_      = rho_B_;
                    m_R_        = m_B_;    
                    et_tilde_R_ = et_tilde_B_;
                    T_R_        = ( et_tilde_R_ - 0.5*magSqr(m_R_)/rho_R_ )/( rho_R_*Cv_ );
                    p_R_        = rho_R_*R_*T_R_;
                    Frho_       = 0.0;  
                    Fm_         = p_R_*n_i_;                                                        
                    Fet_tilde_  = 0.0;
                }
                else if ( BCtype == "transpiration" )
                {
                    // rho, P, et_tilde are needed
                    // REMARK: Transpiration velocity Un_R_ is reinitialized to increase robustness
                    rho_R_      = rho_B_;
                    m_R_        = m_B_;    
                    et_tilde_R_ = et_tilde_B_;
                    U_R_        = (*U).boundaryField()[iPatch][ii];
                    Un_R_       = U_R_.x();
                    T_R_        = ( et_tilde_R_ - 0.5*magSqr(m_R_)/rho_R_ )/( rho_R_*Cv_ );
                    p_R_        = rho_R_*R_*T_R_;
                    Frho_       = rho_R_*Un_R_;  
                    Fm_         = Un_R_*m_R_ + p_R_*n_i_;                                                        
                    Fet_tilde_  = Un_R_*( et_tilde_R_ + p_R_ );
                }
                
                //--------------------------------------------------------------
                // *** Viscous fluxes for NS ***
                //--------------------------------------------------------------
                //IF{ flowType NS RANS
                // WARNING: CHECK CONSISTENCY WITH BCs REINFORCEMENT ABOVE !!!
                 if ( flowType == "NS" )
                {
                    NablaUx_ = (*gradUx).boundaryField()[iPatch][ii];
                    NablaUy_ = (*gradUy).boundaryField()[iPatch][ii];
                    NablaUz_ = (*gradUz).boundaryField()[iPatch][ii];
                    NablaT_  = (*gradT).boundaryField()[iPatch][ii];
                    U_L_     = m_L_/rho_L_;
                    U_R_     = m_R_/rho_R_;
                    U_       = 0.5*( U_L_ + U_R_ );
                    T_L_     = ( et_tilde_L_ - 0.5*magSqr(m_L_)/rho_L_ )/( rho_L_*Cv_ );
                    T_R_     = ( et_tilde_R_ - 0.5*magSqr(m_R_)/rho_R_ )/( rho_R_*Cv_ );
                    T_       = 0.5*( T_L_ + T_R_ );
                    viscous_flux( gamma_, R_, Pr_, NablaUx_, NablaUy_, NablaUz_, NablaT_, U_, T_, 
                                          n_i_, &Frho_, &Fm_, &Fet_tilde_ );  
                }
                //FI} flowType NS RANS
                           
                // Output            
                (*FFrho).boundaryFieldRef()[iPatch][ii]      = Frho_;
                (*FFm).boundaryFieldRef()[iPatch][ii]        = Fm_;
                (*FFet_tilde).boundaryFieldRef()[iPatch][ii] = Fet_tilde_;
            
            }
            
        }
    
    }

    return(0);

}
