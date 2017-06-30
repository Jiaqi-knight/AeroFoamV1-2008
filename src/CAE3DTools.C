//============================================================================//
//                                                                            //
//                        *** READSTRUCTURALMESH ***                          //
//                                                                            //
//============================================================================//
int CAE3DreadStructralMesh( vectorField *xx_s_v, labelField *id_A_s, labelField *id_B_s, labelField *id_C_s )
{
    // Variables definition
    label k, N;
    label ia, ib, ic;
    scalar x, y, z;
    FILE *f1;
    
    // TODO: *** Read structural mesh data from .mail file of Code_Aster *** [???]
    
    // Read StructuralModel.vertices and initialize memory
    f1 = fopen("./Data/StructuralModel.vertices", "r");
    if ( f1 == NULL ) Info << "ERROR: File StructuralModel.vertices not found!" << nl;
    fscanf(f1, "%i", &N);
    (*xx_s_v) = vectorField(N);
    for( k = 0; k < N; k++ )
    {
        fscanf(f1, "%lf %lf %lf\n", &x, &y, &z );
        (*xx_s_v)[k].x() = x;
        (*xx_s_v)[k].y() = y;
        (*xx_s_v)[k].z() = z;
    } 
    fclose(f1);
    
    // Read StructuralModel.elements and initialize memory
    f1 = fopen("./Data/StructuralModel.elements", "r");
    if ( f1 == NULL ) Info << "ERROR: File StructuralModel.elements not found!" << nl;
    fscanf(f1, "%i", &N);
    (*id_A_s) = labelField(N);
    (*id_B_s) = labelField(N);
    (*id_C_s) = labelField(N);
    for( k = 0; k < N; k++ )
    {
        fscanf(f1, "%i %i %i\n", &ia, &ib, &ic );
        (*id_A_s)[k] = ia;
        (*id_B_s)[k] = ib;
        (*id_C_s)[k] = ic;
    } 
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                        *** READSTRUCTURALMODES ***                         //
//                                                                            //
//============================================================================//
int CAE3DreadStructuralModes( label N_s_v, label Nmodes, Matrix<scalar> *UUx_s, Matrix<scalar> *UUy_s, Matrix<scalar> *UUz_s, 
                                                         scalarField *ff_s, scalarField *mm_s, scalarField *cc_s )
{
    // Variables definition
    label im, is;
    FILE *f1;
    char filename[40];
    scalar m_, c_, f_, ux_, uy_, uz_;
    
    // Memory allocation
    (*ff_s) = scalarField(Nmodes, 0.0);
    (*mm_s) = scalarField(Nmodes, 0.0);
    (*cc_s) = scalarField(Nmodes, 0.0);
        
    // Loop on modes
    for( im = 0; im < Nmodes; im++ )
    {
        // Open file and read mode frequency, generalized mass and damping
        sprintf(filename, "./Data/StructuralModel.mode%i", im + 1);
        f1 = fopen(filename, "r");
        if ( f1 == NULL ) Info << "ERROR: File " << filename << " not found!" << nl;
        fscanf(f1, "%lf\n", &f_);
        fscanf(f1, "%lf\n", &m_);
        fscanf(f1, "%lf\n", &c_);
        (*ff_s)[im] = f_;
        (*mm_s)[im] = m_;
        (*cc_s)[im] = c_;
        // Loop on structural nodes
        for ( is = 0; is < N_s_v; is++ )
        {
            fscanf(f1, "%lf %lf %lf\n", &ux_, &uy_, &uz_);
            (*UUx_s)[is][im] = ux_;
            (*UUy_s)[is][im] = uy_;       
            (*UUz_s)[is][im] = uz_;
        }
        fclose(f1);
    }
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** PRINTMODALSTATE ***                          //
//                                                                            //
//============================================================================//
int CAE3DprintModalState( scalarField *qq_s, scalarField *qqp_s )
{
    // Variables definition
    FILE *f1;
    
    // Print modal state for restart
    f1 = fopen("./Log/qqs.sav", "w");
    forAll( (*qq_s), ii ) fprintf(f1, "%g\n", (*qq_s)[ii] );
    forAll( (*qqp_s), ii ) fprintf(f1, "%g\n", (*qqp_s)[ii] );
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** READMODALSTATE ***                           //
//                                                                            //
//============================================================================//
int CAE3DreadModalState( scalarField *qq_s, scalarField *qqp_s )
{
    // Variables definition
    scalar tmp;
    FILE *f1;
    
    // Print modal state for restart
    f1 = fopen("./Log/qqs.sav", "r");
    forAll( (*qq_s), ii ) 
    { 
        fscanf(f1, "%lf\n", &tmp );
        (*qq_s)[ii] = tmp;
    }    
    forAll( (*qqp_s), ii ) 
    { 
        fscanf(f1, "%lf\n", &tmp );
        (*qqp_s)[ii] = tmp;
    }    
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** PRINTNAEMOINF ***                           //
//                                                                            //
//============================================================================//
int CAE3DprintNAEMOInf( label NactiveMode, scalar qMax_s, scalar Uoo_a, scalar rhooo_a, scalar Lref_a, scalar Sref_a, scalar dt )
{
    // Variables definition
    FILE *f1;
    char filename[40];
    
    // Build filename
    // Warning: works only for 9 modes or less!
    sprintf(filename, "./Log/NAEMOGenForces_m%i.inf", NactiveMode);
    
    // Create NAEMO input file NAEMOGenForces_m<Nmode>.inf
    f1 = fopen(filename, "w");
    fprintf(f1, "Applied_scale_factor: %g\n", qMax_s);
    fprintf(f1, "Used_time_step: %g\n", dt);
    fprintf(f1, "Reference_velocity: %g\n", Uoo_a);
    fprintf(f1, "Reference_length: %g\n", Lref_a); // Lref = 0.5*Cref  
    fprintf(f1, "Reference_density: %g\n", rhooo_a);
    fprintf(f1, "Reference_surface: %g\n", Sref_a);
    fprintf(f1, "Iterations: !!!Copy here the number of lines in NAEMOGenForces_m%i.Wrk file!!!\n", NactiveMode);
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** PRINTNAEMOWRK ***                           //
//                                                                            //
//============================================================================//
int CAE3DprintNAEMOWrk( label NactiveMode, scalar qoo_a, scalarField *qq_s, scalarField *QQ0_s, scalarField *QQ_s, label Niter )
{
    // Variables definition
    FILE *f1;
    char filename[40];
    
    // Build filename
    // Warning: works only for 9 modes or less!
    sprintf(filename, "./Log/NAEMOGenForces_m%i.wrk", NactiveMode);
    
    // Create NAEMO input file NAEMOGenForces_m<Nmode>.wrk
    if( Niter <= 1 )
    {
        // If this is the first iteration the existent .wrk file should be erased
        f1 = fopen(filename, "w");
    }
    else
    {
        f1 = fopen(filename, "a");
    }
    fprintf(f1, "%g  ", (*qq_s)[NactiveMode-1]);
    forAll( (*QQ_s), im )
    {
        fprintf(f1, " %g", ( (*QQ_s)[im] - (*QQ0_s)[im] )/qoo_a );
    }
    fprintf(f1, "\n");
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                          *** PRINTDEFORMEDSHAPE ***                        //
//                                                                            //
//============================================================================//
int CAE3DprintDeformedShape( vectorField *xx_s_v, vectorField *uu_s_v, labelField *id_A_s, labelField *id_B_s, labelField *id_C_s )
{
     // Variables definition
    label iA_s, iB_s, iC_s;
    scalar color;
    vector xA_s, xB_s, xC_s;
    vector uA_s, uB_s, uC_s;
    FILE *f1;
    
    // Open files
    f1 = fopen("./Log/StructuralModel.deform", "w");
    
    // Loop on structural elements
    forAll( (*id_A_s), ie )
    {
        iA_s = (*id_A_s)[ie];
        iB_s = (*id_B_s)[ie];
        iC_s = (*id_C_s)[ie];
        xA_s = (*xx_s_v)[iA_s-1];
        xB_s = (*xx_s_v)[iB_s-1];
        xC_s = (*xx_s_v)[iC_s-1]; 
        uA_s = (*uu_s_v)[iA_s-1];
        uB_s = (*uu_s_v)[iB_s-1];
        uC_s = (*uu_s_v)[iC_s-1]; 
        
        // GNUplot format deformed with color
        color = ( uA_s.z() + uB_s.z() + uC_s.z() )/3.0;
        fprintf(f1, "%g %g %g %g\n", xA_s.x() + uA_s.x(), xA_s.y() + uA_s.y(), xA_s.z() + uA_s.z(), color );
        fprintf(f1, "%g %g %g %g\n", xB_s.x() + uB_s.x(), xB_s.y() + uB_s.y(), xB_s.z() + uB_s.z(), color );
        fprintf(f1, "\n");
        fprintf(f1, "%g %g %g %g\n", xC_s.x() + uC_s.x(), xC_s.y() + uC_s.y(), xC_s.z() + uC_s.z(), color );
        fprintf(f1, "%g %g %g %g\n", xC_s.x() + uC_s.x(), xC_s.y() + uC_s.y(), xC_s.z() + uC_s.z(), color );
        fprintf(f1, "\n\n");
    }
    
    // Close files
    fclose(f1);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                             *** PRINTLOADS ***                             //
//                                                                            //
//============================================================================//
int CAE3DprintLoads( scalar t, scalar C_FX, scalar C_FY, scalar C_FZ, scalar C_MX, scalar C_MY, scalar C_MZ )
{
    // Variables definition
    FILE *f1;

    // Print loads on file ( compatible with geplot.out/gnuplot  )
    f1 = fopen("./Log/AerodynamicLoads.txt", "a");
    fprintf(f1, "%g %g %g %g %g %g %g\n", t, C_FX, C_FY, C_FZ, C_MX, C_MY, C_MZ );
    fclose(f1); 
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                               *** PRINTCP ***                              //
//                                                                            //
//============================================================================//
int CAE3DprintCp( vectorField *xx, scalarField *CCp )
{
    // Variables definition
    FILE *f1;

    // Print loads on file ( compatible with geplot.out/gnuplot  )
    f1 = fopen("./Log/CCp.txt", "w");
    forAll( (*CCp), ii )
    {
        fprintf(f1, "%g %g %g %g\n", (*xx)[ii].x(), (*xx)[ii].y(), (*xx)[ii].z(), (*CCp)[ii] );
    }
    fclose(f1); 
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                      *** PRINTGENDISPLACEMENT ***                          //
//                                                                            //
//============================================================================//
int CAE3DprintGenDisplacement( scalar t, scalar Uoo_a, scalar Lref_a, scalarField *qq_s )
{
    // Variables definition
    FILE *f1;
    
    // Print generalized forces on file
    f1 = fopen("./Log/qqs.txt", "a");
    fprintf(f1, "%g ", t*Uoo_a/Lref_a );
    forAll( (*qq_s), im ) fprintf(f1, "%g ", (*qq_s)[im] );
    fprintf(f1, "\n");
    fclose(f1); 
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                          *** PRINTGENFORCES ***                            //
//                                                                            //
//============================================================================//
int CAE3DprintGenForces( scalar t, scalar Uoo_a, scalar Lref_a, scalar qoo_a, scalarField *QQ0_s, scalarField *QQ_s )
{
    // Variables definition
    FILE *f1;
    
    // Print generalized forces on file
    f1 = fopen("./Log/QQa.txt", "a");
    fprintf(f1, "%g ", t*Uoo_a/Lref_a );
    forAll( (*QQ_s), im ) fprintf(f1, "%g ", ( (*QQ_s)[im] - (*QQ0_s)[im] )/qoo_a );
    fprintf(f1, "\n");
    fclose(f1); 
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                             *** ISPINTRI ***                               //
//                                                                            //
//============================================================================//
int CAE3DisPinTri( vector P, vector xA, vector xB, vector xC )
{
    // Variables definition
    label isin;
    vector nnAB, nnBC, nnCA;
    
    // Compute normal vectors (non-normalized)
    nnAB = ( xA - P ) ^ ( xB - P );
    nnBC = ( xB - P ) ^ ( xC - P );
    nnCA = ( xC - P ) ^ ( xA - P );
    // If all the normal vectors point in the same direction, P lays inside the triangle
    // Check for errors if P lays on the triangle boundary
    isin = 0;
    if ( ( nnAB & nnBC ) >= 0.0 )
    {
        if ( ( nnAB & nnCA ) >= 0.0 )
        {
            isin = 1;
        }
    }

    // Return
    return(isin);
}

//============================================================================//
//                                                                            //
//                            *** TRINORMAL ***                               //
//                                                                            //
//============================================================================//
int CAE3DtriNormal( vector xA, vector xB, vector xC, vector *nn )
{   
    // Variables definition
    vector S_hat;

    // Compute normal vector (normalized so that |nn| = 1)
    S_hat = ( xB - xA ) ^ ( xC - xA );
    (*nn) = S_hat/mag(S_hat);
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                             *** BUILDA2S ***                               //
//                                                                            //
//============================================================================//
int CAE3DbuildA2S( fvMesh *mesh, label ia, vectorField *xx_s_v, labelField *id_A_s, labelField *id_B_s, labelField *id_C_s  )
{
    // Variables definition
    label iA_s, iB_s, iC_s, is;
    vector xA_s, xB_s, xC_s, nn;
    vector CG_a, CG_a_s;
    
    // *** Compute Aerodynamic To Structure connectivity matrix
    // Hp: - Structural model is a plane surface with constant normal vector
    //     - Structural mesh is triangular
    //     - Works only for an isolated wing (MLS should be used for complex geometries)
    forAll( (*id_A_s), k ) 
    {
         // Find S element vertices
         iA_s = (*id_A_s)[k];
         iB_s = (*id_B_s)[k];
         iC_s = (*id_C_s)[k];
         xA_s = (*xx_s_v)[iA_s-1];
         xB_s = (*xx_s_v)[iB_s-1];
         xC_s = (*xx_s_v)[iC_s-1]; 
         
         // Compute S mesh normal
         nn = ( xA_s - xC_s ) ^ ( xB_s - xC_s );
         nn = nn/mag(nn);
          
         // Project A node onto S mesh   
         CG_a = (*mesh).Cf()[ia];
         CG_a_s = CG_a - ( CG_a & nn )*nn;        
          
         // Find parent S element 
         is = CAE3DisPinTri( CG_a_s, xA_s, xB_s, xC_s );
         if ( is == 1 ) return(k);
    }
    
    // Return
    return(-1);
} 

//============================================================================//
//                                                                            //
//                             *** BUILDTAS ***                               //
//                                                                            //
//============================================================================//
int CAE3DbuildTas( labelField *id_A_s, labelField *id_B_s, labelField *id_C_s, labelField *id_a2s, 
                   vectorField *xx_s_v, vectorField *xx_a, Matrix<scalar> *Tas )
{
    // Variables definition
    label is, iA_s, iB_s, iC_s;
    scalar det, aA, aB, aC, bA, bB, bC, cA, cB, cC;
    scalar xA_s, yA_s, xB_s, yB_s, xC_s, yC_s;

    // Memory allocation
    (*Tas)  = 0.0;
    forAll( (*xx_a), ia )
    {
        is   = (*id_a2s)[ia];
        iA_s = (*id_A_s)[is];
        iB_s = (*id_B_s)[is];
        iC_s = (*id_C_s)[is];
        xA_s = (*xx_s_v)[iA_s-1].x();
        xB_s = (*xx_s_v)[iB_s-1].x();
        xC_s = (*xx_s_v)[iC_s-1].x();
        yA_s = (*xx_s_v)[iA_s-1].y();
        yB_s = (*xx_s_v)[iB_s-1].y();
        yC_s = (*xx_s_v)[iC_s-1].y();
        det  = ( xB_s*yC_s - xC_s*yB_s ) - xA_s*( yC_s - yB_s ) + yA_s*( xC_s - xB_s );
        aA   =   ( xB_s*yC_s - yB_s*xC_s )/det;
        aB   = - ( xA_s*yC_s - yA_s*xC_s )/det;
        aC   =   ( xA_s*yB_s - yA_s*xB_s )/det;
        bA   = - ( yC_s - yB_s )/det;
        bB   =   ( yC_s - yA_s )/det;
        bC   = - ( yB_s - yA_s )/det;
        cA   =   ( xC_s - xB_s )/det;
        cB   = - ( xC_s - xA_s )/det;
        cC   =   ( xB_s - xA_s )/det;
        (*Tas)[ia][iA_s-1] = aA + bA*(*xx_a)[ia].x() + cA*(*xx_a)[ia].y();
        (*Tas)[ia][iB_s-1] = aB + bB*(*xx_a)[ia].x() + cB*(*xx_a)[ia].y();
        (*Tas)[ia][iC_s-1] = aC + bC*(*xx_a)[ia].x() + cC*(*xx_a)[ia].y();        
    }

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** COMPUTE_CP ***                               //
//                                                                            //
//============================================================================//
int CAE3DcomputeCp( volScalarField *P, scalar Poo, scalar qoo, 
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
//                           *** COMPUTE_LOADS ***                            //
//                                                                            //
//============================================================================//
int CAE3DcomputeLoads( vectorField *xx, scalarField *SS, vectorField *nn, scalarField *CCp, 
                       vector i_wind, vector j_wind, vector k_wind, scalar Sref, scalar Lref, vector Xref, 
                       scalar *C_L, scalar *C_D, scalar *C_X, scalar *C_Y, scalar *C_Z, scalar *C_MX, scalar *C_MY, scalar *C_MZ )
{
    // Variables definition
    vector i_hat, j_hat, k_hat;
    i_hat = vector( 1.0, 0.0, 0.0 );
    j_hat = vector( 0.0, 1.0, 0.0 );
    k_hat = vector( 0.0, 0.0, 1.0 );
    
    // Compute Lift, Drag, Force and Moment coefficients
    (*C_L)  = 0.0; 
    (*C_D)  = 0.0;
    (*C_X)  = 0.0;
    (*C_Y)  = 0.0;
    (*C_Z)  = 0.0;
    (*C_MX) = 0.0;
    (*C_MY) = 0.0;
    (*C_MZ) = 0.0;
    forAll( (*SS), ii )
    {
        // REMARK: nn0, nn vectors point outside the body
        (*C_L)  = (*C_L)  - (*SS)[ii]*( (*nn)[ii] & k_wind )*(*CCp)[ii];
        (*C_D)  = (*C_D)  - (*SS)[ii]*( (*nn)[ii] & i_wind )*(*CCp)[ii]; 
        (*C_X)  = (*C_X)  - (*SS)[ii]*( (*nn)[ii] & i_hat  )*(*CCp)[ii];
        (*C_Y)  = (*C_Y)  - (*SS)[ii]*( (*nn)[ii] & j_hat  )*(*CCp)[ii];
        (*C_Z)  = (*C_Z)  - (*SS)[ii]*( (*nn)[ii] & k_hat  )*(*CCp)[ii];
        // REMARK: in the 3D case C_MP refers to the aerodynamic moment around the P-Axis of the body frame of reference
        (*C_MX) = (*C_MX) - (*SS)[ii]*( ( ( (*xx)[ii] - Xref ) ^ (*nn)[ii]*(*CCp)[ii] ) & i_hat ); 
        (*C_MY) = (*C_MY) - (*SS)[ii]*( ( ( (*xx)[ii] - Xref ) ^ (*nn)[ii]*(*CCp)[ii] ) & j_hat ); 
        (*C_MZ) = (*C_MZ) - (*SS)[ii]*( ( ( (*xx)[ii] - Xref ) ^ (*nn)[ii]*(*CCp)[ii] ) & k_hat );      
    }
    (*C_L)  = (*C_L)/Sref;
    (*C_D)  = (*C_D)/Sref;
    (*C_X)  = (*C_X)/Sref;
    (*C_Y)  = (*C_Y)/Sref;
    (*C_Z)  = (*C_Z)/Sref;
    (*C_MX) = (*C_MX)/( Sref*Lref );  
    (*C_MY) = (*C_MY)/( Sref*Lref );  
    (*C_MZ) = (*C_MZ)/( Sref*Lref );    

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** COMPUTE_R_A ***                             //
//                                                                            //
//============================================================================//
int CAE3DcomputeR_a( scalar qoo_a, scalarField *SS_a, vectorField *nn_a, scalarField *CCp_a, vectorField *RR_a )
{   
    // Loop on aerodynamic elements
    forAll( (*RR_a), ia )
    {
        (*RR_a)[ia] = - qoo_a*(*CCp_a)[ia]*(*SS_a)[ia]*(*nn_a)[ia];
    }
   
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** COMPUTE_R_S_V ***                            //
//                                                                            //
//============================================================================//
int CAE3DcomputeR_s_v( Matrix<scalar> *Tas, vectorField *RR_a, vectorField *RR_s_v )
{   
    // Variables definition
    label ia, is, Na, Ns;
    Na = (*Tas).n(); // Warning!
    Ns = (*Tas).m(); // Warning!

    // Warning: since Tsa = Tas^T loops order should be inverted
    // Loop on Structural vertices
    for( is = 0; is < Ns; is++ )
    {
        // Loop on Aerodynamic elements
        (*RR_s_v)[is] = vector( 0.0, 0.0, 0.0 );
        for( ia = 0; ia < Na; ia++ )
        {
            (*RR_s_v)[is] = (*RR_s_v)[is] + (*Tas)[ia][is]*(*RR_a)[ia];
        }    
    }    
        
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                             *** COMPUTE_Q_S ***                            //
//                                                                            //
//============================================================================//
int CAE3DcomputeQ_s( Matrix<scalar> *UUx_s, Matrix<scalar> *UUy_s, Matrix<scalar> *UUz_s, vectorField *RR_s_v, scalarField *QQ_s )
{
    // Variables definition
    label im, is, Nm, Ns;
    Ns = (*UUx_s).n(); // Warning!
    Nm = (*UUx_s).m(); // Warning!    
    
    // Loops are inverted in order not to store [U]' also
    // Loop on structural modes 
    for ( im = 0; im < Nm; im++ )
    {
        // Loop on structural vertices
        (*QQ_s)[im] = 0.0;
        for ( is = 0; is < Ns; is++ )
        {
            (*QQ_s)[im] = (*QQ_s)[im] + (*UUx_s)[is][im]*(*RR_s_v)[is].x()
                                      + (*UUy_s)[is][im]*(*RR_s_v)[is].y()
                                      + (*UUz_s)[is][im]*(*RR_s_v)[is].z();
        }
    }
        
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                          *** COMPUTE_NN_S_E ***                            //
//                                                                            //
//============================================================================//
int CAE3DcomputeNn_s_e( vectorField *xx_s_v, vectorField *uu_s_v, labelField *id_A_s, labelField *id_B_s, labelField *id_C_s, 
                         vectorField *nn_s_e )
{
    // Variables definition
    label iA_s, iB_s, iC_s;
    vector xA_s, xB_s, xC_s;
    vector uA_s, uB_s, uC_s;
    vector nn;
    
    // Loop on structural elements
    forAll( (*id_A_s), ie )
    {
        iA_s = (*id_A_s)[ie];
        iB_s = (*id_B_s)[ie];
        iC_s = (*id_C_s)[ie];
        xA_s = (*xx_s_v)[iA_s-1];
        xB_s = (*xx_s_v)[iB_s-1];
        xC_s = (*xx_s_v)[iC_s-1]; 
        uA_s = (*uu_s_v)[iA_s-1];
        uB_s = (*uu_s_v)[iB_s-1];
        uC_s = (*uu_s_v)[iC_s-1]; 
        CAE3DtriNormal( xA_s + uA_s, xB_s + uB_s, xC_s + uC_s, &nn );
        (*nn_s_e)[ie] = nn;
    }

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** COMPUTE_NN_A ***                            //
//                                                                            //
//============================================================================//
int CAE3DcomputeNn_a( labelField *id_a2s, labelField *id_A_s, labelField *id_B_s, labelField *id_C_s, vectorField *xx_a, vectorField *xx_s_v,
                      vectorField *nn0_s_e, vectorField *nn_s_e, vectorField *nn0_a, vectorField *nn_a )
{
    // Variables definition
    label is, iA_s, iB_s, iC_s;
    vector xA_s, xB_s, xC_s, CG_s;
    
    // Loop on aerodynamic elements
    forAll( (*nn_a), ia )
    {    
        // Find the parent structural element
        is = (*id_a2s)[ia];
        iA_s = (*id_A_s)[is];
        iB_s = (*id_B_s)[is];
        iC_s = (*id_C_s)[is];
        xA_s = (*xx_s_v)[iA_s-1];
        xB_s = (*xx_s_v)[iB_s-1];
        xC_s = (*xx_s_v)[iC_s-1]; 
        CG_s = ( xA_s + xB_s + xC_s )/3.0;
        
        // Distinguish upper and lower surfaces
        if ( ( ( (*xx_a)[ia] - CG_s ) & (*nn0_s_e)[is] ) >= 0 )
        {
            (*nn_a)[ia] = (*nn0_a)[ia] + ( (*nn_s_e)[is] - (*nn0_s_e)[is] );
        }
        else
        {
            (*nn_a)[ia] = (*nn0_a)[ia] - ( (*nn_s_e)[is] - (*nn0_s_e)[is] );
        }
        // At wing tip do not deform normal vectors
        if ( mag( ( (*nn0_a)[ia] & (*nn0_s_e)[is] ) ) < 1e-5 ) (*nn_a)[ia] = (*nn0_a)[ia];
        // Normalization
        (*nn_a)[ia] = (*nn_a)[ia]/mag( (*nn_a)[ia] );
    }
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** COMPUTE_UU_S_V ***                           //
//                                                                            //
//============================================================================//
int CAE3DcomputeUu_s_v( Matrix<scalar> *UUx_s, Matrix<scalar> *UUy_s, Matrix<scalar> *UUz_s, scalarField *qq_s, vectorField *uu_s_v )
{
    // Variables definition
    label im, is, Nm, Ns;
    Ns = (*UUx_s).n(); // Warning!
    Nm = (*UUx_s).m(); // Warning!
    
    // Loop on structural vertices
    for ( is = 0; is < Ns; is++ )
    {
        // Loop on structural modes
        (*uu_s_v)[is] = vector(0.0, 0.0, 0.0);
        for ( im = 0; im < Nm; im++ )
        {
            (*uu_s_v)[is].x() = (*uu_s_v)[is].x() + (*UUx_s)[is][im]*(*qq_s)[im];
            (*uu_s_v)[is].y() = (*uu_s_v)[is].y() + (*UUy_s)[is][im]*(*qq_s)[im];
            (*uu_s_v)[is].z() = (*uu_s_v)[is].z() + (*UUz_s)[is][im]*(*qq_s)[im];
        }
    }
        
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** COMPUTE_UU_A ***                            //
//                                                                            //
//============================================================================//
int CAE3DcomputeUu_a( Matrix<scalar> *Tas, vectorField *uu_s_v, vectorField *uu_a )
{
    // Variables definition
    label ia, is, Na, Ns;
    Na = (*Tas).n(); // Warning!
    Ns = (*Tas).m(); // Warning!
    
    // Loop on Aerodynamic elements
    for( ia = 0; ia < Na; ia++ )
    {
        // Loop on Structural vertices
        (*uu_a)[ia] = vector( 0.0, 0.0, 0.0 );
        for( is = 0; is < Ns; is++ )
        {
            (*uu_a)[ia] = (*uu_a)[ia] + (*Tas)[ia][is]*(*uu_s_v)[is];
        }    
    }    

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                        *** FORCEDMODALMOTIONPARAM ***                      //
//                                                                            //
//============================================================================//
int CAE3DforcedModalMotionParam( Matrix<scalar> *UUx_s, Matrix<scalar> *UUy_s, Matrix<scalar> *UUz_s, 
                                 label NactiveMode, scalar kMax, scalar Uoo_a, scalar Lref_a, scalar epsilon, 
                                 scalar *UUMax, scalar *tauMax, scalar *qMax )
{
    // Variables definition
    label im, is, Nm, Ns;
    scalar pi = 3.14159265358979;
    scalar Utmp, UMax;
    scalar Omega0;
        
    // Dimensions
    Nm = (*UUx_s).m(); // Warning!
    Ns = (*UUx_s).n(); // Warning!    
      
    // UMax    
    UMax = 0.0;
    im   = NactiveMode-1;
    for( is = 1; is < Ns; is++ )
    {
        Utmp = Foam::sqrt( Foam::pow( (*UUx_s)[is][im], 2) + 
                           Foam::pow( (*UUy_s)[is][im], 2) + 
                           Foam::pow( (*UUz_s)[is][im], 2) );
        if ( Utmp > UMax ) UMax = Utmp;                   
    }
    (*UUMax) = UMax;  

    // tauMax 
    (*tauMax) = 2*pi/kMax;

    // qMax (Linearity should be verified) 
    Omega0 = pi/(*tauMax);
    (*qMax) = 2.0*epsilon*Lref_a/( UMax*Omega0 );

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                           *** FORCEDMODALMOTION ***                        //
//                                                                            //
//============================================================================//
int CAE3DforcedModalMotion( label NactiveMode, scalar qMax_s, scalar tauMax_s, scalar Uoo_a, scalar Lref_a, scalar t, scalar t0,
                            scalarField *qq_s, scalarField *qqp_s )
{
    // Variables definition
    scalar pi = 3.14159265358979;
    scalar Omega0_s, tau_s, q_t, qp_t;
    
    // Smoothed step forced modal motion law
    Omega0_s = pi/tauMax_s;
    tau_s    = Uoo_a/Lref_a*( t - t0 );
    if ( ( tau_s > 0.0 ) && ( tau_s <= tauMax_s ) )
    {
        q_t  = 0.5*qMax_s*( 1.0 - Foam::cos( Omega0_s*tau_s ) );
        qp_t = 0.5*qMax_s*Omega0_s*Foam::sin( Omega0_s*tau_s )*Uoo_a/Lref_a;
    }
    else if ( tau_s > tauMax_s )
    {
        q_t  = qMax_s;
        qp_t = 0.0;
    }
    else
    {
        q_t  = 0.0;
        qp_t = 0.0;
    }
    
    // Initialize active mode
    (*qq_s)[NactiveMode-1]  = q_t;
    (*qqp_s)[NactiveMode-1] = qp_t;   

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                 *** COMPUTE_TRANSPIRATION_VELOCITY ***                     //
//                                                                            //
//============================================================================//
int CAE3DcomputeTranspirationVelocity( vectorField *uup_a, vectorField *nn0_a, vectorField *nn_a, volVectorField *U, labelField *id_bodyCell, 
                                       scalarField *VVbn_a )
{
    // Variables definition
    label id_L;
    
    // Loop on control points (boundary faces centres)
    forAll( (*VVbn_a), ii )
    {
        id_L          = (*id_bodyCell)[ii];      
        // REMARK: nn0, nn vectors point outside the body
        // Warning: in OpenFOAM nn points outside the computational domain
        (*VVbn_a)[ii] = -( (*U)[id_L] & ( (*nn_a)[ii] - (*nn0_a)[ii] ) ) + ( (*uup_a)[ii] & (*nn_a)[ii] ); 
    }
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                 *** SET_TRANSPIRATION_VELOCITY_BCS ***                     //
//                                                                            //
//============================================================================//
int CAE3DsetTranspirationVelocityBCs( fvMesh *mesh, label id_bodyPatch, volVectorField *U, scalarField *VVbn_a, vectorField *nn_a )
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
            (*U).boundaryField()[id_bodyPatch][ii].x() = -(*VVbn_a)[ii];
            // 2) Modified transpiration boundary condition (with deformed noraml vector)
            //(*U).boundaryField()[id_bodyPatch][ii] = (*VVbn_a)[ii]*(*nn_a)[ii]; 
        } 
    }
     
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                        *** FORCED_RIGID_MOTION ***                         //
//                                                                            //
//============================================================================//
int CAE3DforcedRigidMotion( label flag, scalar A0, scalar A1, scalar f, scalar tau, scalar t, 
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
    
    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                         *** COMPUTE_RIGID_NN_A ***                         //
//                                                                            //
//============================================================================//
int CAE3DcomputeRigidNn_a( scalar h, scalar a, scalar h_p, scalar a_p, vector h_hat, vector a_hat,
                           vectorField *xx_a, vector Xref_a, vectorField *nn0_a, vectorField *nn_a )
{
    // Normalization
    a_hat = a_hat/mag( a_hat );
    
    // Loop on aerodynamic elements
    forAll( (*xx_a), ii )
    {
        // Rotation around a_hat axis
        (*nn_a)[ii] = (*nn0_a)[ii] 
                    + Foam::sin( a )*( a_hat ^ (*nn0_a)[ii] ) 
                    + ( 1.0 - Foam::cos( a ) )*( a_hat ^ ( a_hat ^ (*nn0_a)[ii] ) );
    }

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                        *** COMPUTE_RIGID_UUP_A ***                         //
//                                                                            //
//============================================================================//
int CAE3DcomputeRigidUup_a( scalar h, scalar a, scalar h_p, scalar a_p, vector h_hat, vector a_hat,
                            vectorField *xx_a, vector Xref_a, vectorField *uup_a )
{
    // Variables definition
    vector Duup_a;
    
    // Loop on aerodynamic elements
    forAll( (*xx_a), ii )
    {
        Duup_a = ( a_p*a_hat )^( (*xx_a)[ii] - Xref_a ) + h_p*h_hat;
        (*uup_a)[ii] = Duup_a;
    }

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                        *** STATICMODALSOLUTION ***                         //
//                                                                            //
//============================================================================//
int CAE3DstaticModalSolution( scalarField *mm_s, scalarField *cc_s, scalarField *ff_s, scalarField *QQ_s, 
                              scalarField *qq_s, scalarField *qqp_s )
{
    // Variables definition
    scalar pi = 3.14159265358979;
    scalar ki;
    
    // Modal static solution
    forAll( (*qq_s), i )
    {
        ki = ( (*mm_s)[i]*( 2*pi*(*ff_s)[i] )*( 2*pi*(*ff_s)[i] ) );
        (*qq_s)[i] = (*QQ_s)[i]/ki;
    }
    (*qqp_s) = 0.0;

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                        *** DYNAMICMODALSOLUTION ***                        //
//                                                                            //
//============================================================================//
int CAE3DdynamicModalSolution( scalarField *mm_s, scalarField *cc_s, scalarField *ff_s, scalarField *QQ_s, 
                               scalar dt, scalarField *qq_s, scalarField *qqp_s )
{
    // Variables definition
    label N_s;
    scalar pi = 3.14159265358979;
    scalar mi, ci, ki;
    scalarField qqo_s, qqpo_s;
    
    // Dimensions
    N_s = (*qq_s).size();
    
    // Save old solution
    qqo_s  = scalarField(N_s, 0.0);
    qqpo_s = scalarField(N_s, 0.0);
    forAll( (*qq_s), i )
    {
        qqo_s[i]  = (*qq_s)[i];
        qqpo_s[i] = (*qqp_s)[i];
    }
    
    // Modal dynamic solution (EE)
    forAll( (*qq_s), i )
    {
        // 1:N_s state space equations
        (*qq_s)[i] = qqo_s[i] + qqpo_s[i]*dt;
               
        // N_s+1:2*N_s state space equations
        mi = (*mm_s)[i];
        ci = (*cc_s)[i];
        ki = ( (*mm_s)[i]*( 2*pi*(*ff_s)[i] )*( 2*pi*(*ff_s)[i] ) );
        (*qqp_s)[i] = qqpo_s[i] + ( -ki/mi*qqo_s[i] - ci/mi*qqpo_s[i] + 1.0/mi*(*QQ_s)[i] )*dt;
    }

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                            *** PRINTRIGIDMODE ***                          //
//                                                                            //
//============================================================================//
int CAE3DprintRigidMode( scalar h, scalar a, scalar h_p, scalar a_p, vector h_hat, vector a_hat,
                         vectorField *xx_s_v, vector Xref_a )
{
    // Variables definition
    vector xx_para, xx_orto;
    vector UUloc;
    FILE *fid;

    // Normalization
    h_hat = h_hat/mag( h_hat );
    a_hat = a_hat/mag( a_hat );
       
    // Loop on structural nodes
    fid = fopen("./Log/StructuralModel.mode0", "w");
    fprintf(fid, "0.0\n");
    fprintf(fid, "0.0\n");
    fprintf(fid, "0.0\n");
    forAll( (*xx_s_v), ii )
    {
        // Translation along h_hat and Rotation around a_hat axis
        xx_para = ( ( (*xx_s_v)[ii] - Xref_a ) & a_hat )*a_hat;
        xx_orto =   ( (*xx_s_v)[ii] - Xref_a ) - xx_para;
        UUloc = h*h_hat
              + Foam::sin( a )*( a_hat ^ xx_orto ) 
              + ( 1.0 - Foam::cos( a ) )*( a_hat ^ ( a_hat ^ xx_orto ) ); // TODO: verify this 
        fprintf(fid, "%g %g %g\n", UUloc.x(), UUloc.y(), UUloc.z() );   
    }
    fclose(fid);

    // Return
    return(0);
}

//============================================================================//
//                                                                            //
//                       *** UPDATECELLDISPLACEMENT ***                       //
//                                                                            //
//============================================================================//
int CAE3DupdateCellDisplacement( label id_bodyPatch, volVectorField *cellDisplacement, vectorField *uu_a )
{   
    // Update
    forAll( (*uu_a), ii )
    {
        (*cellDisplacement).boundaryField()[id_bodyPatch][ii] = (*uu_a)[ii];    
    }
    
    // Return
    return(0);
}
