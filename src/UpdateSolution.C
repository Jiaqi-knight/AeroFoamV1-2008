//============================================================================//
//                                                                            //
//                           *** UPDATE_SOLUTION ***                          //
//                                                                            //
//============================================================================//
int update_solution( scalar scaleDt, scalarField *localDt, fvMesh *mesh, 
                     volScalarField     *rrho,  volVectorField     *mm,  volScalarField     *eet_tilde, 
                     surfaceScalarField *FFrho, surfaceVectorField *FFm, surfaceScalarField *FFet_tilde )
{

    // Variables definition
    label id_o, id_n;                  // Mesh faces, cells, geometry
    scalar dt_o, dt_n, Sf, Vf_o, Vf_n; //

    // A) Loop on internal faces [index (i)]
    forAll( (*mesh).Sf(), i )
    {
        
        // Owner and neighbour ids
        // Remember: normal vector points from Owner cell to Neigbour cell
        id_o = (*mesh).faceOwner()[i];
        id_n = (*mesh).faceNeighbour()[i];
        Sf   = (*mesh).magSf()[i];
        Vf_o = (*mesh).V()[id_o];
        Vf_n = (*mesh).V()[id_n];
        
        // Local time step
        dt_o = scaleDt*(*localDt)[id_o];
        dt_n = scaleDt*(*localDt)[id_n];
        
        // Explicit Euler
        (*rrho)[id_o]      = (*rrho)[id_o]      - dt_o*Sf/Vf_o*(*FFrho)[i];
        (*rrho)[id_n]      = (*rrho)[id_n]      + dt_n*Sf/Vf_n*(*FFrho)[i];
        (*mm)[id_o].x()    = (*mm)[id_o].x()    - dt_o*Sf/Vf_o*(*FFm)[i].x();
        (*mm)[id_o].y()    = (*mm)[id_o].y()    - dt_o*Sf/Vf_o*(*FFm)[i].y();
        (*mm)[id_o].z()    = (*mm)[id_o].z()    - dt_o*Sf/Vf_o*(*FFm)[i].z();
        (*mm)[id_n].x()    = (*mm)[id_n].x()    + dt_n*Sf/Vf_n*(*FFm)[i].x();
        (*mm)[id_n].y()    = (*mm)[id_n].y()    + dt_n*Sf/Vf_n*(*FFm)[i].y();
        (*mm)[id_n].z()    = (*mm)[id_n].z()    + dt_n*Sf/Vf_n*(*FFm)[i].z();
        (*eet_tilde)[id_o] = (*eet_tilde)[id_o] - dt_o*Sf/Vf_o*(*FFet_tilde)[i];
        (*eet_tilde)[id_n] = (*eet_tilde)[id_n] + dt_n*Sf/Vf_n*(*FFet_tilde)[i];
        
    } 

    // B) Loop on boundary faces and BCs [local index(ii), global index(i)]
    forAll( (*mesh).boundaryMesh(), iPatch )
    {        
    
        forAll( (*mesh).boundaryMesh()[iPatch].faceAreas(), ii )
        {     

            // Boundary type must not be empty
            if ( (*FFrho).boundaryField()[iPatch].size() > 0 )
            {

                // Owner id
                // Remember: normal vector points outside the computational domain
                id_o = (*mesh).boundaryMesh()[iPatch].faceCells()[ii];   
                Sf   = mag( (*mesh).boundaryMesh()[iPatch].faceAreas()[ii] );
                Vf_o = (*mesh).V()[id_o];
                
                // Local time step
                dt_o = scaleDt*(*localDt)[id_o];

                // Explicit Euler                                      
                (*rrho)[id_o]      = (*rrho)[id_o]      - dt_o*Sf/Vf_o*(*FFrho).boundaryField()[iPatch][ii];
                (*mm)[id_o].x()    = (*mm)[id_o].x()    - dt_o*Sf/Vf_o*(*FFm).boundaryField()[iPatch][ii].x();
                (*mm)[id_o].y()    = (*mm)[id_o].y()    - dt_o*Sf/Vf_o*(*FFm).boundaryField()[iPatch][ii].y();
                (*mm)[id_o].z()    = (*mm)[id_o].z()    - dt_o*Sf/Vf_o*(*FFm).boundaryField()[iPatch][ii].z();
                (*eet_tilde)[id_o] = (*eet_tilde)[id_o] - dt_o*Sf/Vf_o*(*FFet_tilde).boundaryField()[iPatch][ii];
                
            }
            
        }
    
    }  
    
    // Return
    return(0);

}
