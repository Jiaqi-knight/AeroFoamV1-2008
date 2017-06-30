//----------------------------------------------------------------------------- Mesh movement 
//IF{ showDisplacement 1 
Info << " Creating displacement fields (presently only for post-processing)\n";
Info << " WARNING: CHECK THAT BCs ARE OF TYPE FIXEDVALUE, OTHERWISE DISPLACEMENT WON'T BE STORED!!!\n";
// Presently mesh movement is only performed during postprocessing, while 
// transpiration BCs ( geometric and kinematic ) are imposed
// cellDisplacement (--absolute--)
volVectorField cellDisplacement
(
    IOobject
    (
        "cellDisplacement",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// pointDisplacement (--absolute--)
// REMARK: using moveMesh utility pointDisplacement can be left unchanged
pointMesh pointMesh(mesh);
pointVectorField pointDisplacement
(
    IOobject
    (
        "pointDisplacement",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pointMesh,
    vector(vector::zero)
);
//FI} showDisplacement 1
