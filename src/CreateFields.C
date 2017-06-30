//==============================================================================
// *** CREATE FIELDS ***
//==============================================================================
Info << "========================================" << nl;
Info << " Reading Initial Conditions..."           << nl;
Info << "========================================" << nl;
//------------------------------------------------------------------------------ Primitive variables timestep t^(k)
// p
Info << " Reading field p\n";
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// T
Info << " Reading field T\n";
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// U
Info << " Reading field U\n";
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info << "----------------------------------------" << nl;
//------------------------------------------------------------------------------ Conservative variables timestep t^(k)
Info << " Creating conservative variables at t^(k)\n";
// rrho
volScalarField rrho
(
    IOobject
    (
        "rrho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p/( R*T ),
    p.boundaryField().types()
);
// mm
volVectorField mm
(
    IOobject
    (
        "mm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rrho*U,
    U.boundaryField().types()
);
// eet_tilde
volScalarField eet_tilde
(
    IOobject
    (
        "eet_tilde",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rrho*Cv*T + 0.5*rrho*magSqr(U),
    T.boundaryField().types()
);
//------------------------------------------------------------------------------ Conservative variables timestep t^(k-1)
Info << " Creating conservative variables at t^(k-1)\n";
// rrho_o
volScalarField rrho_o
(
    IOobject
    (
        "rrho_o",
        runTime.timeName(),
        mesh
    ),
    rrho,
    p.boundaryField().types()
);
// mm_o
volVectorField mm_o
(
    IOobject
    (
        "mm_o",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mm,
    U.boundaryField().types()
);
// eet_tilde_o
volScalarField eet_tilde_o
(
    IOobject
    (
        "eet_tilde_o",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eet_tilde,
    T.boundaryField().types()
);
//------------------------------------------------------------------------------ Conservative variables timestep t^(k-2)
Info << " Creating conservative variables at t^(k-2)\n";
// rrho_oo
volScalarField rrho_oo
(
    IOobject
    (
        "rrho_oo",
        runTime.timeName(),
        mesh
    ),
    rrho,
    p.boundaryField().types()
);
// mm_oo
volVectorField mm_oo
(
    IOobject
    (
        "mm_oo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mm,
    U.boundaryField().types()
);
// eet_tilde_oo
volScalarField eet_tilde_oo
(
    IOobject
    (
        "eet_tilde_oo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eet_tilde,
    T.boundaryField().types()
);
//------------------------------------------------------------------------------ Fluxes
Info << " Creating conservative variables fluxes at t^(k)\n";
// FFrho
surfaceScalarField FFrho
(
    IOobject
    (
        "FFrho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(0*rrho)
);
// FFm
surfaceVectorField FFm
(
    IOobject
    (
        "FFm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(0*mm)
);
// FFet_tilde
surfaceScalarField FFet_tilde
(
    IOobject
    (
        "FFet_tilde",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(0*eet_tilde)
);
//------------------------------------------------------------------------------ GradientFields for viscous fluxes
Info << " Creating primitive variables gradient fields (optional)\n";
// GradT
volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad( T )
);
// GradUx
volVectorField gradUx
(
    IOobject
    (
        "gradUx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad( U.component(0) )
);
// GradUy
volVectorField gradUy
(
    IOobject
    (
        "gradUy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad( U.component(1) )
);
// GradUz
volVectorField gradUz
(
    IOobject
    (
        "gradUz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad( U.component(2) )
);
//----------------------------------------------------------------------------- Deltat adaptivity
// localCo
scalarField localCo(mesh.V().size(), 0.0);
// localDt (initialized to dtsim)
scalarField localDt(mesh.V().size(), runTime.deltaT().value());
Info << "----------------------------------------" << nl << nl;
