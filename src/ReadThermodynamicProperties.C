//==============================================================================
// *** READ THERMODYNAMIC PROPERTIES ***
//==============================================================================
Info << "========================================"           << nl;
Info << " Reading Thermodynamic Properties..."               << nl;
Info << "========================================"           << nl;
IOdictionary thermodynamicProperties
(
    IOobject
    (
        "thermodynamicProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
// R 
dimensionedScalar R
(
    thermodynamicProperties.lookup("R")
);
// Cv
dimensionedScalar Cv
(
    thermodynamicProperties.lookup("Cv")
);
// Pr
dimensionedScalar Pr
(
    thermodynamicProperties.lookup("Pr")
);
    
// gamma
dimensionedScalar gamma = ( 1 + R/Cv );
Info << " gamma  [-] = " << gamma.value()                    << nl;
Info << " Pr     [-] = " << Pr.value()                       << nl;
Info << " R  [J/KgK] = " << R.value()                        << nl;
Info << " Cp [J/KgK] = " << Cv.value()                       << nl;
Info << " Cv [J/KgK] = " << Cv.value() + R.value()           << nl;
Info << "----------------------------------------"           << nl << nl;

