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

volScalarField xi
(
    IOobject
    (
        "xi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

IOdictionary relPermProps
(
    IOobject
    (
        "relPermProps",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

relPerm CPerm("relPerm",relPermProps,xi,p,K);

const dimensionedScalar& nf = CPerm.nf_;
const volScalarField& rhog = CPerm.rhog_;
const dimensionedScalar& rhol = CPerm.rhol_;
dimensionedScalar& kg = CPerm.kg_;

if (incompressible)
{
    kg *= 0;
}


volVectorField Ug
(
    IOobject
    (
        "Ug",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Ul
(
    IOobject
    (
        "Ul",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phig
(
    IOobject
    (
        "phig",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ug) & mesh.Sf()
);

surfaceScalarField phil
(
    IOobject
    (
        "phil",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(Ul) & mesh.Sf()
);
