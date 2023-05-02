#include "relPerm.H"

Foam::relPerm::relPerm(
	const word& name,
	const dictionary& dict,
	const volScalarField& xi,
	const volScalarField& p,
	const volScalarField& K)
	:
	name_(name), dict_(dict), xi_(xi), p_(p), K_(K),
	nf_
    (
        dimensionedScalar
        (
            "nf",
            dimless,
            dict_.lookupOrDefault("nf",1)
        )
    ),
	nb_(dict_.lookup("nb")),
	mug_(dict_.lookup("mug")),
	mul_(dict_.lookup("mul")),
	krg_
    (
        IOobject
        (
            "krg",
            xi_.time().timeName(),
            xi_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        xi_.mesh(),
        dimensionedScalar("krg",dimless,0)
    ),
	krl_
    (
        IOobject
        (
            "krl",
            xi_.time().timeName(),
            xi_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        xi_.mesh(),
        dimensionedScalar("krl",dimless,0)
    ),
	Mg_
    (
        IOobject
        (
            "Mg",
            xi_.time().timeName(),
            xi_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        xi_.mesh(),
        dimensionedScalar("Mg",dimLength*dimLength/dimPressure/dimTime, Zero)
    ),
	Ml_
    (
        IOobject
        (
            "Ml",
            xi_.time().timeName(),
            xi_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        xi_.mesh(),
        dimensionedScalar("Ml",dimLength*dimLength/dimPressure/dimTime, Zero)
    ),
	rhog_
    (
        IOobject
        (
            "rhog",
            p_.time().timeName(),
            p_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("rhog",dimDensity,0)
    ),
	rhol_(
        dimensionedScalar
        (
            "rhol", 
            dimDensity, 
            dict_.lookupOrDefault("rhol",0)
        )
    ),
	kg_
    (
        dimensionedScalar
        (
            "kg", 
            dimTime*dimTime/(dimLength*dimLength), 
            dict_.lookupOrDefault("kg",0)
        )
    )

{
	correctRelPerm();
	correctDensity();
}

void Foam::relPerm::correctRelPerm()
{
    krl_ = pow(xi_,3+2/nb_);
    krg_ = (1-xi_)*(1-xi_)*(1-pow(xi_,1+2/nb_));
    
    Mg_ = K_ * krg_ / mug_;
    Ml_ = K_ * krl_ / mul_;
}

void Foam::relPerm::correctDensity()
{
    rhog_ = kg_ * p_;
}
