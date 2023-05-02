#include "relPerm.H"

Foam::relPerm::relPerm(
	const word& name,
	const dictionary& dict,
	const volScalarField& xi,
	const volScalarField& p,
	const volScalarField& K)
	:
	name_(name), dict_(dict), xi_(xi), p_(p), K_(K),
	nf_(dict_.lookup("nf")),
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
        ),
        xi_.mesh(),
        dimensionedScalar("Ml",dimLength*dimLength/dimPressure/dimTime, Zero)
    ),
	rhog_
    (
        IOobject
        (
            "rhog",
            xi_.time().timeName(),
            xi_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        xi_.mesh(),
        dimensionedScalar("rhog", dimDensity, 0)
    ),
	rhol_(dict_.lookup("rhol")),
	kg_(dict_.lookup("kg"))
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
