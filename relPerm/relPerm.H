#ifndef RELPERM_H
#define RELPERM_H

#include "volFields.H"

namespace Foam
{

struct relPerm {
	word name_;
	dictionary dict_;
	const volScalarField& xi_;
	const volScalarField& p_;
	const volScalarField& K_;
	dimensionedScalar nf_;
	dimensionedScalar nb_;
	dimensionedScalar mug_;
	dimensionedScalar mul_;
	volScalarField krg_;
	volScalarField krl_;
	volScalarField Mg_;
	volScalarField Ml_;
	volScalarField rhog_;
	dimensionedScalar rhol_;
	dimensionedScalar kg_;

	relPerm(
		const word& name,
		const dictionary& dict,
		const volScalarField& xi,
		const volScalarField& p,
		const volScalarField& K);

	~relPerm(){};

	void correctRelPerm();
	void correctDensity();
};

}

#endif
