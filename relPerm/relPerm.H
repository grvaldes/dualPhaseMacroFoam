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
	const volTensorField& K_;
	dimensionedScalar nf_;
	dimensionedScalar nb_;
	dimensionedScalar mug_;
	dimensionedScalar mul_;
	volScalarField krg_;
	volScalarField krl_;
	volTensorField Mg_;
	volTensorField Ml_;
	volScalarField rhog_;
	dimensionedScalar rhol_;
	dimensionedScalar kg_;

	relPerm(
		const word& name,
		const dictionary& dict,
		const volScalarField& xi,
		const volScalarField& p,
		const volTensorField& K);

	~relPerm(){};

	void correctRelPerm();
	void correctDensity();
};

}

#endif
