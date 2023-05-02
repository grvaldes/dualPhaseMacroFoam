#ifndef IRELPERM_H
#define IRELPERM_H

#include "volFields.H"

namespace Foam
{

struct iRelPerm {
	word name_;
	dictionary dict_;
	const volScalarField& xi_;
	const volScalarField& K_;
	dimensionedScalar nf_;
	dimensionedScalar nb_;
	dimensionedScalar mug_;
	dimensionedScalar mul_;
	volScalarField krg_;
	volScalarField krl_;
	volScalarField Mg_;
	volScalarField Ml_;

	iRelPerm(
		const word& name,
		const dictionary& dict,
		const volScalarField& xi,
		const volScalarField& K_);

	~iRelPerm(){};

	void correctRelPerm();
};

}

#endif
