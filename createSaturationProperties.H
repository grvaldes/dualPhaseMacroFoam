const volScalarField& Mg = CPerm.Mg_;
const volScalarField& Ml = CPerm.Ml_;

surfaceScalarField rhogf("rhogf",fvc::interpolate(rhog,"Mg"));
surfaceScalarField Mgf("Mgf",fvc::interpolate(Mg,"Mg"));
surfaceScalarField Mlf("Mlf",fvc::interpolate(Ml,"Ml"));
