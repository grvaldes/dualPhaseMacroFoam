const volTensorField& Mg = CPerm.Mg_;
const volTensorField& Ml = CPerm.Ml_;

surfaceScalarField rhogf("rhogf",fvc::interpolate(rhog,"Mg"));
surfaceTensorField Mgf("Mgf",fvc::interpolate(Mg,"Mg"));
surfaceTensorField Mlf("Mlf",fvc::interpolate(Ml,"Ml"));