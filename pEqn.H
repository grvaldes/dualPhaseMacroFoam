{
    p.storeOldTime();
    xi.storeOldTime();

    Info<< "P avant = " << average(p) << endl;

    fvScalarMatrix pEqn
    (
        fvm::laplacian(-rhogf*Mgf-rhol*Mlf,p) + 
        nf*(rhol - rhog)*fvc::ddt(xi) + 
        nf*(1 - xi)*kg*fvm::ddt(p)
    );

    pEqn.solve();

    Info<< "P apres = " << average(p) << endl;

    surfaceScalarField phiP = pEqn.flux();

    phig = phiP * Mgf / (rhogf*Mgf + rhol*Mlf);
    phil = phiP * Mlf / (rhogf*Mgf + rhol*Mlf);

    Info<< "phiP = " << average(phiP) << nl
        << "phig = " << average(phig) << nl
        << "phil = " << average(phil) << endl;

    Ug = fvc::reconstruct(phig);
    Ul = fvc::reconstruct(phil);

    Ug.correctBoundaryConditions();
    Ul.correctBoundaryConditions();

    CPerm.correctDensity();
}
