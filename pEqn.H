{
    p.storeOldTime();
    xi.storeOldTime();

    fvScalarMatrix pEqn
    (
        fvm::laplacian(-rhogf*Mgf-rhol*Mlf,p) + 
        nf*(rhol - rhog)*fvc::ddt(xi) + 
        nf*(1 - xi)*kg*fvm::ddt(p)
    );

    Info<< "Ml after" << endl;
    Info<< Ml << endl;
	Info<< "Mg after" << endl;
    Info<< Mg << endl;

    Info<< "p before" << endl;
    Info<< p << endl;

    pEqn.solve();

    Info<< "p after" << endl;
    Info<< p << endl;

    surfaceVectorField phiP = linearInterpolate(fvc::grad(p));

    Info<< "phiP" << endl;
    Info<< phiP << endl;

    Info<< "phil before" << endl;
    Info<< phil << endl;

    phig = -(Mgf & phiP) & p.mesh().Sf();
    phil = -(Mlf & phiP) & p.mesh().Sf();

    Info<< "phil after" << endl;
    Info<< phil << endl;

    Ug = fvc::reconstruct(phig);
    Ul = fvc::reconstruct(phil);

    Ug.correctBoundaryConditions();
    Ul.correctBoundaryConditions();

    if (!incompressible)
    {
        CPerm.correctDensity();
    }
}
