{
    p.storeOldTime();
    xi.storeOldTime();

    fvScalarMatrix pEqn
    (
        fvm::laplacian(-rhogf*Mgf-rhol*Mlf,p) + 
        nf*(rhol - rhog)*fvc::ddt(xi) + 
        nf*(1 - xi)*kg*fvm::ddt(p)
    );

    pEqn.solve();

    surfaceScalarField phiP = pEqn.flux();
    surfaceVectorField gradP = linearInterpolate(fvc::grad(p));

    phig = (phiP - rhogf * ((Mgf & gradP) & p.mesh().Sf())) / rhogf;
    phil = (phiP - rhol * ((Mlf & gradP) & p.mesh().Sf()))/ rhol;

    Ug = fvc::reconstruct(phig);
    Ul = fvc::reconstruct(phil);

    Ug.correctBoundaryConditions();
    Ul.correctBoundaryConditions();

    if (!incompressible)
    {
        CPerm.correctDensity();
    }
}
