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

    surfaceVectorField phiP = linearInterpolate(fvc::grad(p));
    
    phig = -(Mgf & phiP) & p.mesh().Sf();
    phil = -(Mlf & phiP) & p.mesh().Sf();

    Ug = fvc::reconstruct(phig);
    Ul = fvc::reconstruct(phil);

    Ug.correctBoundaryConditions();
    Ul.correctBoundaryConditions();

    if (!incompressible)
    {
        CPerm.correctDensity();
    }
}
