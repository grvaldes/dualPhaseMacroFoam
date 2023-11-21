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

    // surfaceVectorField gradP = linearInterpolate(fvc::grad(p));

    // phig = -(Mgf * gradP) & p.mesh().Sf();
    // phil = -(Mlf * gradP) & p.mesh().Sf();

    surfaceScalarField phiP = pEqn.flux();
    
    phig = phiP * Mgf / (rhogf*Mgf + rhol*Mlf);
    phil = phiP * Mlf / (rhogf*Mgf + rhol*Mlf);
    
    Ug = fvc::reconstruct(phig);
    Ul = fvc::reconstruct(phil);

    Ug.correctBoundaryConditions();
    Ul.correctBoundaryConditions();

    if (!incompressible)
    {
        CPerm.correctDensity();
    }
}
