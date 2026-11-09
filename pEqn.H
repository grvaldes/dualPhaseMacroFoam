{
    // Update pressure gradient at boundaries where Ul is prescribed with fixedGradient
    forAll(mesh.boundary(), patchi)
    {
        if (Ul.boundaryField()[patchi].fixesValue())
        {
            if (isA<fixedGradientFvPatchScalarField>(p.boundaryFieldRef()[patchi]))
            {
                fixedGradientFvPatchScalarField& pPatch =
                    refCast<fixedGradientFvPatchScalarField>(p.boundaryFieldRef()[patchi]);

                const fvsPatchScalarField& MlfPatch = Mlf.boundaryField()[patchi];
                const fvsPatchVectorField& SfPatch = mesh.Sf().boundaryField()[patchi];
                const fvsPatchScalarField& magSfPatch = mesh.magSf().boundaryField()[patchi];
                const fvPatchVectorField& UlPatch = Ul.boundaryField()[patchi];

                scalarField MlSafe(MlfPatch);
                forAll(MlSafe, facei)
                {
                    if (mag(MlSafe[facei]) < 1e-15)
                    {
                        MlSafe[facei] = 1e-15;
                    }
                }

                // Darcy's law: snGrad(p) = - (Ul & n) / Ml
                pPatch.gradient() = - (UlPatch & (SfPatch / magSfPatch)) / MlSafe;
                pPatch.evaluate();
            }
        }
    }

    fvScalarMatrix pEqn
    (
        fvm::laplacian(-rhogf*Mgf-rhol*Mlf,p) + 
        nf*(rhol - rhog)*fvc::ddt(xi) + 
        nf*(1 - xi)*kg*fvm::ddt(p)
    );
    
    // If user specified zeroGradient on p at a velocity-prescribed boundary,
    // inject the boundary mass flux directly into the matrix source
    forAll(mesh.boundary(), patchi)
    {
        if (Ul.boundaryField()[patchi].fixesValue() && !p.boundaryField()[patchi].fixesValue()
            && !isA<fixedGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            const labelList& faceCells = mesh.boundary()[patchi].faceCells();
            const vectorField& Sf = mesh.Sf().boundaryField()[patchi];
            const vectorField& Up = Ul.boundaryField()[patchi];
            forAll(faceCells, facei)
            {
                scalar phiM = rhol.value() * (Up[facei] & Sf[facei]);
                pEqn.source()[faceCells[facei]] -= phiM;
            }
        }
    }

    pEqn.solve();

    surfaceScalarField phiP = pEqn.flux();

    // On boundaries with prescribed velocity, ensure mass flux matches the prescription
    forAll(mesh.boundary(), patchi)
    {
        if (Ul.boundaryField()[patchi].fixesValue())
        {
            phiP.boundaryFieldRef()[patchi] =
                rhol.value() * (Ul.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]);
        }
    }
    
    phig = phiP * Mgf / (rhogf*Mgf + rhol*Mlf);
    phil = phiP * Mlf / (rhogf*Mgf + rhol*Mlf);
    
    // Explicitly enforce face fluxes at boundaries with prescribed velocity
    forAll(mesh.boundary(), patchi)
    {
        if (Ul.boundaryField()[patchi].fixesValue())
        {
            phil.boundaryFieldRef()[patchi] =
                Ul.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (Ug.boundaryField()[patchi].fixesValue())
        {
            phig.boundaryFieldRef()[patchi] =
                Ug.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }

    Ug = fvc::reconstruct(phig);
    Ul = fvc::reconstruct(phil);

    Ug.correctBoundaryConditions();
    Ul.correctBoundaryConditions();

    if (!incompressible)
    {
        CPerm.correctDensity();
    }
}
