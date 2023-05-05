{
	fvScalarMatrix xiEqn
	(
        nf*fvm::ddt(xi) + fvc::div(Ul)
	);

	xiEqn.solve();

	Info<< "Saturation: " << " Min(xi) = " << gMin(xi)
		<< " Max(xi) = " << gMax(xi) << endl;
}
