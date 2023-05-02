{
	fvScalarMatrix xiEqn
	(
        	nf*fvm::ddt(xi) + fvc::div(phil)
	);

	xiEqn.solve();

	Info<< "Saturation: " << " Min(xi) = " << gMin(xi)
		<< " Max(xi) = " << gMax(xi) << endl;
}
