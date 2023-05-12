{
	fvScalarMatrix xiEqn
	(
        	nf*fvm::ddt(xi) + fvc::div(phil)
	);

	Info<< "phil before" << endl;
    Info<< phil << endl;

	Info<< "xi before" << endl;
    Info<< xi << endl;

	xiEqn.solve();

    Info<< "phil after" << endl;
    Info<< phil << endl;

	Info<< "xi after" << endl;
    Info<< xi << endl;

	Info<< "Ml before" << endl;
    Info<< Ml << endl;
	Info<< "Mg before" << endl;
    Info<< Mg << endl;

	Info<< "Saturation: " << " Min(xi) = " << gMin(xi)
		<< " Max(xi) = " << gMax(xi) << endl;
}
