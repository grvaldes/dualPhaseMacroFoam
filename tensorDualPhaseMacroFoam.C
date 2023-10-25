#include "fvCFD.H"
#include "relPerm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "incompressible",
        "Solves the equation assuming an incompressible gas."
    );
    argList::addBoolOption
    (
        "i"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool incompressible = args.optionFound("i") ||
    args.optionFound("incompressible");
    
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createSaturationProperties.H"
    #include "readTimeControls.H"
    
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime.setDeltaT(runTime.deltaTValue());

        runTime++;

        Info<< nl << "Time = " << runTime.timeName() << endl;

        #include "xiEqn.H"
        #include "updateSaturationProperties.H"
        #include "pEqn.H"
        
	    runTime.write();
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
