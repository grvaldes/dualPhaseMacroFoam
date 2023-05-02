#include "fvCFD.H"
#include "cRelPerm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
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

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "xiEqn.H"
        #include "updateSaturationProperties.H"
        #include "pEqn.H"
        
	runTime.write();
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
