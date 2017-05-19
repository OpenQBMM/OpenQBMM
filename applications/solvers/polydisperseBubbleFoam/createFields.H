#include "readGravitationalAcceleration.H"
#include "readhRef.H"

    // Create fields for 2-fluid solver
    IOdictionary phaseDict
    (
        IOobject
        (
            "phaseProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //- Total U
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, Zero)
    );

    //- Total flux
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    );

    wordList phases(phaseDict.lookup("phases"));
    word dispersed(phases[0]);
    word continuous(phases[1]);

    // Create dispersed phase fields (uses quadrature)
    pdPhaseModel phase1(mesh,phaseDict,dispersed);
    volScalarField& alpha1 = phase1;
    volVectorField& U1 = phase1.U();
    surfaceScalarField& phi1 = phase1.phi();
    surfaceScalarField& alphaPhi1 = phase1.alphaPhi();
    surfaceScalarField& alphaRhoPhi1 = phase1.alphaRhoPhi();

    const label nNodes = phase1.nNodes();

    // Create continuous phase fields
    phaseModel phase2(mesh,phaseDict,continuous);
    volScalarField& alpha2 = phase2;
    alpha2 == scalar(1) - alpha1;
    volVectorField& U2 = phase2.U();
    surfaceScalarField& phi2 = phase2.phi();
    surfaceScalarField& alphaPhi2 = phase2.alphaPhi();
    surfaceScalarField& alphaRhoPhi2 = phase2.alphaRhoPhi();

    // Set information used in phase pair class
    phasePair::dictTable aspectRatioTable(phaseDict.lookup("aspectRatio"));
    phasePair::scalarTable sigmaTable(phaseDict.lookup("sigma"));
    HashTable<autoPtr<blendingMethod>, word, word::hash> blendingMethods;

    forAllConstIter(dictionary, phaseDict.subDict("blending"), iter)
    {
        blendingMethods.insert
        (
            iter().dict().dictName(),
            blendingMethod::New
            (
                iter().dict(),
                wordList(phaseDict.lookup("phases"))
            )
        );
    }

    // Create phase pairs
    orderedPhasePair pair1In2
    (
        phase1,
        phase2,
        g,
        sigmaTable,
        aspectRatioTable
    );
    orderedPhasePair pair2In1
    (
        phase2,
        phase1,
        g,
        sigmaTable,
        aspectRatioTable
    );
    phasePair pair1And2
    (
        phase1,
        phase2,
        g,
        sigmaTable
    );


    // Create interfacial models
    autoPtr<BlendedInterfacialModel<dragModel> > drag
    (
        new BlendedInterfacialModel<dragModel>
        (
            phaseDict.lookup("drag"),
            (
                blendingMethods.found("drag")
              ? blendingMethods["drag"]
              : blendingMethods["default"]
            ),
            pair1And2,
            pair1In2,
            pair2In1,
            false // Do not zero drag coefficent at fixed-flux BCs
        )
    );

    autoPtr<BlendedInterfacialModel<turbulentDispersionModel> > dispersion
    (
        new BlendedInterfacialModel<turbulentDispersionModel>
        (
            phaseDict.lookup("turbulentDispersion"),
            (
                blendingMethods.found("turbulentDispersion")
              ? blendingMethods["turbulentDispersion"]
              : blendingMethods["default"]
            ),
            pair1And2,
            pair1In2,
            pair2In1
        )
    );

    autoPtr<BlendedInterfacialModel<virtualMassModel> > virtualMass
    (
        new BlendedInterfacialModel<virtualMassModel>
        (
            phaseDict.lookup("virtualMass"),
            (
                blendingMethods.found("virtualMass")
              ? blendingMethods["virtualMass"]
              : blendingMethods["default"]
            ),
            pair1And2,
            pair1In2,
            pair2In1
        )
    );

    autoPtr<BlendedInterfacialModel<liftModel> > lift
    (
        new BlendedInterfacialModel<liftModel>
        (
            phaseDict.lookup("lift"),
            (
                blendingMethods.found("lift")
              ? blendingMethods["lift"]
              : blendingMethods["default"]
            ),
            pair1And2,
            pair1In2,
            pair2In1
        )
    );

    autoPtr<BlendedInterfacialModel<wallLubricationModel> > wallLubrication
    (
        new BlendedInterfacialModel<wallLubricationModel>
        (
            phaseDict.lookup("wallLubrication"),
            (
                blendingMethods.found("wallLubrication")
              ? blendingMethods["wallLubrication"]
              : blendingMethods["default"]
            ),
            pair1And2,
            pair1In2,
            pair2In1
        )
    );

    // Rhos
    volScalarField rho1 = phase1.rho();
    volScalarField rho2 = phase2.rho();

    // Nus
    volScalarField nu1 = phase1.nu();
    volScalarField nu2 = phase2.nu();

    //- Total U
    U == alpha1*U1 + alpha2*U2;

    //- Total flux
    phi == fvc::interpolate(alpha1)*phi1 + fvc::interpolate(alpha2)*phi2;

    //- Total rho
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    );

    // Pressure
    dimensionedScalar pMin
    (
        "pMin",
        dimPressure,
        0.0
    );

    #include "gh.H"

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    p_rgh = p - gh*rho;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );
    mesh.setFluxRequired(p_rgh.name());
    mesh.setFluxRequired(alpha1.name());


#include "createMRF.H"
