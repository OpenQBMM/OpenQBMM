    runTime.write();

    // The zero-order moment of the distribution, whatever the number of
    // its dimensions
    const wordList moment0Names
    (
        mesh.sortedNames<volScalarField>
        (
            wordRe("moment\\.0(_0)*\\.mixing", wordRe::REGEX)
        )
    );

    if (moment0Names.size())
    {
        Info<< moment0Names[0] << " = "
            << mesh.lookupObject<volScalarField>(moment0Names[0])[0] << endl;
    }

