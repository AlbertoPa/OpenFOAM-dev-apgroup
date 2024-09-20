/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Copyright (C) 2023-2024 Alberto Passalacqua (apcfd@outlook.com)
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pseudoTurbulentHeatTransferModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pseudoTurbulentHeatTransferModel>
Foam::pseudoTurbulentHeatTransferModel::New
(
    const dictionary& dict,
    const phaseModel& phase,
    const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
)
{
    word pseudoTurbulentHeatTransferModelType
    (
        dict.lookup("pseudoTurbulentHeatTransferModel")
    );

    Info<< "Selecting pseudoTurbulentHeatTransferModel "
        << pseudoTurbulentHeatTransferModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pseudoTurbulentHeatTransferModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown pseudoTurbulentHeatTransferModel type "
            << pseudoTurbulentHeatTransferModelType << endl << endl
            << "Valid pseudoTurbulentHeatTransferModel types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pseudoTurbulentHeatTransferModel>
    (
        cstrIter()(dict, phase, pseudoTurbulentMomentum)
    );
}


// ************************************************************************* //
