/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "Higbie.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusiveMassTransferModels
{
    defineTypeNameAndDebug(Higbie, 0);
    addToRunTimeSelectionTable
    (
        diffusiveMassTransferModel,
        Higbie,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::Higbie::Higbie
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    diffusiveMassTransferModel(dict, interface),
    interface_
    (
        interface.modelCast
        <
            diffusiveMassTransferModel,
            dispersedPhaseInterface
        >()
    ),
    Le_("Le", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::Higbie::~Higbie()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diffusiveMassTransferModels::Higbie::K() const
{
    const volScalarField Sh
    (
        1.13*sqrt(interface_.Re()*Le_*interface_.Pr())
    );

    return 6.0*interface_.dispersed()*Sh/sqr(interface_.dispersed().d());
}


// ************************************************************************* //