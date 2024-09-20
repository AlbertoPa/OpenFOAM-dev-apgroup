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

#include "Jajuee.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusiveMassTransferModels
{
    defineTypeNameAndDebug(Jajuee, 0);
    addToRunTimeSelectionTable
    (
        diffusiveMassTransferModel,
        Jajuee,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::Jajuee::Jajuee
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

Foam::diffusiveMassTransferModels::Jajuee::~Jajuee()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diffusiveMassTransferModels::Jajuee::K() const
{
    const phaseCompressible::momentumTransportModel& turbulence =
        interface_.dispersed().mesh().lookupType<phaseCompressible::momentumTransportModel>
        (
            interface_.continuous().name()
        );

    const volScalarField kL
    (
        2.0*pow(turbulence.epsilon(), 0.25)
       *pow(interface_.continuous().rho()/interface_.continuous().fluidThermo().mu(), 0.75)
       *sqrt(Le_*interface_.Pr()/constant::mathematical::pi)
    );

    return 6.0*interface_.dispersed()*kL/interface_.dispersed().d();
}


// ************************************************************************* //