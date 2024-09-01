/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 Alberto Passalacqua
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

#include "Panicker.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(Panicker, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        Panicker,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Panicker::Panicker
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedTurbulentDispersionModel(dict, interface),
    Cdis_("Cdis", dimless, dict, 4.544)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Panicker::~Panicker()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::Panicker::D() const
{
    // Finding total volume fraction of dispersed phase(s)
    volScalarField alpha1(1.0 - interface_.continuous());

    const dragModels::dispersedDragModel& drag =
        interface_.mesh().lookupObject<dragModels::dispersedDragModel>
        (
            IOobject::groupName
            (
                dragModel::typeName,
                interface_.name()
            )
        );

    scalar b = 0.5;
    scalar a = 1.0 + b - (1.0/3.0);

    return
        drag.Ki()
       *Cdis_
       *interface_.continuous().fluidThermo().nu()
       *interface_.Re()
       *pos0(alpha1 - 0.001)
       *alpha1*(1.0 - a*alpha1 + b*sqr(alpha1));
}


// ************************************************************************* //
