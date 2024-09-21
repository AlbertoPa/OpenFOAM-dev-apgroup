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
#include "fvc.H"
#include "fvm.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudoTurbulentHeatTransferModel, 0);
    defineRunTimeSelectionTable(pseudoTurbulentHeatTransferModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudoTurbulentHeatTransferModel::pseudoTurbulentHeatTransferModel
(
    const dictionary& dict,
    const phaseModel& phase,
    const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
)
:
    dict_(dict),
    phase_(phase),
    alphapt_
    (
         IOobject
        (
            IOobject::groupName("alphapt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_().mesh(),
        dimensionedTensor("alphaptZero", dimKinematicViscosity, tensor::zero)
    ),
    pseudoTurbulentMomentum_(pseudoTurbulentMomentum)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudoTurbulentHeatTransferModel::~pseudoTurbulentHeatTransferModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volTensorField& Foam::pseudoTurbulentHeatTransferModel
::alphapt() const
{
    return alphapt_;
}

Foam::tmp<Foam::volScalarField> Foam::pseudoTurbulentHeatTransferModel
::divqpt(volScalarField& he) const
{
    return
       -fvc::laplacian
        (
            phase_*phase_.rho()*phase_.thermo().Cpv()*alphapt_, 
            phase_.thermo().T()
        );
       //-fvm::laplacianCorrection
       // (
       //     phase_*phase_.rho()*alphapt_,  //thermo.kappa()/thermo.Cpv()
       //     he
       // );
}


// ************************************************************************* //
