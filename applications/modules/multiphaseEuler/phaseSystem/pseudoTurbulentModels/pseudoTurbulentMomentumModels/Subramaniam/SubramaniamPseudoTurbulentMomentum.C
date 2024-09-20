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

#include "SubramaniamPseudoTurbulentMomentum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pseudoTurbulentMomentumModels
{
    defineTypeNameAndDebug(Subramaniam, 0);

    addToRunTimeSelectionTable
    (
        pseudoTurbulentMomentumModel,
        Subramaniam,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudoTurbulentMomentumModels::Subramaniam
::Subramaniam
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    pseudoTurbulentMomentumModel(dict, phase)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudoTurbulentMomentumModels::Subramaniam
::~Subramaniam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::pseudoTurbulentMomentumModels::Subramaniam
::correct()
{
    calcTransformationTensor();

    const volScalarField alphac(phase_);
    const phaseModel& disperse = dispersePhase();
    const volScalarField& alphad(alphaDisperse());

    volScalarField Ef(0.5*sqr(magUr_)); 

    Rem_ = (alphac*disperse.d()*magUr_*phase_.rho())/phase_.fluidThermo().mu();
    
    scalar minRem = min(Rem_).value();
    scalar maxRem = max(Rem_).value();

    if (minRem < 1 || maxRem > 300.0)
    {
        WarningInFunction
            << "Slip Reynolds number larger out of validity range " 
            << "(0.01 < Rem < 300). "
            << "    Min. slip Reynolds number = " << minRem << endl
            << "    Max. slip Reynolds number = " << maxRem << endl;
    }

    // Calculate pseudoturbulent kinetic energy
    kpt_ = 
        Ef
       *pos(alphac - 0.5) // Remove where alphac is too small
       *(
            2.0*alphad + 2.5*alphad*(pow3(alphac))*exp(-alphad*sqrt(Rem_))
        );

    volScalarField bParallel
    (
        kpt_*
        (
            (0.523/(1.0 + 0.305*exp(-0.114*Rem_)))
           *exp((-3.511*(alphad))/(1.0 + 1.801*exp(-0.005*Rem_)))
        )
    );

    forAll (Rpt_, celli)
    {
        RptDiagonal_[celli].component(symmTensor::XX) = 
            2.0*bParallel[celli] + (2.0/3.0)*kpt_[celli];

        RptDiagonal_[celli].component(symmTensor::YY) = 
            -bParallel[celli] + (2.0/3.0)*kpt_[celli];

        RptDiagonal_[celli].component(symmTensor::ZZ) = 
            -bParallel[celli] + (2.0/3.0)*kpt_[celli];

        Rpt_[celli] = 
            transform(transformationTensor_[celli], RptDiagonal_[celli]);
    }

    Rpt_.correctBoundaryConditions();
}


// ************************************************************************* //
