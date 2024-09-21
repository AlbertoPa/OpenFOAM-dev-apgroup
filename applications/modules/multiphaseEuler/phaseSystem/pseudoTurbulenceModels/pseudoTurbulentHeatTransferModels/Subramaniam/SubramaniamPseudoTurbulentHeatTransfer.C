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

#include "SubramaniamPseudoTurbulentHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pseudoTurbulentHeatTransferModels
{
    defineTypeNameAndDebug(Subramaniam, 0);

    addToRunTimeSelectionTable
    (
        pseudoTurbulentHeatTransferModel,
        Subramaniam,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudoTurbulentHeatTransferModels::Subramaniam
::Subramaniam
(
    const dictionary& dict,
    const phaseModel& phase,
    const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
)
:
    pseudoTurbulentHeatTransferModel(dict, phase, pseudoTurbulentMomentum)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudoTurbulentHeatTransferModels::Subramaniam
::~Subramaniam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::pseudoTurbulentHeatTransferModels::Subramaniam
::correct()
{
    const volScalarField& alphac(phase_);
    const phaseModel& disperse(pseudoTurbulentMomentum_.dispersePhase());
    const volScalarField& alphad(pseudoTurbulentMomentum_.alphaDisperse());

    scalar minRem = min(pseudoTurbulentMomentum_.Rem()).value();
    scalar maxRem = max(pseudoTurbulentMomentum_.Rem()).value();

    if (minRem < 0.01 || maxRem > 100.0)
    {
        WarningInFunction
            << "Slip Reynolds number larger out of validity range " 
            << "(0.01 < Rem < 300). " << endl
            << "    Min. slip Reynolds number = " << minRem << endl
            << "    Max. slip Reynolds number = " << maxRem << endl;
    }

    const volScalarField& Rem(pseudoTurbulentMomentum_.Rem());

    // Using alphacBound and alphadBound to ensure robustness when dividing for
    // the phase fractions when they are small. Limited versions are used only
    // where needed.
    volScalarField alphacBound(max(alphac, phase_.residualAlpha()));
    volScalarField alphadBound(max(alphad, disperse.residualAlpha()));

    volScalarField Pr
    (
        "Pr", 
        phase_.fluidThermo().mu()*phase_.thermo().Cpv()
       /phase_.fluidThermo().kappa()
    );

    volScalarField alphacSqr("alphacSqr", alphac);

    volScalarField Nu
    (
        "Nu",
        (-0.46 + 1.77*alphac + 0.69*alphacSqr)/pow3(alphacBound)
      + (1.37 - 2.4*alphac + 1.2*alphacSqr)*pow(Rem, 0.7)*pow(Pr, 1.0/3.0)
    );

    volScalarField alphaptParallel
    (
        "alphaptParallel",
        pos(alphac - 0.5)  // Remove when alphac is less than 50%
       *phase_.fluidThermo().kappa()/(phase_.rho()*phase_.thermo().Cpv())
       *(
            2.0*Rem*(Rem + 1.4)*sqr(Pr)
           *(
                alphac*(-5.11*alphad + 10.1*sqr(alphad) - 10.85*pow3(alphad)) 
              + 1.0 - exp(-10.96*alphad)
            )
           *exp(-0.002089*Rem)
        )
       /(
            3.0*M_PI*sqr(alphacBound)*Nu
           *(
                1.17*alphadBound 
              - 0.2021*sqrt(alphadBound) 
              + 0.08568*pow(alphadBound, 1.0/4.0)
            )
           *(
                1.0 
              - 1.6*alphacBound*alphadBound 
              - 3.0*alphadBound*pow4(alphacBound)*exp(-alphadBound*pow(Rem, 0.4))
            )
        )
    );

    const volSymmTensorField& RptDiagonal(pseudoTurbulentMomentum_.RptDiagonal());

    const volTensorField& tranformationTensor
    (
        pseudoTurbulentMomentum_.transformationTensor()
    );

    forAll(alphapt_, celli)
    {
        if (mag(RptDiagonal[celli].component(tensor::XX)) > small)
        {
            alphapt_[celli].component(tensor::XX) = alphaptParallel[celli];

            alphapt_[celli].component(tensor::YY) = 
                alphaptParallel[celli]*RptDiagonal[celli].component(tensor::YY)
               /RptDiagonal[celli].component(tensor::XX);
                    
            alphapt_[celli].component(tensor::ZZ) = 
                alphapt_[celli].component(tensor::YY);
                
            alphapt_[celli] = 
                transform(tranformationTensor[celli], alphapt_[celli]);
        }
    }
}


// ************************************************************************* //
