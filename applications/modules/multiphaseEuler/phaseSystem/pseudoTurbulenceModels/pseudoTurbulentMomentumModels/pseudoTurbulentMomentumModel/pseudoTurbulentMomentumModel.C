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

#include "pseudoTurbulentMomentumModel.H"
#include "fvc.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudoTurbulentMomentumModel, 0);
    defineRunTimeSelectionTable(pseudoTurbulentMomentumModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudoTurbulentMomentumModel::pseudoTurbulentMomentumModel
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    dict_(dict),
    phase_(phase),
    dispersePhaseName_
    (
        dict_.lookupOrDefault("dispersePhase", word::null)
    ),
    alphaDisperse_
    (
        1.0 - phase_
    ),
    Ur_
    (
        IOobject
        (
            IOobject::groupName("Urpt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_().mesh(),
        dimensionedVector("UrZero", dimVelocity, vector::zero)
    ),
    magUr_ 
    (
        IOobject
        (
            IOobject::groupName("magUrpt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_().mesh(),
        dimensionedScalar("UrZero", dimVelocity, scalar(0))
    ),
    Rem_  // Should recover from the momentum exchange if possible
    (
        IOobject
        (
            IOobject::groupName("Rempt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_().mesh(),
        dimensionedScalar("RemZero", dimless, scalar(0))
    ),
    kpt_
    (
        IOobject
        (
            IOobject::groupName("kpt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_().mesh(),
        dimensionedScalar("kptZero", sqr(dimVelocity), scalar(0))
    ),
    RptDiagonal_
    (
        IOobject
        (
            IOobject::groupName("Rpt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_().mesh(),
        dimensionedSymmTensor("RptZero", sqr(dimVelocity), symmTensor::zero)
    ),
    Rpt_
    (
        IOobject
        (
            IOobject::groupName("Rpt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_().mesh(),
        dimensionedTensor("RptZero", sqr(dimVelocity), tensor::zero)
    ),
    transformationTensor_
    (
        IOobject
        (
            IOobject::groupName("TransformationTensorPT", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_().mesh(),
        dimensionedTensor("TransformationTensorPT", dimless, tensor::zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudoTurbulentMomentumModel::~pseudoTurbulentMomentumModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseModel& 
Foam::pseudoTurbulentMomentumModel::dispersePhase() const
{
    const phaseSystem& fluid = phase_.fluid();

    if (dispersePhaseName_ == word::null)
    {
        if (fluid.movingPhases().size() != 2)
        {
            FatalIOErrorInFunction(dict_)
                << "Disperse phase name must be specified "
                << "when there are more than two moving phases."
                << exit(FatalIOError);
        }

        forAll(fluid.movingPhases(), movingPhasei)
        {
            const phaseModel& otherPhase = fluid.movingPhases()[movingPhasei];

            if (&otherPhase != &phase_)
            {
                return otherPhase;
            }
        }
    }

    return fluid.phases()[dispersePhaseName_];
}

void Foam::pseudoTurbulentMomentumModel::calcTransformationTensor()
{
    // Retrieve the disperse phase
    const phaseModel& disperse = dispersePhase();

    // Correct fields used in the pseudoturbulence model
    Ur_ = disperse.U() - phase_.U();
    magUr_ = mag(Ur_);
    alphaDisperse_ = 1.0 - phase_;

    // Define othonormal base aligned with the axes
    vector e1(1, 0, 0);
    vector e2(0, 1, 0);
    vector e3(0, 0, 1);  
    
    volVectorField e1Parallel
    (
        IOobject
        (
            IOobject::groupName("kpt", phase_.name()),
            phase_().time().name(),
            phase_().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_().mesh(),
        dimensionedVector("eZero", dimless, vector::zero)
    );

    volVectorField e2Parallel(e1Parallel);
    volVectorField e3Parallel(e2Parallel);

    // Looping to manage the case of small slip velocity, when necessary
    forAll(e1Parallel, celli)
    {
        if (magUr_[celli] > 0.001)
        {
            e1Parallel[celli] = Ur_[celli]/magUr_[celli];

            vector u2 = e2 - (e1Parallel[celli] & e2)
                /magSqr(e1Parallel[celli])*e1Parallel[celli];

            e2Parallel[celli] = u2/mag(u2);

            vector u3 = e3 - (e1Parallel[celli] & e3)
                /magSqr(e1Parallel[celli])*e1Parallel[celli] 
              - (u2 & e3)/magSqr(u2)*u2;
           
            e3Parallel[celli] = u3/mag(u3);

            transformationTensor_[celli] = Foam::tensor
            (
                e1 & e1Parallel[celli], 
                e1 & e2Parallel[celli], 
                e1 & e3Parallel[celli], 
                e2 & e1Parallel[celli], 
                e2 & e2Parallel[celli], 
                e2 & e3Parallel[celli], 
                e3 & e1Parallel[celli],
                e3 & e2Parallel[celli],
                e3 & e3Parallel[celli]
            );
        }
    }

    e1Parallel.correctBoundaryConditions();
    e2Parallel.correctBoundaryConditions();
    e3Parallel.correctBoundaryConditions();
    transformationTensor_.correctBoundaryConditions();
}

const Foam::volScalarField& Foam::pseudoTurbulentMomentumModel
::alphaDisperse() const
{
    return alphaDisperse_;
}

const Foam::volScalarField& Foam::pseudoTurbulentMomentumModel
::magUr() const
{
    return magUr_;
}

const Foam::volScalarField& Foam::pseudoTurbulentMomentumModel
::Rem() const
{
    return Rem_;
}

const Foam::volScalarField& Foam::pseudoTurbulentMomentumModel
::kpt() const
{
    return kpt_;
}

const Foam::volSymmTensorField& Foam::pseudoTurbulentMomentumModel
::RptDiagonal() const
{
    return RptDiagonal_;
}

const Foam::volTensorField& Foam::pseudoTurbulentMomentumModel
::Rpt() const
{
    return Rpt_;
}


const Foam::volTensorField& Foam::pseudoTurbulentMomentumModel
::transformationTensor() const
{
    return transformationTensor_;
}

Foam::tmp<Foam::volVectorField> Foam::pseudoTurbulentMomentumModel
::divRpt() const
{
    return fvc::div(phase_*phase_.rho()*Rpt_);
}



const Foam::word& Foam::pseudoTurbulentMomentumModel::dispersePhaseName() const
{
    return dispersePhaseName_;
}


// ************************************************************************* //
