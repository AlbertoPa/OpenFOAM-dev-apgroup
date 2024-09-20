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

Class
    Foam::pseudoTurbulentMomentumModels::pseudoTurbulentMomentumModel

SourceFiles
    pseudoTurbulentMomentumModel.C

\*---------------------------------------------------------------------------*/

#ifndef pseudoTurbulentMomentumModel_H
#define pseudoTurbulentMomentumModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class pseudoTurbulentMomentumModel Declaration
\*---------------------------------------------------------------------------*/

class pseudoTurbulentMomentumModel
{
protected:

    // Protected data

        const dictionary& dict_;

        //- Phase to which the pseudotubulence model is applied
        const phaseModel& phase_;

        //- Name of the continuous phase
        word dispersePhaseName_;

        //- Total disperse phase volume fraction
        volScalarField alphaDisperse_;

        //- Relative velocity magnitude
        volVectorField Ur_;

        //- Slip velocity magnitude
        volScalarField magUr_;

        //- Slip Reynolds number
        volScalarField Rem_;

        //- Pseudoturbulent kinetic energy
        volScalarField kpt_;

        //- Pseudoturbulent stress tensor before transoformation
        volSymmTensorField RptDiagonal_;
       
        //- Pseudoturbulent stress tensor after transformation
        volTensorField Rpt_;

        //- Transformation tensor field
        volTensorField transformationTensor_;

        //- Update the transformation tensor field
        void calcTransformationTensor();

public:

    //- Runtime type information
    TypeName("pseudoTurbulentMomentumModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pseudoTurbulentMomentumModel,
        dictionary,
        (
            const dictionary& dict,
            const phaseModel& phase
        ),
        (dict, phase)
    );


    // Constructors

        //- Construct from components
        pseudoTurbulentMomentumModel
        (
            const dictionary& dict,
            const phaseModel& phase
        );

        //- Disallow default bitwise copy construction
        pseudoTurbulentMomentumModel
        (
            const pseudoTurbulentMomentumModel&
        ) = delete;


    // Selectors

        static autoPtr<pseudoTurbulentMomentumModel> New
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~pseudoTurbulentMomentumModel();


    // Member Functions

        const volScalarField& alphaDisperse() const;

        const volScalarField& magUr() const;

        //- Pseudoturbulent kinetic energy field
        virtual const volScalarField& kpt() const;

        //- Pseudoturbulent stress tensor field before transformation
        virtual const volSymmTensorField& RptDiagonal() const;

        //- Pseudoturbulent stress tensor field after tranformation
        virtual const volTensorField& Rpt() const;

        //- Divergence of the pseudoturbulent stress tensor field
        virtual tmp<volVectorField> divRpt() const;

        //- Transformation matrix
        virtual const volTensorField& transformationTensor() const;

        //- Return the disperse phase model
        //  which for two-phases is the "other" phase
        //  and for more than two phases must be specified
        const phaseModel& dispersePhase() const;

        //- Dispersed phase name
        const word& dispersePhaseName() const;

        //- Slip Reynolds number
        const volScalarField& Rem() const;

        //- Update pseudoturbulent model
        virtual void correct() = 0;

        virtual bool read()
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pseudoTurbulentMomentumModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
