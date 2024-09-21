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
    Foam::pseudoTurbulenceHeatTransferModels::pseudoTurbulenceHeatTransferModel

SourceFiles
    pseudoTurbulenceHeatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef pseudoTurbulentHeatTransferModel_H
#define pseudoTurbulentHeatTransferModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "pseudoTurbulentMomentumModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class pseudoTurbulentHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class pseudoTurbulentHeatTransferModel
{
protected:

    // Protected data

        const dictionary& dict_;

        //- Phase to which the pseudotubulence model is applied
        const phaseModel& phase_;

        //- Pseudoturbulent thermal diffusivity energy
        volTensorField alphapt_;

        //- Pseudoturbulent momentum model
        const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum_;

public:

    //- Runtime type information
    TypeName("pseudoTurbulentHeatTransferModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pseudoTurbulentHeatTransferModel,
        dictionary,
        (
            const dictionary& dict,
            const phaseModel& phase,
            const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
        ),
        (dict, phase, pseudoTurbulentMomentum)
    );


    // Constructors

        //- Construct from components
        pseudoTurbulentHeatTransferModel
        (
            const dictionary& dict,
            const phaseModel& phase,
            const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
        );

        //- Disallow default bitwise copy construction
        pseudoTurbulentHeatTransferModel
        (
            const pseudoTurbulentHeatTransferModel&
        ) = delete;


    // Selectors

        static autoPtr<pseudoTurbulentHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseModel& phase,
            const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
        );


    //- Destructor
    virtual ~pseudoTurbulentHeatTransferModel();


    // Member Functions

        //- Pseudoturbulent kinetic energy field
        virtual const volTensorField& alphapt() const;

        //- Pseudoturbulent heat flux
        virtual tmp<volScalarField> divqpt(volScalarField& he) const;

        //- Update pseudoturbulence model
        virtual void correct() = 0;

        virtual bool read()
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pseudoTurbulentHeatTransferModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
