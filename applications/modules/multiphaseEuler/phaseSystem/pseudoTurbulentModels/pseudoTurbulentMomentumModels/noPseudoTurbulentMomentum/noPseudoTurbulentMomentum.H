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
    Foam::pseudoTurbulentMomentumModels::pseudoTurbulentMomentumModel::nopseudoTurbulentMomentumModel

Description

SourceFiles
    nopseudoTurbulentMomentumModel.C

\*---------------------------------------------------------------------------*/

#ifndef noPseudoTurbulentMomentumModel_H
#define noPseudoTurbulentMomentumModel_H

#include "pseudoTurbulentMomentumModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pseudoTurbulentMomentumModels
{

/*---------------------------------------------------------------------------*\
                Class noPseudoTurbulentMomentumModel Declaration
\*---------------------------------------------------------------------------*/

class noPseudoTurbulentMomentumModel
:
    public pseudoTurbulentMomentumModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noPseudoTurbulentMomentumModel
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~noPseudoTurbulentMomentumModel();


    // Member Functions

        //- Update pseudoturbulent model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pseudoTurbulentMomentumModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
