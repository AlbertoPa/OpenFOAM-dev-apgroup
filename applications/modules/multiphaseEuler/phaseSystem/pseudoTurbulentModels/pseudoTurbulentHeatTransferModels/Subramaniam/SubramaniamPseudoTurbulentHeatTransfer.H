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
    Foam::pseudoTurbulentHeatTransferModels::pseudoTurbulenceHeatTransferModel::Subramaniam

Description
    Pseudoturbulent heat transfer model of Zhou et al. (2024).

    References:
    \verbatim
        Zhou, J., Sun, B., & Subramaniam, S. (2024). 
        Modeling pseudo-turbulent heat flux in gas-solid heat transfer. 
        Chemical Engineering Science, 283, 119371. 
    \verbatim

SourceFiles
    Subramaniam.C

\*---------------------------------------------------------------------------*/

#ifndef Subramaniam_H
#define Subramaniam_H

#include "pseudoTurbulentHeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pseudoTurbulentHeatTransferModels
{

/*---------------------------------------------------------------------------*\
                      Class Subramaniam Declaration
\*---------------------------------------------------------------------------*/

class Subramaniam
:
    public pseudoTurbulentHeatTransferModel
{
public:

    //- Runtime type information
    TypeName("Subramaniam");


    // Constructors

        //- Construct from components
        Subramaniam
        (
            const dictionary& dict,
            const phaseModel& phase,
            const pseudoTurbulentMomentumModel& pseudoTurbulentMomentum
        );


    //- Destructor
    virtual ~Subramaniam();


    // Member Functions

        //- Update pseudoturbulent model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pseudoTurbulentHeatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
