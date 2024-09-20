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

Class
    Foam::diffusiveMassTransferModels::Jajuee

Description
    Jajuee correlation for mass transfer from the surface of a sphere to the 
    surrounding fluid.

    Reference:
    \verbatim
        Jajuee, B., Margaritis, A., Karamanev, D., Bergougnou, M.A., (2006).
        Application of surface-renewal-stretch model for interface mass transfer.
        Chem. Eng. Sci. 61, pp. 3917–3929.
    \endverbatim

SourceFiles
    Jajuee.C

\*---------------------------------------------------------------------------*/

#ifndef Jajuee_H
#define Jajuee_H

#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace diffusiveMassTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class Jajuee Declaration
\*---------------------------------------------------------------------------*/

class Jajuee
:
    public diffusiveMassTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;

        //- Lewis number
        const dimensionedScalar Le_;


public:

    //- Runtime type information
    TypeName("Jajuee");


    // Constructors

        //- Construct from a dictionary and an interface
        Jajuee
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~Jajuee();


    // Member Functions

        //- The implicit mass transfer coefficient
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusiveMassTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //