/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "MaxwellianThermal.H"
#include "constants.H"
#include "standardNormal.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MaxwellianThermal<CloudType>::MaxwellianThermal
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallInteractionModel<CloudType>(cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MaxwellianThermal<CloudType>::~MaxwellianThermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MaxwellianThermal<CloudType>::correct
(
    typename CloudType::parcelType& p
)
{
    const polyMesh& mesh = this->owner().mesh();

    vector& U = p.U();

    scalar& Ei = p.Ei();

    label typeId = p.typeId();

    const label wppIndex = p.patch(mesh);

    const polyPatch& wpp = mesh.boundaryMesh()[wppIndex];

    label wppLocalFace = wpp.whichFace(p.face());

    const vector nw = p.normal(mesh);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    CloudType& cloud(this->owner());

    randomGenerator& rndGen = cloud.rndGen();
    distributions::standardNormal& stdNormal = cloud.stdNormal();

    while (mag(Ut) < small)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.scalar01()),
            U.y()*(0.8 + 0.2*rndGen.scalar01()),
            U.z()*(0.8 + 0.2*rndGen.scalar01())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    scalar T = cloud.boundaryT().boundaryField()[wppIndex][wppLocalFace];

    scalar mass = cloud.constProps(typeId).mass();

    direction iDof = cloud.constProps(typeId).internalDegreesOfFreedom();

    U =
        sqrt(physicoChemical::k.value()*T/mass)
       *(
            stdNormal.sample()*tw1
          + stdNormal.sample()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.scalar01(), vSmall)))*nw
        );

    U += cloud.boundaryU().boundaryField()[wppIndex][wppLocalFace];

    Ei = cloud.equipartitionInternalEnergy(T, iDof);
}


// ************************************************************************* //
