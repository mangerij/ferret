// ===================================================================================
// Copyright ScalFmm 2014 I
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

#ifndef FUSERLEAFCONTAINER_HPP
#define FUSERLEAFCONTAINER_HPP

/**
 * @file This file contains a class that is an implementation of
 * FParticleContainer, designed to be used by the User Kernel API.
 */

#include <Components/FBasicParticleContainer.hpp>


/**
 * @author Piacibello
 *
 * @brief This class define another Particle Container, with dynamic
 * (i.e. no template static) storage in order to store the user
 * particles informations
 */
template<class FReal>
class FUserLeafContainer : public FP2PParticleContainerIndexed<FReal>{

    void * userAttributes;

public:
    FUserLeafContainer(const FUserLeafContainer&) = delete;
    FUserLeafContainer& operator =(const FUserLeafContainer&) = delete;

    FUserLeafContainer() : userAttributes(nullptr){

    }

    void setContainer(void * inputPtr){
        userAttributes = inputPtr;
    }
    void * getContainer() const {
        return userAttributes;
    }

};



#endif // FUSERLEAFCONTAINER_HPP
