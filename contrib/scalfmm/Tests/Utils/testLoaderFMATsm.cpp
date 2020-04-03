// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
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

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Components/FTypedLeaf.hpp"

#include "../../Src/Files/FFmaTsmLoader.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"

#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


int main(int argc, char ** argv ){
    FHelpDescribeAndExit(argc, argv, "Load a file and put the particles in a tree with the TSM (target source model)",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight);

    typedef double FReal;
    typedef FBasicParticleContainer<FReal,1,FReal>     ContainerClass;
    typedef FTypedLeaf< FReal, ContainerClass >                     LeafClass;
    typedef FOctree< FReal, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Tsm loader\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderCreate.exe to create this file
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.tsm.fma");
    std::cout << "Opening : " << filename << "\n";

    // open basic particles loader
    FFmaTsmLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }
    {
        // otree
        OctreeClass tree(FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5), FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3),
                         loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        FPoint<FReal> particlePosition;
        FReal physicalValue = 0.0;
        FParticleType particleType;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition,&physicalValue, &particleType);
            tree.insert(particlePosition, particleType, physicalValue);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;

        // -----------------------------------------------------
        std::cout << "Deleting particles ..." << std::endl;
        counter.tic();
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------

    return 0;
}



