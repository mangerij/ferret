// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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

#include <limits>
#include <iostream>


#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Files/FTreeBuilder.hpp"

int main(int argc, char** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test the insertion of particle using TreeBuilder.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;

    static const int P = 9;

    typedef FRotationCell<FReal,P>               CellClass;
    typedef FP2PParticleContainer<FReal>          ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");

    // -----------------------------------------------------

    FFmaGenericLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    {
        ContainerClass particles;
        particles.reserve(loader.getNumberOfParticles());

        FTic timer;

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> particlePosition;
            FReal physicalValue;
            loader.fillParticle(&particlePosition,&physicalValue);
            particles.push(particlePosition, physicalValue );
        }
        std::cout << "Load the file in " << timer.tacAndElapsed() << "s\n";

        timer.tic();
        {
            const FReal*const partX = particles.getPositions()[0];
            const FReal*const partY = particles.getPositions()[1];
            const FReal*const partZ = particles.getPositions()[2];
            const FReal*const physicalValues = particles.getPhysicalValues();

            OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

            for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                const FPoint<FReal> particlePosition(
                            partX[idxPart], partY[idxPart], partZ[idxPart]
                            );
                tree.insert(particlePosition, physicalValues[idxPart] );
            }
        }
        std::cout << "Create the tree in sequential " << timer.tacAndElapsed() << "s\n";

        timer.tic();
        {
            OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
            // This will modify particles!
            FTreeBuilder<FReal,OctreeClass, LeafClass>::BuildTreeFromArray(&tree, particles);
        }
        std::cout << "Create the tree in parallel " << timer.tacAndElapsed() << "s\n";

    }

    // -----------------------------------------------------

    return 0;
}

