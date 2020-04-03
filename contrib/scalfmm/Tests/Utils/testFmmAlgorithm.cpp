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


#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/** This program show an example of use of the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test FMM algorithm by counting the nb of interactions each particle receive.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbParticles);

    typedef double FReal;
    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart       = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000000));
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoader<FReal> loader(NbPart, 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    {
        FPoint<FReal> particlePosition;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition);
            tree.insert(particlePosition);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgo<OctreeClass, CellClass, ContainerClass, LeafClass>(&tree);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



