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

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    {
        FHelpDescribeAndExit(argc, argv, "Show the memory usage (if mem stats is turned on at the compilation)",
                             FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                             FParameterDefinitions::NbParticles, FParameterDefinitions::SHDevelopment);

        typedef double FReal;
        typedef FSphericalCell<FReal>                 CellClass;

        typedef FP2PParticleContainer<FReal>      ContainerClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
        //typedef FTestKernels< CellClass, ContainerClass >         KernelClass;
        typedef FSphericalKernel<FReal, CellClass, ContainerClass >          KernelClass;

        typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
        ///////////////////////What we do/////////////////////////////
        std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
        //////////////////////////////////////////////////////////////

        const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
        const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
        const FSize NbPart       = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000000));
        const int DevP         = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 5);
        const FPoint<FReal> centerOfBox = FPoint<FReal>(0.5,0.5,0.5);
        FTic counter;

        srand48 ( 1 ); // volontary set seed to constant

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        CellClass::Init(DevP);
        const FReal boxWidth = 1.0;
        OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, centerOfBox);

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
        std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
        counter.tic();

        {
            FPoint<FReal> particlePosition;
            FReal physicalValue = 0.10;
            for(FSize idxPart = 0 ; idxPart < NbPart ; ++idxPart){
                particlePosition.setPosition(FReal(drand48()),FReal(drand48()),FReal(drand48()));
                tree.insert(particlePosition, physicalValue);
            }
        }

        counter.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Working on particles ..." << std::endl;
        counter.tic();

        // FTestKernels FBasicKernels
        //KernelClass kernels;
        KernelClass kernels(DevP, NbLevels,boxWidth, centerOfBox);
        FmmClass algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////
    }

    std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
    std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
    std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";

    return 0;
}



