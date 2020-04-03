// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
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
#include "../../Src/Components/FTypedLeaf.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"
#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Arranger/FBasicParticleContainerIndexedMover.hpp"
#include "../../Src/Arranger/FParticleTypedIndexedMover.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Arranger/FOctreeArranger.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


// Simply create particles and try the kernels
int main(int argc, char ** argv){

    FHelpDescribeAndExit(argc, argv,
                         "Put the particles into a tree, then change the position of some particles and update the tree.\n"
                         "This method should be used to avoid the tree reconstruction.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;

    typedef FBasicCell                      CellClass;
    typedef FP2PParticleContainerIndexed<FReal>      ContainerClass;

    typedef FTypedLeaf< FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    typedef FParticleTypedIndexedMover<FReal, OctreeClass, ContainerClass> MoverClass;
    typedef FOctreeArranger<FReal, OctreeClass, ContainerClass, MoverClass> ArrangerClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart_Source     = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(200000));
    const int NbPart_Target     = 10000;

    FTic counter;

    srand48 ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal BoxWidth = 1.0;
    const FReal BoxCenter = 0.5;

    OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, FPoint<FReal>(BoxCenter,BoxCenter,BoxCenter));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart_Source << " particles sources and "<< NbPart_Target <<" particles target" << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();


    {

        FPoint<FReal> particleToFill;
        for(FSize idxPart = 0 ; idxPart < NbPart_Source; ++idxPart){
            particleToFill.setPosition(
                                       (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)),
                                       (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)),
                                       (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)));
            tree.insert(particleToFill,FParticleTypeSource,idxPart);
        }

        for(FSize idxPart = 0 ; idxPart < NbPart_Target; ++idxPart){
            particleToFill.setPosition(
                                       (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)),
                                       (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)),
                                       (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)));
            tree.insert(particleToFill,FParticleTypeTarget,idxPart);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;



    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    ArrangerClass arrange(&tree);

    // In order to test multiple displacement of particles
    //    for(int ite=0 ; ite<10 ; ++ite){
    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    { // Create new position for each particles
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass* particles_sources = octreeIterator.getCurrentListTargets();
            ContainerClass* particles_targets = octreeIterator.getCurrentListSrc();
            for(FSize idxPart = 0; idxPart < particles_sources->getNbParticles() ; ++idxPart){
                particles_sources->getWPositions()[0][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
                particles_sources->getWPositions()[1][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
                particles_sources->getWPositions()[2][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
            }
            for(FSize idxPart = 0; idxPart < particles_targets->getNbParticles() ; ++idxPart){
                particles_targets->getWPositions()[0][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
                particles_targets->getWPositions()[1][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
                particles_targets->getWPositions()[2][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
            }
        } while(octreeIterator.moveRight());
    }

    counter.tac();
    std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Arrange ..." << std::endl;
    counter.tic();

    //FOctreeArranger<FReal,OctreeClass, ContainerClass, TestParticle, Converter<TestParticle> > arrange(&tree);
    arrange.rearrange();


    counter.tac();
    std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //    }
    std::cout << "Test ..." << std::endl;
    counter.tic();

    { // Check that each particle has been put into the right leaf
        long counterPart_source = 0;
        long counterPart_target = 0;
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            const MortonIndex leafIndex = octreeIterator.getCurrentGlobalIndex();

            ContainerClass* particles_targets = octreeIterator.getCurrentLeaf()->getTargets();
            ContainerClass* particles_sources = octreeIterator.getCurrentLeaf()->getSrc();;
            //Check for current sources
            for(FSize idxPart = 0; idxPart < particles_sources->getNbParticles() ; ++idxPart){
                const FPoint<FReal> particlePosition( particles_sources->getWPositions()[0][idxPart],
                                               particles_sources->getWPositions()[1][idxPart],
                                               particles_sources->getWPositions()[2][idxPart]);
                const MortonIndex particleIndex = tree.getMortonFromPosition( particlePosition );
                if( leafIndex != particleIndex){
                    std::cout << "Index problem, should be " << leafIndex <<
                        " particleIndex "<< particleIndex << std::endl;
                }

            }
            //Check for current targets
            for(FSize idxPart = 0; idxPart < particles_targets->getNbParticles() ; ++idxPart){
                const FPoint<FReal> particlePosition( particles_targets->getWPositions()[0][idxPart],
                                               particles_targets->getWPositions()[1][idxPart],
                                               particles_targets->getWPositions()[2][idxPart]);

                const MortonIndex particleIndex = tree.getMortonFromPosition( particlePosition );
                if( leafIndex != particleIndex){
                    std::cout << "Index problem, should be " << leafIndex <<
                        " particleIndex "<< particleIndex << std::endl;
                }

            }

            //Sum the parts
            counterPart_target += octreeIterator.getCurrentListTargets()->getNbParticles();
            counterPart_source += octreeIterator.getCurrentListSrc()->getNbParticles();
            //Check
            if(octreeIterator.getCurrentListTargets()->getNbParticles() == 0
               && octreeIterator.getCurrentListSrc()->getNbParticles() == 0){
                std::cout << "Problem, leaf is empty at index " << leafIndex << std::endl;
            }
        } while(octreeIterator.moveRight());

        if( counterPart_source != NbPart_Source || counterPart_target != NbPart_Target ){
            std::cout <<"Wrong particles number, src : should be " << NbPart_Source << " but is " << counterPart_source << "\n tgt : should be " << NbPart_Target << " but is " << counterPart_target << std::endl;
        }
    }

    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator(&tree);
        OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = NbLevels - 1;
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            // for each cells
            do{
                int countChild = 0;
                CellClass** const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                    if( child[idxChild] ){
                        countChild += 1;
                    }
                }

                if(countChild == 0){
                    std::cout << "Problem at level " << idxLevel << " cell has no child " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Test = " << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}
