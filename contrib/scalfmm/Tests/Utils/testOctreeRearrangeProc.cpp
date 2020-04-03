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

// ==== CMAKE =====
// @FUSE_MPI
// ================

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

#include "../../Src/Arranger/FOctreeArrangerProc.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

template <class FReal>
struct TestParticle{
    FPoint<FReal> position;
    const FPoint<FReal>& getPosition(){
        return position;
    }
};

template <class FReal, class ParticleClass>
class Converter {
public:
    template <class ContainerClass>
    static ParticleClass GetParticleAndRemove(ContainerClass* container, const FSize idxExtract){
        TestParticle<FReal> part;
        part.position.setPosition(
                    container->getPositions()[0][idxExtract],
                container->getPositions()[1][idxExtract],
                container->getPositions()[2][idxExtract]);
        container->removeParticles(&idxExtract, 1);
        return part;
    }

    template <class OctreeClass>
    static void Insert(OctreeClass* tree, const ParticleClass& part){
        tree->insert(part.position);
    }
};



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "In distributed!\n"
                         "Put the particles into a tree, then change the position of some particles and update the tree.\n"
                         "This method should be used to avoid the tree reconstruction.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;

    typedef FTestCell                   CellClass;
    typedef FBasicParticleContainer<FReal,0,FReal>      ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app(argc, argv);

    const int NbLevels          = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart            = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20000));

    FTic counter;

    srand48 ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal BoxWidth = 1.0;
    const FReal BoxCenter = 0.5;

    OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, FPoint<FReal>(BoxCenter,BoxCenter,BoxCenter));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    counter.tic();


    {
        TestParticle<FReal>* particles = new TestParticle<FReal>[NbPart];
        for(FSize idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particles[idxPart].position.setPosition(
                        (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/FReal(2.0))),
                        (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/FReal(2.0))),
                        (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/FReal(2.0))));
        }

        FVector<TestParticle<FReal>> finalParticles;
        FLeafBalance balancer;
        // FMpiTreeBuilder< FReal,TestParticle<FReal> >::ArrayToTree(app.global(), particles, NbPart,
        // 					     FPoint<FReal>(BoxCenter,BoxCenter,BoxCenter),
        // 					     BoxWidth, tree.getHeight(), &finalParticles,&balancer);
        FMpiTreeBuilder< FReal,TestParticle<FReal> >::DistributeArrayToContainer(app.global(),particles,
                                                                                 NbPart,
                                                                                 FPoint<FReal>(BoxCenter,BoxCenter,BoxCenter),
                                                                                 BoxWidth,NbLevels,
                                                                                 &finalParticles, &balancer);
        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            tree.insert(finalParticles[idx].position);
        }

        delete[] particles;
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass* particles = octreeIterator.getCurrentListTargets();
            for(FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){
                particles->getWPositions()[0][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
                particles->getWPositions()[1][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
                particles->getWPositions()[2][idxPart] = (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2));
            }
        } while(octreeIterator.moveRight());
    }

    counter.tac();
    std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Arrange ..." << std::endl;
    counter.tic();

    FOctreeArrangerProc<FReal, OctreeClass, ContainerClass, TestParticle<FReal>, Converter<FReal, TestParticle<FReal>> > arrange(&tree);
    arrange.rearrange(app.global());

    counter.tac();
    std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Test ..." << std::endl;
    counter.tic();

    { // Check that each particle has been put into the right leaf
        long counterPart = 0;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        MortonIndex interval[2] = {0,0};
        interval[0] = octreeIterator.getCurrentGlobalIndex();

        do{
            const MortonIndex leafIndex = octreeIterator.getCurrentGlobalIndex();

            ContainerClass* particles = octreeIterator.getCurrentListTargets();
            for(FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){
                const FPoint<FReal> particlePosition( particles->getWPositions()[0][idxPart],
                        particles->getWPositions()[1][idxPart],
                        particles->getWPositions()[2][idxPart]);

                const MortonIndex particleIndex = tree.getMortonFromPosition( particlePosition );
                if( leafIndex != particleIndex){
                    std::cout << "Index problem, should be " << leafIndex <<
                                 " particleIndex "<< particleIndex << std::endl;
                }

            }

            counterPart += octreeIterator.getCurrentListTargets()->getNbParticles();
            if(octreeIterator.getCurrentListTargets()->getNbParticles() == 0){
                std::cout << "Problem, leaf is empty at index " << leafIndex << std::endl;
            }
        } while(octreeIterator.moveRight());
        interval[1] = octreeIterator.getCurrentGlobalIndex();

        counterPart = app.global().reduceSum(counterPart);
        if(app.global().processId() == 0 && counterPart != NbPart  * app.global().processCount() ){
            std::cout <<"Wrong particles number, should be " << (NbPart  * app.global().processCount()) << " but is " << counterPart << std::endl;
        }

        MortonIndex*const allintervals = new MortonIndex[ 2 * app.global().processCount() ];
        MPI_Allgather( interval, sizeof(MortonIndex) * 2, MPI_BYTE, allintervals, sizeof(MortonIndex) * 2, MPI_BYTE, app.global().getComm());
        if(app.global().processId() == 0){
            for(int idxProc = 1 ; idxProc < app.global().processCount() ; ++idxProc){
                if( allintervals[idxProc*2-1] > allintervals[idxProc*2] ){
                    std::cout << "Interval problem for [" << idxProc-1 << "].max = " << allintervals[idxProc*2-1] <<
                                 "  [" << idxProc << "].min = "<< allintervals[idxProc*2] << std::endl;
                }
            }
        }
        delete[] allintervals;
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



