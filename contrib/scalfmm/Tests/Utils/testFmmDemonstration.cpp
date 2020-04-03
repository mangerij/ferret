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
#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"


#include "../../Src/Components/FAbstractParticleContainer.hpp"
#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

// My cell is actually a basic cell => minimum of data
class MyCell : public FBasicCell {
};

// My fack Container is simply saving the indexes of each
// particles (nothing more!)
template <class FReal>
class MyContainer : public FAbstractParticleContainer<FReal> {
    FVector<FSize> indexes;
public:
    template<typename... Args>
    void push(const FPoint<FReal>& /*inParticlePosition*/, const FSize newParticleIndex, Args ... /*args*/){
        indexes.push(newParticleIndex);
    }

    FSize getSize() const {
        return indexes.getSize();
    }

    const FVector<FSize>& getIndexes() const{
        return indexes;
    }
};

// My leaf process the particles and save only
// those where keepIt is true (during the push method)
template <class FReal>
class MyLeaf : public FAbstractLeaf< FReal, MyContainer<FReal> > {
    MyContainer<FReal> particles;

public:
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const bool keepIt, Args ... args){
        if(keepIt) particles.push(inParticlePosition, args...);
    }
    MyContainer<FReal>* getSrc(){
        return &particles;
    }
    MyContainer<FReal>* getTargets(){
        return &particles;
    }
};

// My kernel actually does nothing except showing how to retreive data from an
// array from the indexes vector giving by the leaf in the P2M
template< class CellClass, class ContainerClass>
class MyKernel : public FAbstractKernels<CellClass,ContainerClass>{
    MortonIndex* indexForEachParticle;
public:
    MyKernel(const FSize inNbParticles): indexForEachParticle(new MortonIndex[inNbParticles]) {
        memset(indexForEachParticle,0,sizeof(MortonIndex)*inNbParticles);
    }

    ~MyKernel(){
        delete[] indexForEachParticle;
    }

    void P2M(CellClass* const cell, const ContainerClass* const particles)  override  {
        for(FSize idxPart = 0 ; idxPart < particles->getSize() ; ++idxPart){
            // save the current morton index for each particles
            indexForEachParticle[ particles->getIndexes()[idxPart] ] = cell->getMortonIndex();
        }
    }

    void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int )  override  {
    }

    void M2L(CellClass* const FRestrict , const CellClass* [], const int [], const int , const int )  override  {
    }

    void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) override  {
    }

    void L2P(const CellClass* const , ContainerClass* const ) override {
    }

    void P2P(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [], const int [], const int ) override {

    }


    void P2POuter(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict ,
                     ContainerClass* const [],
                    const int [], const int ) override {

    }

};


int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Explains how to use ScalFMM (only the code is interesting).",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbParticles);

    typedef double FReal;
    // Custom data structure here
    typedef MyCell            CellClass;
    typedef MyContainer<FReal>       ContainerClass;
    typedef MyLeaf<FReal>            LeafClass;

    // Standard things here
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef MyKernel< CellClass, ContainerClass >         KernelClass;
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    //////////////////////////////////////////////////////////////
    ///////////////////////What we do/////////////////////////////

    std::cout << ">> This executable has to be used to demonstrate the use of scalfmm.\n";

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoader<FReal> loader(NbPart, 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    FPoint<FReal>*const realsParticlesPositions = new FPoint<FReal>[NbPart];
    {
        FPoint<FReal> particlePosition;
        bool keepIt;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            // get a random position
            loader.fillParticle(&particlePosition);
            // let say we remove 1/5 particles
            keepIt = (idxPart%5);
            // insert in the tree
            tree.insert(particlePosition, keepIt, idxPart);
            // save the position
            realsParticlesPositions[idxPart] = particlePosition;
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Which particles are in wich leafs ..." << std::endl;
    counter.tic();

    OctreeClass::Iterator octreeIterator(&tree);
    octreeIterator.gotoBottomLeft();
    do{
        FVector<FSize>::ConstBasicIterator iter(octreeIterator.getCurrentListTargets()->getIndexes());
        const MortonIndex indexAtThisLeaf = octreeIterator.getCurrentGlobalIndex();

        while( iter.hasNotFinished() ){
            std::cout << "Particles with index " << iter.data() <<
                         " has a morton index of " << indexAtThisLeaf << std::endl;
            std::cout << " it real position was " << realsParticlesPositions[iter.data()] << std::endl;
            iter.gotoNext();
        }
    } while(octreeIterator.moveRight());

    counter.tac();
    std::cout << "Done  " << "(@Counting = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(NbPart);
    FmmClass algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    delete [] realsParticlesPositions;

    return 0;
}



