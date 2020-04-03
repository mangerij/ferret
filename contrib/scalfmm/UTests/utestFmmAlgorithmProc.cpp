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

#include "FUTester.hpp"

#include "Utils/FMpi.hpp"
#include "Utils/FTic.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FGlobal.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Utils/FPoint.hpp"

#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"
#include "Components/FTestParticleContainer.hpp"

#include "Core/FFmmAlgorithmThreadProc.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

#include "Files/FMpiFmaGenericLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "Components/FBasicKernels.hpp"

#include "BalanceTree/FLeafBalance.hpp"

#include <iostream>
#include <cstdio>
#include <cstdlib>


/** this class test the bool array container */
class TestFmmAlgoProc : public FUTesterMpi<TestFmmAlgoProc> {

    // Check if tree is built correctly
    template<class OctreeClass>
    void ValidateTree(OctreeClass& realTree, OctreeClass& treeValide){

        typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
        octreeIteratorValide.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIterator(&realTree);
        octreeIterator.gotoBottomLeft();

        while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            if(octreeIteratorValide.moveRight() == false){
                uassert(false);
                return;
            }
        }

        while(true){
            if(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                uassert(false);
                return;
            }

            if(octreeIteratorValide.getCurrentListSrc()->getNbParticles() != octreeIterator.getCurrentListSrc()->getNbParticles()){
                uassert(false);
                return;
            }

            if(octreeIterator.moveRight() == false){
                break;
            }

            uassert(octreeIteratorValide.moveRight() != false);
        }
    }



    /** This function tests the octree to be sure that the fmm algorithm
     * has worked completly.
     */
    template<class OctreeClass, class ContainerClass, class FmmClassProc>
    void ValidateFMMAlgoProc(OctreeClass* const badTree,
                             OctreeClass* const valideTree,
                             FmmClassProc* const fmm){
        const int OctreeHeight = badTree->getHeight();
        {
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator octreeIteratorValide(valideTree);
            octreeIteratorValide.gotoBottomLeft();

            for(int level = OctreeHeight - 1 ; level > 0 && fmm->hasWorkAtLevel(level) ; --level){

                while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()) {
                    octreeIteratorValide.moveRight();
                }

                while(octreeIteratorValide.getCurrentGlobalIndex() != fmm->getWorkingInterval(level).leftIndex){
                    octreeIteratorValide.moveRight();
                    octreeIterator.moveRight();
                }

                FSize countCheck = 0;
                do{
                    if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                        uassert(false);
                    }
                    else{
                        uassert(octreeIterator.getCurrentCell()->getDataUp() == octreeIteratorValide.getCurrentCell()->getDataUp());
                        uassert(octreeIterator.getCurrentCell()->getDataDown() == octreeIteratorValide.getCurrentCell()->getDataDown());
                    }
                    ++countCheck;
                } while(octreeIteratorValide.moveRight() && octreeIterator.moveRight());

                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                octreeIteratorValide.moveUp();
                octreeIteratorValide.gotoLeft();
            }
        }
        {
            FSize NbPart = 0;
            FSize NbLeafs = 0;
            {
                typename OctreeClass::Iterator octreeIterator(valideTree);
                octreeIterator.gotoBottomLeft();
                do{
                    NbPart += octreeIterator.getCurrentListSrc()->getNbParticles();
                    ++NbLeafs;
                } while(octreeIterator.moveRight());
            }
            {
                typename OctreeClass::Iterator octreeIterator(badTree);
                octreeIterator.gotoBottomLeft();

                do {
                    const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

                    ContainerClass* container = (octreeIterator.getCurrentListTargets());
                    const long long int*const dataDown = container->getDataDown();

                    for(FSize idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                        uassert((!isUsingTsm && dataDown[idxPart] == NbPart - 1) || (isUsingTsm && dataDown[idxPart] == NbPart));
                    }
                } while( octreeIterator.moveRight());
            }
        }
        {
            {
                // Check that each particle has been summed with all other
                typename OctreeClass::Iterator octreeIterator(badTree);
                octreeIterator.gotoBottomLeft();

                do {
                    uassert(octreeIterator.getCurrentListSrc()->getNbParticles() == octreeIterator.getCurrentCell()->getDataUp());
                } while( octreeIterator.moveRight() );
            }
        }
        {
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator valideOctreeIterator(valideTree);
            valideOctreeIterator.gotoBottomLeft();
            while(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                valideOctreeIterator.moveRight();
            }

            do {
                if(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                    uassert(false);
                    break;
                }

                if(octreeIterator.getCurrentListTargets()->getNbParticles() != valideOctreeIterator.getCurrentListTargets()->getNbParticles()){
                    uassert(false);
                }
                else {
                    ContainerClass* container = (octreeIterator.getCurrentListTargets());
                    const long long int*const dataDown = container->getDataDown();

                    ContainerClass* containerValide = (valideOctreeIterator.getCurrentListTargets());
                    const long long int*const dataDownValide = containerValide->getDataDown();

                    for(FSize idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                        uassert(dataDown[idxPart] == dataDownValide[idxPart]);
                    }
                }

            }while( octreeIterator.moveRight() && valideOctreeIterator.moveRight());
        }
    }

    typedef double FReal;
    typedef FTestCell                  CellClass;
    typedef FTestParticleContainer<FReal>     ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    typedef FFmmAlgorithmThreadProc<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassProc;


    void TestAlgo(){
        const int NbLevels = 7;
        const int SizeSubLevels = 3;
        const char* const filename = "../../Data/unitCubeXYZQ20k.bfma";
        FMpiFmaGenericLoader<FReal> loader(filename,app.global());

        OctreeClass realTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

        if( app.global().processCount() != 1){
            struct TestParticle{
                FPoint<FReal> position;
                const FPoint<FReal>& getPosition(){
                    return position;
                }
            };

            TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
            memset(particles, 0, sizeof(TestParticle) * loader.getNumberOfParticles());
            FReal physicalValue;
            for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(&particles[idxPart].position,&physicalValue);
            }

            FVector<TestParticle> finalParticles;
            FLeafBalance balancer;
            FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                        loader.getMyNumberOfParticles(),
                                                                        realTree.getBoxCenter(),
                                                                        realTree.getBoxWidth(),realTree.getHeight(),
                                                                        &finalParticles, &balancer);
            for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
                realTree.insert(finalParticles[idx].position);
            }

            delete[] particles;
        }
        else{
            FPoint<FReal> position;
            FReal physicalValue;
            const FSize nbParticles = loader.getNumberOfParticles();
            for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
                loader.fillParticle(&position,&physicalValue);
                realTree.insert(position);
            }
        }

        OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
        {
            FFmaGenericLoader<FReal> loaderSeq(filename);
            FPoint<FReal> position;
            FReal physicalValue;
            for(FSize idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
                loaderSeq.fillParticle(&position,&physicalValue);
                treeValide.insert(position);
            }
        }

        ValidateTree(realTree, treeValide);


        KernelClass kernels;
        FmmClassProc algo(app.global(),&realTree,&kernels);
        algo.execute();

        FmmClass algoValide(&treeValide,&kernels);
        algoValide.execute();

        ValidateFMMAlgoProc<OctreeClass,ContainerClass, FmmClassProc>(&realTree,&treeValide,&algo);
    }


    // set test
    void SetTests(){
        AddTest(&TestFmmAlgoProc::TestAlgo,"Test Algorithm");
    }
public:
    TestFmmAlgoProc(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }
};

// You must do this
TestClassMpi(TestFmmAlgoProc)


