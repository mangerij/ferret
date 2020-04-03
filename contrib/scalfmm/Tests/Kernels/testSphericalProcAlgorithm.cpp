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

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMpi.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include <iostream>

#include <cstdio>
#include <cstdlib>

// Uncoment to validate the FMM
#define VALIDATE_FMM

/** This program show an example of use of
  * the fmm basic algo it also check that eachh particles is little or longer
  * related that each other
  */

typedef double FReal;

#ifdef VALIDATE_FMM

static const FReal Epsilon = FReal(0.0005);

///////////////////////////////////////////////////////
// to test equality between good and potentialy bad solution
///////////////////////////////////////////////////////
/** To compare data */
template <class CellClass>
bool isEqualPole(const CellClass& me, const CellClass& other, FReal*const cumul){
    FMath::FAccurater<FReal> accurate;
    for(int idx = 0; idx < CellClass::GetPoleSize(); ++idx){
        accurate.add(me.getMultipole()[idx].getImag(),other.getMultipole()[idx].getImag());
        accurate.add(me.getMultipole()[idx].getReal(),other.getMultipole()[idx].getReal());
    }
    *cumul = accurate.getInfNorm()+ accurate.getL2Norm();
    return accurate.getInfNorm() < Epsilon && accurate.getL2Norm() < Epsilon;//FMath::LookEqual(cumul,FReal(0.0));
}

/** To compare data */
bool isEqualLocal(const FSphericalCell<FReal>& me, const FSphericalCell<FReal>& other,FReal*const cumul){
    FMath::FAccurater<FReal> accurate;
    for(int idx = 0; idx < FSphericalCell<FReal>::GetLocalSize(); ++idx){
        accurate.add(me.getLocal()[idx].getImag(),other.getLocal()[idx].getImag());
        accurate.add(me.getLocal()[idx].getReal(),other.getLocal()[idx].getReal());
    }
    *cumul = accurate.getInfNorm()+ accurate.getL2Norm();
    return accurate.getInfNorm() < Epsilon && accurate.getL2Norm() < Epsilon;//FMath::LookEqual(cumul,FReal(0.0));
}


template<class OctreeClass, class ContainerClass>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree){
    std::cout << "Check Result\n";
    {
        const int OctreeHeight = valideTree->getHeight();
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 1 ; --level){
            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                octreeIteratorValide.moveRight();
            }

            do {
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    FReal cumul;
                    if( !isEqualPole(*octreeIterator.getCurrentCell(),*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Pole Data are different. Cumul " << cumul << " at level " << level << " index is " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                    if( !isEqualLocal(*octreeIterator.getCurrentCell(),*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Local Data are different. Cumul " << cumul << " at level " << level << " index is " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                }

            } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    {
        // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            octreeIteratorValide.moveRight();
        }

        do {

            if( octreeIterator.getCurrentListSrc()->getNbParticles() != octreeIteratorValide.getCurrentListSrc()->getNbParticles()){
                std::cout << " Particules numbers is different " << std::endl;
            }
            if( octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                std::cout << " Index are differents " << std::endl;
            }

            ContainerClass* firstLeaf = octreeIterator.getCurrentListTargets();
            ContainerClass* valideLeaf = octreeIteratorValide.getCurrentListTargets();

            const FReal*const potentials = firstLeaf->getPotentials();
            const FReal*const forcesX = firstLeaf->getForcesX();
            const FReal*const forcesY = firstLeaf->getForcesY();
            const FReal*const forcesZ = firstLeaf->getForcesZ();
            const FReal*const positionX = firstLeaf->getPositions()[0];
            const FReal*const positionY = firstLeaf->getPositions()[1];
            const FReal*const positionZ = firstLeaf->getPositions()[2];
            const FReal*const validePositionX = valideLeaf->getPositions()[0];
            const FReal*const validePositionY = valideLeaf->getPositions()[1];
            const FReal*const validePositionZ = valideLeaf->getPositions()[2];
            const FReal*const validePotentials = valideLeaf->getPotentials();
            const FReal*const valideForcesX = valideLeaf->getForcesX();
            const FReal*const valideForcesY = valideLeaf->getForcesY();
            const FReal*const valideForcesZ = valideLeaf->getForcesZ();

            for(FSize idxLeaf = 0 ; idxLeaf < firstLeaf->getNbParticles() ; ++idxLeaf){

                int idxValideLeaf = 0;
                for(; idxValideLeaf < valideLeaf->getNbParticles() ; ++idxValideLeaf){
                    if( FMath::LookEqual(validePositionX[idxValideLeaf],positionX[idxLeaf]) &&
                            FMath::LookEqual(validePositionY[idxValideLeaf],positionY[idxLeaf]) &&
                            FMath::LookEqual(validePositionZ[idxValideLeaf],positionZ[idxLeaf]) ){
                        break;
                    }
                }

                if( idxValideLeaf < valideLeaf->getNbParticles() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    bool error = false;
                    if( FMath::RelatifDiff(validePotentials[idxValideLeaf] , potentials[idxLeaf])  > Epsilon ){
                        std::cout << " Potential error : " << validePotentials[idxValideLeaf] << " " << potentials[idxLeaf] << "\n";
                        error = true;
                    }
                    if( FMath::RelatifDiff(valideForcesX[idxValideLeaf],forcesX[idxLeaf]) > Epsilon
                            || FMath::RelatifDiff(valideForcesY[idxValideLeaf],forcesY[idxLeaf]) > Epsilon
                            || FMath::RelatifDiff(valideForcesZ[idxValideLeaf],forcesZ[idxLeaf]) > Epsilon){
                        std::cout << " Forces error : x " << valideForcesX[idxValideLeaf] << " " << forcesX[idxLeaf]
                                     << " y " << valideForcesY[idxValideLeaf]  << " " << forcesY[idxLeaf]
                                        << " z " << valideForcesZ[idxValideLeaf]  << " " << forcesZ[idxLeaf] << "\n";
                        error = true;
                    }
                    if( error ){
                        std::cout << "At position " << FPoint<FReal>(validePositionX[idxValideLeaf],validePositionY[idxValideLeaf],validePositionZ[idxValideLeaf])
                                  << " == " << FPoint<FReal>(positionX[idxLeaf],positionY[idxLeaf],positionZ[idxLeaf]) << std::endl;
                    }
                }
                else{
                    std::cout << "Particle not found " << FPoint<FReal>(positionX[idxLeaf],positionY[idxLeaf],positionZ[idxLeaf]) << std::endl;
                }
            }

        } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());
    }

    std::cout << "Done\n";
}
#endif


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Run a Spherical Harmonic (Old Implementation) FMM kernel with mpir parallelization.\n"
                         "The input file should be in binary format to enable distributed access.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::SHDevelopment);

    typedef FSphericalCell<FReal>         CellClass;
    typedef FP2PParticleContainer<FReal>         ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<FReal, CellClass, ContainerClass >     KernelClass;

    typedef FFmmAlgorithmThreadProc<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithmThread<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClassNoProc;


    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int DevP = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 8);
    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    FTic counter;
    const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
                "../Data/test20k.bin.fma.single":
                "../Data/test20k.bin.fma.double";
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFilename);

    std::cout << "Opening : " << filename << "\n";

    FMpiFmaGenericLoader<FReal> loader(filename, app.global());
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    CellClass::Init(DevP);


    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;
        counter.tic();

        struct TestParticle{
            FPoint<FReal> position;
            FReal physicalValue;
            const FPoint<FReal>& getPosition(){
                return position;
            }
        };

        TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
        memset(particles, 0, sizeof(TestParticle) * loader.getNumberOfParticles());

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
        }

        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;
        // FMpiTreeBuilder< FReal,TestParticle >::ArrayToTree(app.global(), particles, loader.getNumberOfParticles(),
        // 						 tree.getBoxCenter(),
        // 						 tree.getBoxWidth(),
        // 						 tree.getHeight(), &finalParticles,&balancer);
        FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                    loader.getMyNumberOfParticles(),
                                                                    tree.getBoxCenter(),
                                                                    tree.getBoxWidth(),tree.getHeight(),
                                                                    &finalParticles, &balancer);

        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            tree.insert(finalParticles[idx].position,finalParticles[idx].physicalValue);

        }

        delete[] particles;

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        FPoint<FReal> position;
        FReal physicalValue;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&position,&physicalValue);
            tree.insert(position, physicalValue);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    std::cout << "Create kernel..." << std::endl;

    KernelClass kernels(DevP, NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());

    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    FmmClass algo(app.global(),&tree,&kernels);

    counter.tic();
    algo.execute();
    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential

        FReal potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "My potential is " << potential << std::endl;

        potential = app.global().reduceSum(potential);
        fx = app.global().reduceSum(fx);
        fy = app.global().reduceSum(fy);
        fz = app.global().reduceSum(fz);


        if(app.global().processId() == 0){
            std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
            std::cout << "Potential Sum = " << potential << std::endl;
        }
    }

#ifdef VALIDATE_FMM
    {
        OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
        {
            FFmaGenericLoader<FReal> loaderSeq(filename);
            FPoint<FReal> position;
            FReal physicalValue;
            for(FSize idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
                loaderSeq.fillParticle(&position,&physicalValue);
                treeValide.insert(position,physicalValue);
            }
        }

        std::cout << "Working on particles ..." << std::endl;
        FmmClassNoProc algoValide(&treeValide,&kernels);
        counter.tic();
        algoValide.execute();
        counter.tac();
        std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

        FReal potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
        std::cout << "Potential = " << potential << std::endl;

        ValidateFMMAlgoProc<OctreeClass,ContainerClass>(&tree,&treeValide);
    }
#endif


    // -----------------------------------------------------

    return 0;
}



