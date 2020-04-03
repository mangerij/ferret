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
// @FUSE_BLAS
// @FUSE_MPI
// ================


#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Containers/FOctree.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FMpiFmaGenericLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "Core/FFmmAlgorithmThreadProc.hpp"
#include "BalanceTree/FLeafBalance.hpp"

#include "FUTester.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/*
 * In this test, we compare the results of the Chebyshev FMM over
 * multiple processus and the direct results.
 */


class TestChebyshevMpiDirect : public FUTesterMpi<TestChebyshevMpiDirect>{

    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
              class LeafClass, class OctreeClass, class FmmClassProc>
    void RunTest(){
        const std::string parFile( (sizeof(FReal) == sizeof(float))?
                                       "Test/DirectFloatbfma":
                                       "UTest/DirectDouble.bfma");
        std::string filename(SCALFMMDataPath+parFile);
        //    std::string filename("./sphere120Solved.bfma");

        FMpiFmaGenericLoader<FReal> loader(filename,app.global());
        Print("Number of particles :");
        Print(loader.getNumberOfParticles());

        const int nbLevels = 4;
        const int sizeOfSubLevel = 2;

        // Create Matrix Kernel
        const MatrixKernelClass MatrixKernel; // FUKernelTester is only designed to work with 1/R, i.e. matrix kernel ctor takes no argument.

        // Create octree

        struct TestParticle : public FmaRWParticle<FReal, 8,8>{
            FSize index;
            // const FPoint<FReal>& getPosition(){
            // 	return position;
            // }
        };

        FSize nbParticles = loader.getMyNumberOfParticles();
        TestParticle* const particles = new TestParticle[nbParticles];
        memset(particles,0,sizeof(TestParticle)*nbParticles);

        //idx (in file) of the first part that will be used by this proc.
        FSize idxStart = loader.getStart();

        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            //Storage of the index (in the original file) of each part.
            particles[idxPart].index = idxPart + idxStart;
            // Read particles from file
            loader.fillParticle(particles[idxPart]);
        }

        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;
        OctreeClass tree(nbLevels,sizeOfSubLevel,loader.getBoxWidth(),loader.getCenterOfBox());
        // FMpiTreeBuilder< FReal,TestParticle >::ArrayToTree(app.global(), particles, loader.getMyNumberOfParticles(),
        // 						 tree.getBoxCenter(),
        // 						 tree.getBoxWidth(),
        // 						 tree.getHeight(), &finalParticles,&balancer);
        FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                    loader.getMyNumberOfParticles(),
                                                                    tree.getBoxCenter(),
                                                                    tree.getBoxWidth(),tree.getHeight(),
                                                                    &finalParticles, &balancer);
        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            tree.insert(finalParticles[idx].getPosition(),int(finalParticles[idx].index),finalParticles[idx].getPhysicalValue());
        }


        KernelClass* kernels= new KernelClass(nbLevels, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        FmmClassProc algorithm(app.global(),&tree, kernels);
        algorithm.execute();

        //Check datas
        {
            Print("Comput Differences with direct computation\n");
            FMath::FAccurater<FReal> potentialDiff;
            FMath::FAccurater<FReal> fx, fy, fz;
            FReal energy = 0.0;
            FReal * datas = new FReal[loader.getNbRecordPerline()];
            memset(datas,0,loader.getNbRecordPerline()*sizeof(FReal));
            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                const FReal*const physicalValues    = leaf->getTargets()->getPhysicalValues();
                const FReal*const forcesX           = leaf->getTargets()->getForcesX();
                const FReal*const forcesY           = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ           = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf         = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = FSize(indexes[idxPart])-idxStart;
                    //It's a proc on my left that used to keep this part
                    if(indexPartOrig < 0){
                        loader.fill1Particle(datas,indexes[idxPart]);
                        potentialDiff.add(datas[4],potentials[idxPart]);
                        fx.add(datas[5],forcesX[idxPart]);
                        fy.add(datas[6],forcesY[idxPart]);
                        fz.add(datas[7],forcesZ[idxPart]);
                        energy   += potentials[idxPart]*physicalValues[idxPart];
                        memset(datas,0,loader.getNbRecordPerline()*sizeof(FReal));
                    }
                    else{
                        //It's a proc on my right that used to keep this part
                        if(indexPartOrig >= loader.getMyNumberOfParticles()){
                            loader.fill1Particle(datas,FSize(indexes[idxPart]));
                            potentialDiff.add(datas[4],potentials[idxPart]);
                            fx.add(datas[5],forcesX[idxPart]);
                            fy.add(datas[6],forcesY[idxPart]);
                            fz.add(datas[7],forcesZ[idxPart]);
                            energy   += potentials[idxPart]*physicalValues[idxPart];
                            if(datas[0] != leaf->getTargets()->getPositions()[0][idxPart]){
                                printf("- %d - Problem %lld !! \t [%e,%e,%e,%e] [%e,%e,%e,%e] \n",
                                       app.global().processId(),indexes[idxPart],
                                       datas[0],datas[1],datas[2],datas[3],
                                        potentials[idxPart],forcesX[idxPart],
                                        forcesY[idxPart],forcesZ[idxPart]);
                            }
                            memset(datas,0,loader.getNbRecordPerline()*sizeof(FReal));
                        }
                        //I already have this part
                        else{
                            potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                            fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                            fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                            fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                            energy   += potentials[idxPart]*physicalValues[idxPart];
                            //if(particles[indexPartOrig].getPosition().getX() != leaf->getTargets()->getPositions()[0][idxPart]){
                            // printf("%d - Problem %d !! [%e,%e,%e,%e] [%e,%e,%e,%e] \n",
                            //        app.global().processId(),indexPartOrig,particles[indexPartOrig].forces[0],particles[indexPartOrig].forces[1],
                            //        particles[indexPartOrig].forces[2],particles[indexPartOrig].potential,
                            //        forcesX[idxPart],forcesY[idxPart],
                            //        forcesZ[idxPart],potentials[idxPart]);
                            //}
                        }
                    }
                }
            });

            //Compute Direct Energy
            FReal energyD = 0.0;
            for(FSize idx = 0 ; idx <  loader.getMyNumberOfParticles()  ; ++idx){
                energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
            }

            FReal energyDTot = 0.0;
            FReal energyFMMTot = 0.0;
            //Double reduce !! could be made into one ...
            MPI_Reduce(&energyD,&energyDTot,1,FMpi::GetType(energyD),MPI_SUM,0,app.global().getComm());
            MPI_Reduce(&energy,&energyFMMTot,1,FMpi::GetType(energy),MPI_SUM,0,app.global().getComm());


            //Summarize
            FMath::FAccurater<FReal>* FXYZP = new FMath::FAccurater<FReal>[app.global().processCount()*4];
            FMath::FAccurater<FReal>* fxyzp = new FMath::FAccurater<FReal>[4];
            fxyzp[0] = fx;
            fxyzp[1] = fy;
            fxyzp[2] = fz;
            fxyzp[3] = potentialDiff;
            MPI_Gather(fxyzp,4*sizeof(FMath::FAccurater<FReal>),MPI_BYTE,FXYZP,4*sizeof(FMath::FAccurater<FReal>),MPI_BYTE,0,app.global().getComm());
            if(app.global().processId() == 0){
                for(int k=1 ; k<app.global().processCount(); ++k){
                    FXYZP[0].add(FXYZP[k*4+0]);
                    FXYZP[1].add(FXYZP[k*4+1]);
                    FXYZP[2].add(FXYZP[k*4+2]);
                    FXYZP[3].add(FXYZP[k*4+3]);
                }
                printf("Potential diff is = \n");
                printf("         Pot L2Norm     %e\n",FXYZP[3].getL2Norm());
                printf("         Pot RL2Norm    %e\n",FXYZP[3].getRelativeL2Norm());
                printf("         Pot RMSError   %e\n",FXYZP[3].getRMSError());
                printf("Fx diff is = \n");
                printf("         Fx L2Norm     %e\n",FXYZP[0].getL2Norm());
                printf("         Fx RL2Norm    %e\n",FXYZP[0].getRelativeL2Norm());
                printf("         Fx RMSError   %e\n",FXYZP[0].getRMSError());
                printf("Fy diff is = \n");
                printf("        Fy L2Norm     %e\n",FXYZP[1].getL2Norm());
                printf("        Fy RL2Norm    %e\n",FXYZP[1].getRelativeL2Norm());
                printf("        Fy RMSError   %e\n",FXYZP[1].getRMSError());
                printf("Fz diff is = \n");
                printf("        Fz L2Norm     %e\n",FXYZP[2].getL2Norm());
                printf("        Fz RL2Norm    %e\n",FXYZP[2].getRelativeL2Norm());
                printf("        Fz RMSError   %e\n",FXYZP[2].getRMSError());
                FReal L2error = (FXYZP[0].getRelativeL2Norm()*FXYZP[0].getRelativeL2Norm()
                        + FXYZP[1].getRelativeL2Norm()*FXYZP[1].getRelativeL2Norm()
                        + FXYZP[2].getRelativeL2Norm() *FXYZP[2].getRelativeL2Norm());

                printf("  Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
                printf("  Energy Error  =   %.12e\n",FMath::Abs(energyFMMTot-energyDTot));
                printf("  Energy FMM    =   %.12e\n",FMath::Abs(energyFMMTot));
                printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyDTot));

                // Assert
                const FReal MaximumDiffPotential = FReal(9e-3);
                const FReal MaximumDiffForces     = FReal(9e-2);

                Print("Test1 - Error Relative L2 norm Potential ");
                uassert(FXYZP[3].getRelativeL2Norm() < MaximumDiffPotential);    //1
                Print("Test2 - Error RMS L2 norm Potential ");
                uassert(FXYZP[3].getRMSError() < MaximumDiffPotential);  //2
                Print("Test3 - Error Relative L2 norm FX ");
                uassert(FXYZP[0].getRelativeL2Norm()  < MaximumDiffForces);                       //3
                Print("Test4 - Error RMS L2 norm FX ");
                uassert(FXYZP[0].getRMSError() < MaximumDiffForces);                      //4
                Print("Test5 - Error Relative L2 norm FY ");
                uassert(FXYZP[1].getRelativeL2Norm()  < MaximumDiffForces);                       //5
                Print("Test6 - Error RMS L2 norm FY ");
                uassert(FXYZP[1].getRMSError() < MaximumDiffForces);                      //6
                Print("Test7 - Error Relative L2 norm FZ ");
                uassert(FXYZP[2].getRelativeL2Norm()  < MaximumDiffForces);                      //8
                Print("Test8 - Error RMS L2 norm FZ ");
                uassert(FXYZP[2].getRMSError() < MaximumDiffForces);                                           //8
                Print("Test9 - Error Relative L2 norm F ");
                uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
                Print("Test10 - Relative error Energy ");
                uassert(FMath::Abs(energyFMMTot-energyDTot) /energyDTot < MaximumDiffPotential);                     //10  Total Energy

            }
        }
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << app.global().processId() << "-> Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                      << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << app.global().processId() << "-> Max memory used " << FMemStats::controler.getMaxAllocated()
                      << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << app.global().processId() << "-> Total memory used " << FMemStats::controler.getTotalAllocated()
                      << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }





    /** TestChebSymKernel */
    void TestChebSymKernel(){
        const unsigned int ORDER = 6;
        typedef double FReal;
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClassProc;
        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClassProc>();
    }

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestChebyshevMpiDirect::TestChebSymKernel,"Test Chebyshev Kernel with 16 small SVDs and symmetries");
    }

public:
    TestChebyshevMpiDirect(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }

};

TestClassMpi(TestChebyshevMpiDirect);
