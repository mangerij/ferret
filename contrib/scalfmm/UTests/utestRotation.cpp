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

#include "Utils/FGlobal.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FTemplate.hpp"
#include "Utils/FMath.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Kernels/Rotation/FRotationCell.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"


/** the test class
 *
 */
class TestRotationDirect : public FUTester<TestRotationDirect> {
    /** The test method to factorize all the test based on different kernels */
    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class LeafClass,
              class OctreeClass, class FmmClass, int ORDER>
    void RunTest(){
        //
        // Load particles
        //
        if(sizeof(FReal) == sizeof(float) ) {
            std::cerr << "No input data available for Float "<< std::endl;
            exit(EXIT_FAILURE);
        }
        const std::string parFile( (sizeof(FReal) == sizeof(float))?
                                       "Test/DirectFloat.bfma":
                                       "UTest/DirectDouble.bfma");
        //
        std::string filename(SCALFMMDataPath+parFile);
        //
        FFmaGenericLoader<FReal> loader(filename);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            uassert(false);
            return;
        }
        //
        FSize nbParticles = loader.getNumberOfParticles() ;
        Print("Number of particles:");
        Print(nbParticles);

        std::cout << "Approximation order: " << ORDER << std::endl;

        const int NbLevels      = 4;
        const int SizeSubLevels = 2;
        //
        FmaRWParticle<FReal, 8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];

        loader.fillParticle(particles,nbParticles);
        //
        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        //   Insert particle in the tree
        //
        FReal  sum = 0.0, a= 0.0;

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> position(particles[idxPart].getPosition() );
            //
            tree.insert(particles[idxPart].getPosition() , idxPart, particles[idxPart].getPhysicalValue());
            //
            sum += FMath::Abs(particles[idxPart].getPhysicalValue()) ;
            a = std::max(a,position.getX()*position.getX()+position.getY()*position.getY()+position.getZ()*position.getZ());

        }
        double CorErr = sum/a;


        // Run FMM
        Print("Fmm...");
        //KernelClass kernels(NbLevels,loader.getBoxWidth());
        KernelClass* kernels = new KernelClass(NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&tree,kernels);
        FTic timer;
        algo.execute();
        timer.tac();
        std::cout << "Computation Time: " << timer.elapsed() << " for order "<<  ORDER <<std::endl;
        //
        FReal energy= 0.0 , energyD = 0.0 ;
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Compute direct energy
        /////////////////////////////////////////////////////////////////////////////////////////////////

        for(FSize idx = 0 ; idx < loader.getNumberOfParticles()  ; ++idx){
            energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Compare
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Compute Diff...");
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
        { // Check that each particle has been summed with all other

            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                    potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                    fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                    fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                    fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                    energy   += potentials[idxPart]*physicalValues[idxPart];
                }
            });
        }
        delete kernels;
        delete[] particles;

        // Print for information

        Print("Potential diff is = ");
        printf("         Pot L2Norm     %e\n",potentialDiff.getL2Norm());
        printf("         Pot RL2Norm   %e\n",potentialDiff.getRelativeL2Norm());
        printf("         Pot RMSError   %e\n",potentialDiff.getRMSError());
        Print("Fx diff is = ");
        printf("         Fx L2Norm     %e\n",fx.getL2Norm());
        printf("         Fx RL2Norm   %e\n",fx.getRelativeL2Norm());
        printf("         Fx RMSError   %e\n",fx.getRMSError());
        Print("Fy diff is = ");
        printf("        Fy L2Norm     %e\n",fy.getL2Norm());
        printf("        Fy RL2Norm   %e\n",fy.getRelativeL2Norm());
        printf("        Fy RMSError   %e\n",fy.getRMSError());
        Print("Fz diff is = ");
        printf("        Fz L2Norm     %e\n",fz.getL2Norm());
        printf("        Fz RL2Norm   %e\n",fz.getRelativeL2Norm());
        printf("        Fz RMSError   %e\n",fz.getRMSError());
        FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
        printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
        printf("  Energy Error  =   %.12e\n",FMath::Abs(energy-energyD));
        printf("  Energy FMM    =   %.12e\n",FMath::Abs(energy));
        printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyD));

        // Assert
        //	const FReal MaximumDiffPotential = FReal(9e-3);
        //	const FReal MaximumDiffForces     = FReal(9e-2);
        double epsilon = 1.0/FMath::pow2(ORDER);
        const FReal MaximumDiffPotential = FReal(CorErr*epsilon);
        const FReal MaximumDiffForces     = FReal(10*CorErr*epsilon);
        printf(" Criteria error - Epsilon  %e  \n",epsilon);

        Print("Test1 - Error Relative L2 norm Potential ");
        uassert(potentialDiff.getRelativeL2Norm() < MaximumDiffPotential);    //1
        Print("Test2 - Error RMS L2 norm Potential ");
        FReal CoerrRMS = potentialDiff.getl2Dot()/FMath::Sqrt(static_cast<FReal>(nbParticles));

        uassert(potentialDiff.getRMSError() < CoerrRMS*MaximumDiffPotential);  //2
        Print("Test3 - Error Relative L2 norm FX ");
        uassert(fx.getRelativeL2Norm()  < MaximumDiffForces);
        CoerrRMS = fx.getl2Dot()/FMath::Sqrt(static_cast<FReal>(nbParticles));
        //3
        Print("Test4 - Error RMS L2 norm FX ");
        uassert(fx.getRMSError() < CoerrRMS*MaximumDiffForces);                      //4
        Print("Test5 - Error Relative L2 norm FY ");
        uassert(fy.getRelativeL2Norm()  < MaximumDiffForces);                       //5
        Print("Test6 - Error RMS L2 norm FY ");
        CoerrRMS = fy.getl2Dot()/FMath::Sqrt(static_cast<FReal>(nbParticles));
        uassert(fy.getRMSError() < CoerrRMS*MaximumDiffForces);                      //6
        Print("Test7 - Error Relative L2 norm FZ ");
        uassert(fz.getRelativeL2Norm()  < MaximumDiffForces);                      //8
        Print("Test8 - Error RMS L2 norm FZ ");
        CoerrRMS = fz.getl2Dot()/FMath::Sqrt(static_cast<FReal>(nbParticles));
        uassert(fz.getRMSError() < CoerrRMS*MaximumDiffForces);                                           //8
        Print("Test9 - Error Relative L2 norm F ");
        uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
        Print("Test10 - Relative error Energy ");
        uassert(FMath::Abs(energy-energyD) /FMath::Abs(energyD)< MaximumDiffPotential);                     //10  Total Energy


    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }

    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    //    static const int P = 9;

    /** Rotation */
    template<int P>
    void TestRotation(){
        typedef double FReal;
        typedef FRotationCell<FReal,P>              CellClass;
        typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

        typedef FRotationKernel<FReal,CellClass, ContainerClass, P >          KernelClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<FReal,CellClass, ContainerClass, KernelClass, LeafClass, OctreeClass, FmmClass,P>();
    }

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        FForAllThis::For<int, 6, 40, 2, TestRotationDirect>(this);
    }

public:
    template <int P>
    void For(){
        AddTest(&TestRotationDirect::TestRotation<P>,"Test Rotation Kernel");
    }
};


// You must do this
TestClass(TestRotationDirect)
