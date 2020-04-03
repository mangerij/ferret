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
// ================
#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Containers/FOctree.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Core/FFmmAlgorithmThread.hpp"

#include "FUTester.hpp"

#include "Components/FSimpleLeaf.hpp"


#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
/*
  In this test we compare the Chebyschev fmm results and the direct results.
 */


/** the test class
 *
 */
class TestChebyshevDirect : public FUTester<TestChebyshevDirect> {

    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
              class LeafClass, class OctreeClass, class FmmClass>
    void RunTest()	{
        //
#ifdef _OPENMP
        std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
        std::cout << "\n>> OpenMP test !!!\n" << std::endl;
        exit(EXIT_FAILURE);
#endif
        //
        // Load particles
        //
        if(sizeof(FReal) == sizeof(float) ) {
            std::cerr << "No input data available for Float "<< std::endl;
            exit(EXIT_FAILURE);
        }
        const std::string parFile( (sizeof(FReal) == sizeof(float))?
                                       "Test/DirectFloatbfma":
                                       "UTest/DirectDouble.bfma");
        //
        std::string filename(SCALFMMDataPath+parFile);
        //
        FFmaGenericLoader<FReal> loader(filename);
        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels        = 4;
        const int SizeSubLevels = 2;

        // Create Matrix Kernel
        const MatrixKernelClass MatrixKernel; // FUKernelTester is only designed to work with 1/R, i.e. matrix kernel ctor takes no argument.

        // Load particles
        FSize nbParticles = loader.getNumberOfParticles() ;
        FmaRWParticle<FReal, 8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];

        loader.fillParticle(particles,nbParticles);

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        //
        //   Insert particle in the tree
        //
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            tree.insert(particles[idxPart].getPosition() , idxPart, particles[idxPart].getPhysicalValue() );
        }
        //
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run FMM computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Fmm...");
        KernelClass kernels(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        FmmClass algo(&tree,&kernels);
        algo.execute();
        //0
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
        const FReal MaximumDiffPotential = FReal(9e-3);
        const FReal MaximumDiffForces     = FReal(9e-2);

        Print("Test1 - Error Relative L2 norm Potential ");
        uassert(potentialDiff.getRelativeL2Norm() < MaximumDiffPotential);    //1
        Print("Test2 - Error RMS L2 norm Potential ");
        uassert(potentialDiff.getRMSError() < MaximumDiffPotential);  //2
        Print("Test3 - Error Relative L2 norm FX ");
        uassert(fx.getRelativeL2Norm()  < MaximumDiffForces);                       //3
        Print("Test4 - Error RMS L2 norm FX ");
        uassert(fx.getRMSError() < MaximumDiffForces);                      //4
        Print("Test5 - Error Relative L2 norm FY ");
        uassert(fy.getRelativeL2Norm()  < MaximumDiffForces);                       //5
        Print("Test6 - Error RMS L2 norm FY ");
        uassert(fy.getRMSError() < MaximumDiffForces);                      //6
        Print("Test7 - Error Relative L2 norm FZ ");
        uassert(fz.getRelativeL2Norm()  < MaximumDiffForces);                      //8
        Print("Test8 - Error RMS L2 norm FZ ");
        uassert(fz.getRMSError() < MaximumDiffForces);                                           //8
        Print("Test9 - Error Relative L2 norm F ");
        uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
        Print("Test10 - Relative error Energy ");
        uassert(FMath::Abs(energy-energyD) /energyD< MaximumDiffPotential);                     //10  Total Energy


    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                      << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated()
                      << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated()
                      << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }


    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////


    /** TestChebKernel */
    void TestChebKernel(){
        typedef double FReal;
        const unsigned int ORDER = 6;
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>();
    }

    /** TestChebSymKernel */
    void TestChebSymKernel(){
        typedef double FReal;
        const unsigned int ORDER = 6;
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>();
    }



    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestChebyshevDirect::TestChebKernel,"Test Chebyshev Kernel with one big SVD");
        AddTest(&TestChebyshevDirect::TestChebSymKernel,"Test Chebyshev Kernel with 16 small SVDs and symmetries");
    }
};


// You must do this
TestClass(TestChebyshevDirect)




