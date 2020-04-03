
// ===================================================================================
// Copyright ScalFmm 2011 INRIA,
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

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Kernels/Rotation/FRotationCell.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Components/FTypedLeaf.hpp"
#include "Extensions/FExtendCellType.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FFmaGenericLoader.hpp"

#include "Core/FFmmAlgorithmThreadTsm.hpp"
#include "Core/FFmmAlgorithmTsm.hpp"

#include "FUTester.hpp"


/** the test class the rotation and target source model.
 *
 */
class TestRotationDirectTsm : public FUTester<TestRotationDirectTsm> {
	/** The test method to factorize all the test based on different kernels */
    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class LeafClass,
	class OctreeClass, class FmmClass>
	void RunTest(){
		// Warning in make test the exec dir it Build/UTests
		// Load particles
		const int nbSources = 5000;
		const int nbTargets = 5000;

		FRandomLoader<FReal> loader(nbSources + nbTargets);

		Print("Number of particles:");
		Print(loader.getNumberOfParticles());

		const int NbLevels      = 4;
		const int SizeSubLevels = 3;

		// Create octree
		OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

		const FReal physicalValue = 0.10;
		//
		FmaRWParticle<FReal, 8,8>* const particlesTargets = new FmaRWParticle<FReal, 8,8>[nbTargets];
		for(FSize idxPart = 0 ; idxPart < nbTargets ; ++idxPart){
            FPoint<FReal> position;
			loader.fillParticle(&position);
			// put in tree
			tree.insert(position, FParticleTypeTarget, idxPart, physicalValue);
			// get copy
			particlesTargets[idxPart].setPosition(position);
			*(particlesTargets[idxPart].setPhysicalValue()) = physicalValue;
			*(particlesTargets[idxPart].setPotential())        = 0.0;
			particlesTargets[idxPart].setForces()[0]        = 0.0;
			particlesTargets[idxPart].setForces()[1]        = 0.0;
			particlesTargets[idxPart].setForces()[2]        = 0.0;
		}

		FmaRWParticle<FReal, 8,8>* const particlesSources = new FmaRWParticle<FReal, 8,8>[nbSources];
		for(FSize idxPart = 0 ; idxPart < nbSources ; ++idxPart){
            FPoint<FReal> position;
			loader.fillParticle(&position);
			// put in tree
			tree.insert(position, FParticleTypeSource, idxPart, physicalValue);
			// get copy
			particlesSources[idxPart].setPosition(position);
			*(particlesSources[idxPart].setPhysicalValue()) = physicalValue;
		}


		// Run FMM
		Print("Fmm...");
		//KernelClass kernels(NbLevels,loader.getBoxWidth());
		KernelClass kernels(NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
		FmmClass algo(&tree,&kernels);
		algo.execute();
		//
		//

		// Run direct computation
        Print("Direct...");
		for(int idxTarget = 0 ; idxTarget < nbTargets ; ++idxTarget){
			for(int idxOther = 0 ; idxOther < nbSources ; ++idxOther){
                FP2PR::NonMutualParticles(
						particlesSources[idxOther].getPosition().getX(), particlesSources[idxOther].getPosition().getY(),
						particlesSources[idxOther].getPosition().getZ(),particlesSources[idxOther].getPhysicalValue(),
						particlesTargets[idxTarget].getPosition().getX(), particlesTargets[idxTarget].getPosition().getY(),
						particlesTargets[idxTarget].getPosition().getZ(),particlesTargets[idxTarget].getPhysicalValue(),
						&particlesTargets[idxTarget].setForces()[0],&particlesTargets[idxTarget].setForces()[1],
                        &particlesTargets[idxTarget].setForces()[2],particlesTargets[idxTarget].setPotential());
			}
		}

		//
		// Assert
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Compute direct energy
		/////////////////////////////////////////////////////////////////////////////////////////////////
		FReal energy= 0.0 , energyD = 0.0 ;
		for(int idx = 0 ; idx <  nbTargets  ; ++idx){
			energyD +=  particlesTargets[idx].getPotential()*particlesTargets[idx].getPhysicalValue() ;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Compare
		/////////////////////////////////////////////////////////////////////////////////////////////////
		Print("Compute Diff...");
		FMath::FAccurater<FReal> potentialDiff;
		FMath::FAccurater<FReal> fx, fy, fz;
		{ // Check that each particle has been summed with all other

			tree.forEachLeaf([&](LeafClass* leaf){
				if( leaf->getTargets()->getNbParticles() ){
					const FReal*const potentials = leaf->getTargets()->getPotentials();
					const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
					const FReal*const forcesX = leaf->getTargets()->getForcesX();
					const FReal*const forcesY = leaf->getTargets()->getForcesY();
					const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
					const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
					const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

					for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
						const FSize indexPartOrig = indexes[idxPart];
						potentialDiff.add(particlesTargets[indexPartOrig].getPotential(),potentials[idxPart]);
						fx.add(particlesTargets[indexPartOrig].getForces()[0],forcesX[idxPart]);
						fy.add(particlesTargets[indexPartOrig].getForces()[1],forcesY[idxPart]);
						fz.add(particlesTargets[indexPartOrig].getForces()[2],forcesZ[idxPart]);
						energy   += potentials[idxPart]*physicalValues[idxPart];
					}
				}
			});
		}

		delete[] particlesTargets;
		delete[] particlesSources;
		//
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
		//
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
			std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
			std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
			std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
		}
	}

	///////////////////////////////////////////////////////////
	// The tests!
	///////////////////////////////////////////////////////////

	static const int P = 9;

	/** Rotation */
	void TestRotation(){
        typedef double FReal;
        typedef FTypedRotationCell<FReal,P>    CellClass;
		typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

        typedef FRotationKernel<FReal, CellClass, ContainerClass, P >          KernelClass;

        typedef FTypedLeaf<FReal,ContainerClass >                     LeafClass;
		typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

		typedef FFmmAlgorithmTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<FReal, CellClass, ContainerClass, KernelClass, LeafClass, OctreeClass, FmmClass>();
	}

	void TestRotationThread(){
        typedef double FReal;
        typedef FTypedRotationCell<FReal,P>    CellClass;
		typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

        typedef FRotationKernel<FReal, CellClass, ContainerClass, P >          KernelClass;

        typedef FTypedLeaf<FReal,ContainerClass >                     LeafClass;
		typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

		typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<FReal, CellClass, ContainerClass, KernelClass, LeafClass, OctreeClass, FmmClass>();
	}

	///////////////////////////////////////////////////////////
	// Set the tests!
	///////////////////////////////////////////////////////////

	/** set test */
	void SetTests(){
		AddTest(&TestRotationDirectTsm::TestRotation,"Test Rotation Kernel TSM");
		AddTest(&TestRotationDirectTsm::TestRotationThread,"Test Rotation Kernel TSM thread");
	}
};


// You must do this
TestClass(TestRotationDirectTsm)



