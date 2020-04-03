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

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Kernels/Spherical/FSphericalCell.hpp"
#include "Kernels/Spherical/FSphericalKernel.hpp"
#include "Components/FSimpleLeaf.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Files/FTreeIO.hpp"

#include "Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/**
 * This test compare a previous FMM result with a previous simulation result.
 */

typedef double FReal;
typedef FSphericalCell<FReal>           CellClass;
typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

typedef FSphericalKernel< FReal, CellClass, ContainerClass >          KernelClass;

typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

/** To check if a value is correct */
bool IsSimilar(const FReal good, const FReal other){
	const FReal Epsilon = FReal(0.0001);
	return (FMath::Abs(good-other)/FMath::Abs(good)) < Epsilon;
}

/** The test class */
class TestSphericalWithPrevious : public FUTester<TestSphericalWithPrevious> {
	/** the test */
	void TestTree(){
		if(sizeof(FReal) == sizeof(float) ) {
			std::cerr << "No input data available for Float "<< std::endl;
			uassert(false);
		}
		//
		//  Load a Tree
		const std::string DataFile = (sizeof(FReal) == sizeof(float))?
				SCALFMMDataPath+"UTest/SphericalPrevious.data.single":
				SCALFMMDataPath+"UTest/SphericalPrevious.data.double";
		//
		// Load particles
		//

		const std::string parFile( (sizeof(FReal) == sizeof(float))?
				"UTest/DirectFloat.bfma":
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
		FSize nbParticles = loader.getNumberOfParticles() ;
		//
		const int NbLevels      = 5;
		const int SizeSubLevels = 3;
		const int DevP = 9;
		//
		// Create octree
		//
		FSphericalCell<FReal>::Init(DevP);
		//
		OctreeClass testTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
		//
		for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            FPoint<FReal> position;
			FReal physicalValue = 0.0;
			loader.fillParticle(&position,&physicalValue);
			// put in tree
			testTree.insert(position, idxPart, physicalValue);
		}
		//
		//  Run simulation 1
		KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
		FmmClass algo(&testTree,&kernels);
		Print("Run simulation 1 ...");

		algo.execute();

		// If needed save the result
        //FTreeIO<FReal>::Save<OctreeClass, CellClass, LeafClass, ContainerClass >(DataFile.c_str(), testTree);

		// Load previous result
		OctreeClass goodTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FTreeIO<FReal>::Load<OctreeClass, CellClass, LeafClass, ContainerClass >(DataFile.c_str(), goodTree);

		// Compare the two simulations
		Print("Check the particles...");
		{ // Check that each particle has been summed with all other
			OctreeClass::Iterator testOctreeIterator(&testTree);
			OctreeClass::Iterator goodOctreeIterator(&goodTree);

			testOctreeIterator.gotoBottomLeft();
			goodOctreeIterator.gotoBottomLeft();

			do{
				if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
					uassert(false);
					break;
				}

				if(testOctreeIterator.getCurrentListSrc()->getNbParticles() != goodOctreeIterator.getCurrentListSrc()->getNbParticles()){
					uassert(false);
					break;
				}

				const ContainerClass* testLeaf = testOctreeIterator.getCurrentListSrc();
				const ContainerClass* goodLeaf = goodOctreeIterator.getCurrentListSrc();

				for(FSize idxPart = 0 ; idxPart < testLeaf->getNbParticles() ; ++idxPart ){
					uassert( IsSimilar(goodLeaf->getPotentials()[idxPart], testLeaf->getPotentials()[idxPart]) );
					uassert( IsSimilar(goodLeaf->getForcesX()[idxPart], testLeaf->getForcesX()[idxPart]) );
					uassert( IsSimilar(goodLeaf->getForcesY()[idxPart], testLeaf->getForcesY()[idxPart]) );
					uassert( IsSimilar(goodLeaf->getForcesZ()[idxPart], testLeaf->getForcesZ()[idxPart]) );
				}

				if(!testOctreeIterator.moveRight()){
					if(goodOctreeIterator.moveRight()){
						uassert(false);
					}
					break;
				}
				if(!goodOctreeIterator.moveRight()){
					uassert(false);
					break;
				}

			} while(true);
		}
		Print("Check the leaves...");
		{ // Ceck if there is number of NbPart summed at level 1
			OctreeClass::Iterator testOctreeIterator(&testTree);
			OctreeClass::Iterator goodOctreeIterator(&goodTree);

			testOctreeIterator.gotoBottomLeft();
			goodOctreeIterator.gotoBottomLeft();

			for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
				do{
					if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
						uassert(false);
						break;
					}

					for(int idxLocal = 0 ; idxLocal < CellClass::GetLocalSize() ; ++idxLocal){
						IsSimilar(testOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getReal(),
								goodOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getReal());
						IsSimilar(testOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getImag(),
								goodOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getImag());
					}

					for(int idxPole = 0 ; idxPole < CellClass::GetPoleSize() ; ++idxPole){
						IsSimilar(testOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getReal(),
								goodOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getReal());
						IsSimilar(testOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getImag(),
								goodOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getImag());
					}

					if(!testOctreeIterator.moveRight()){
						if(goodOctreeIterator.moveRight()){
							uassert(false);
						}
						break;
					}
					if(!goodOctreeIterator.moveRight()){
						uassert(false);
						break;
					}

				} while(true);

				testOctreeIterator.moveUp();
				testOctreeIterator.gotoLeft();

				goodOctreeIterator.moveUp();
				goodOctreeIterator.gotoLeft();
			}
		}
		Print("Over...");
	}


	// set test
	void SetTests(){
		AddTest(&TestSphericalWithPrevious::TestTree,"Test Simu and compare tree");
	}
};



// You must do this
TestClass(TestSphericalWithPrevious)



