// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"


#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
* In this file we show how to use octree
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Test the octree with the Chebyshev kernel.");

	const int ORDER = 5;

    typedef double FReal;
    typedef FP2PParticleContainer<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
	typedef FChebCell<FReal,ORDER> CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
	
	///////////////////////What we do/////////////////////////////
	std::cout << ">> This executable is useless to execute.\n";
	std::cout << ">> It is only interesting to understand the code\n";
	std::cout << ">> and how to use the Octree\n";
	//////////////////////////////////////////////////////////////
	
    const long NbPart = 100000;
	FTic counter;
	
    srand48( static_cast<unsigned int>(time(NULL)) );

	const FReal BoxWidth = 1.;
    const FPoint<FReal> BoxCenter(.5, .5, .5);
	const unsigned int TreeHeight = 10;
	OctreeClass tree(TreeHeight, 3, BoxWidth, BoxCenter);

	// -----------------------------------------------------
	std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
	counter.tic();
    {
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            tree.insert(FPoint<FReal>(FReal(drand48()),FReal(drand48()),FReal(drand48())));
		}
	}
	counter.tac();
	std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
	// -----------------------------------------------------


	
	// Check if particles are strictly within its containing leaf cells
	{
		const FReal BoxWidthLeaf = BoxWidth / FReal(FMath::pow(2, TreeHeight-1));
        tree.forEachCellLeaf([&](CellClass* LeafCell, LeafClass* leaf){
            const FPoint<FReal> Origin(BoxCenter - BoxWidth / FReal(2.));
            const FPoint<FReal> LeafCellCenter(Origin.getX() + (FReal(LeafCell->getCoordinate().getX()) + FReal(.5)) * BoxWidthLeaf,
																			 Origin.getY() + (FReal(LeafCell->getCoordinate().getY()) + FReal(.5)) * BoxWidthLeaf,
																			 Origin.getZ() + (FReal(LeafCell->getCoordinate().getZ()) + FReal(.5)) * BoxWidthLeaf);

            const ContainerClass *const Particles = leaf->getSrc();
            const FReal*const positionsX = Particles->getPositions()[0];
            const FReal*const positionsY = Particles->getPositions()[1];
            const FReal*const positionsZ = Particles->getPositions()[2];

            for(FSize idxPart = 0 ; idxPart < Particles->getNbParticles() ; ++idxPart){
                const FPoint<FReal> distance(LeafCellCenter-FPoint<FReal>(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]));
				if (std::abs(distance.getX())>BoxWidthLeaf/FReal(2.) ||
						std::abs(distance.getY())>BoxWidthLeaf/FReal(2.) ||
						std::abs(distance.getZ())>BoxWidthLeaf/FReal(2.)) {
					std::cout << "Particle (center - particle = " << distance << " < " << BoxWidthLeaf/FReal(2.) << ") is out of cell. STOP"
										<< std::endl;
					//exit(-1);
                }
			}
        });
	}

	return 0;
}


// [--END--]
