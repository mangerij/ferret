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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdexcept>

#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymM2LHandler.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


/**
* In this file we show how to use octree
*/
int main(int argc, char* argv[])
{ 
    FHelpDescribeAndExit(argc, argv, "Generate and store Chebyshev M2L matrices for several orders.");

	// typedefs   
    typedef double FReal;
	typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

	// instantiations
    FTic time;
	MatrixKernelClass MatrixKernel;

	/*
	const unsigned int ORDER = 3;
	const FReal EPSILON = FReal(1e-3);

	const int nnodes = ORDER*ORDER*ORDER;
	ComputeAndCompressAndStoreInBinaryFile<ORDER>(&MatrixKernel, EPSILON);
	
	FReal* K[343];
	int LowRank[343];
	ReadFromBinaryFile<ORDER>(EPSILON, K, LowRank);
	
	for (unsigned int idx=0; idx<343; ++idx) {
		if (K[idx] != NULL) {
			const int rank = LowRank[idx];
			std::cout << "\nM2L " << idx << " has rank " << rank << std::endl;
			for (int i=0; i<nnodes; ++i) {
				for (int r=0; r<rank; ++r)
					std::cout << K[idx][r*nnodes + i] << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
			for (int i=0; i<nnodes; ++i) {
				for (int r=0; r<rank; ++r)
					std::cout << K[idx][rank*nnodes + r*nnodes + i] << " ";
				std::cout << std::endl;
			}
		}
	}

	for (unsigned int idx=0; idx<343; ++idx)
		if (K[idx]!=NULL) delete [] K[idx];
	*/

    ComputeAndCompressAndStoreInBinaryFile<FReal,3>(&MatrixKernel, FReal(1e-3));
    ComputeAndCompressAndStoreInBinaryFile<FReal,5>(&MatrixKernel, FReal(1e-5));
    ComputeAndCompressAndStoreInBinaryFile<FReal,7>(&MatrixKernel, FReal(1e-7));
    ComputeAndCompressAndStoreInBinaryFile<FReal,9>(&MatrixKernel, FReal(1e-9));


	return 0;
}


// [--END--]
