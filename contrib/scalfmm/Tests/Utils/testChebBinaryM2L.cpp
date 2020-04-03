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

#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebM2LHandler.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


/**
* In this file we show how to use octree
*/
int main(int argc, char* argv[])
{ 
    FHelpDescribeAndExit(argc, argv, "Just generate Chebyshev M2L matrices for several orders.");

	// typedefs   
    typedef double FReal;
	typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
	

	// instantiations
	FTic time;
	MatrixKernelClass MatrixKernel;

	/*
	// constants
  const FReal epsilon = FReal(atof(argv[1]));
	const unsigned int order = 4;

	// write precomputed compressed M2l operators to binary file 
	std::cout << "\nCompute compressed M2L operators of ACC("
						<< order << ", " << epsilon << ") and write them to a binary file ..."
						<< std::endl;
	time.tic();
    typedef FChebM2LHandler<FReal,order,MatrixKernel> M2LHandlerClass;
	M2LHandlerClass::ComputeAndCompressAndStoreInBinaryFile(epsilon);
	time.tac();
	std::cout << " in " << time.elapsed() << "sec." << std::endl;

	// read precomputed compressed M2l operators from binary file 
	std::cout << "\nRead compressed M2L operators of ACC(" << order << ", " << epsilon << ") from a binary file ..." << std::endl;
	time.tic();
	M2LHandlerClass M2L(epsilon);
	M2L.ReadFromBinaryFileAndSet();
	time.tac();
	std::cout << " in " << time.elapsed() << "sec." << std::endl;
	*/

	// order 2
	time.tic();
    FChebM2LHandler<FReal,2,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-1));
    FChebM2LHandler<FReal,2,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-2));
	// order 3
    FChebM2LHandler<FReal,3,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-2));
    FChebM2LHandler<FReal,3,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-3));
	// order 4
    FChebM2LHandler<FReal,4,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-3));
    FChebM2LHandler<FReal,4,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-4));
	// order 5
    FChebM2LHandler<FReal,5,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-4));
    FChebM2LHandler<FReal,5,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-5));
	// order 6
    FChebM2LHandler<FReal,6,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-5));
    FChebM2LHandler<FReal,6,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-6));
	// order 7
    FChebM2LHandler<FReal,7,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-6));
    FChebM2LHandler<FReal,7,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-7));
	// order 8
    FChebM2LHandler<FReal,8,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-7));
    FChebM2LHandler<FReal,8,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-8));
	// order 9
    FChebM2LHandler<FReal,9,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-8));
    FChebM2LHandler<FReal,9,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-9));
	// order 10
    FChebM2LHandler<FReal,10,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-9));
    FChebM2LHandler<FReal,10,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(&MatrixKernel,FReal(1e-10));


	return 0;
}


// [--END--]
