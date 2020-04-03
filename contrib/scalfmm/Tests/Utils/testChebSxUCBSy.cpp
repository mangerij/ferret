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

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"


#include "../../Src/Kernels/Chebyshev/FChebInterpolator.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


template <class FReal>
void applyM2M(FReal *const S,	FReal *const w, const unsigned int n,	FReal *const W, const unsigned int N)
{ FBlas::gemtva(n, N, FReal(1.), S,	w, W); }


/**
* In this file we show how to use octree
*/

int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv, "Test Chebyshev interation computations.");

    typedef double FReal;
    typedef FP2PParticleContainer<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
	typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

	///////////////////////What we do/////////////////////////////
	std::cout << "\nTask: Compute interactions between source particles in leaf Y and target\n";
	std::cout << " particles in leaf X. Compare the fast summation K ~ Sx K Sy' with the\n";
	std::cout << " direct computation.\n" << std::endl;
	//////////////////////////////////////////////////////////////

    MatrixKernelClass MatrixKernel;
	FTic time;

	
	// Leaf size
    FReal width(FReal(drand48()) * FReal(10.));

	////////////////////////////////////////////////////////////////////
	LeafClass X;
    FPoint<FReal> cx(0., 0., 0.);
	const long M = 1000;
	std::cout << "Fill the leaf X of width " << width
						<< " centered at cx=[" << cx.getX() << "," << cx.getY() << "," << cx.getZ()
						<< "] with M=" << M << " target particles" << std::endl;
    {
		for(long i=0; i<M; ++i){
            FReal x = (FReal(drand48()) - FReal(.5)) * width + cx.getX();
            FReal y = (FReal(drand48()) - FReal(.5)) * width + cx.getY();
            FReal z = (FReal(drand48()) - FReal(.5)) * width + cx.getZ();
            X.push(FPoint<FReal>(x,y,z));
		}
	}


	////////////////////////////////////////////////////////////////////
	LeafClass Y;
    FPoint<FReal> cy(FReal(2.)*width, 0., 0.);
	const long N = 1000;
	std::cout << "Fill the leaf Y of width " << width
						<< " centered at cy=[" << cy.getX() << "," << cy.getY() << "," << cy.getZ()
						<< "] with N=" << N << " target particles" << std::endl;
    {
		for(long i=0; i<N; ++i){
            FReal x = (FReal(drand48()) - FReal(.5)) * width + cy.getX();
            FReal y = (FReal(drand48()) - FReal(.5)) * width + cy.getY();
            FReal z = (FReal(drand48()) - FReal(.5)) * width + cy.getZ();
            Y.push(FPoint<FReal>(x, y, z),FReal(drand48()));
		}
	}


	////////////////////////////////////////////////////////////////////
	// approximative computation
	const unsigned int ORDER = 10;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
	typedef FChebInterpolator<FReal,ORDER,MatrixKernelClass> InterpolatorClass;
	InterpolatorClass S;
	

	std::cout << "\nCompute interactions approximatively, ORDER = " << ORDER << std::endl;
	FReal W[nnodes]; // multipole expansion
	FReal F[nnodes]; // local expansion
	for (unsigned int n=0; n<nnodes; ++n)
		W[n] = F[n] = FReal(0.);

	// Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
	time.tic();
	S.applyP2M(cy, width, W, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M
	std::cout << "P2M done in " << time.tacAndElapsed() << "s" << std::endl;

	// M2L (direct)
    FPoint<FReal> rootsX[nnodes], rootsY[nnodes];
	FChebTensor<FReal,ORDER>::setRoots(cx, width, rootsX);
	FChebTensor<FReal,ORDER>::setRoots(cy, width, rootsY);

	{
		for (unsigned int i=0; i<nnodes; ++i) {
			F[i] = FReal(0.);
			for (unsigned int j=0; j<nnodes; ++j)
				F[i] += MatrixKernel.evaluate(rootsX[i], rootsY[j]) * W[j];
		}
	}

//	{
//		for (unsigned int ix=0; ix<ORDER; ++ix)
//			for (unsigned int jx=0; jx<ORDER; ++jx)
//				for (unsigned int kx=0; kx<ORDER; ++kx)  {
//					const unsigned int idx = kx*ORDER*ORDER + jx*ORDER + ix;
//					F[idx] = FReal(0.);
//					for (unsigned int iy=0; iy<ORDER; ++iy)
//						for (unsigned int jy=0; jy<ORDER; ++jy)
//							for (unsigned int ky=0; ky<ORDER; ++ky) {
//								const unsigned int idy = ky*ORDER*ORDER + jy*ORDER + iy;
//								F[idx] += MatrixKernel.evaluate(rootsX[idx], rootsY[idy]) * W[idy];
//							}
//				}
//	}
		

	// Interpolate f_i = \sum_m^L S(x_i,\bar x_m) * F_m
	time.tic();
	//S.applyL2PTotal(cx, width, F, X.getTargets());
	S.applyL2PTotal(cx, width, F, X.getTargets());
	std::cout << "L2P done in " << time.tacAndElapsed() << "s" << std::endl;

	// -----------------------------------------------------

	////////////////////////////////////////////////////////////////////
	// direct computation
	std::cout << "Compute interactions directly ..." << std::endl;
	time.tic();

	FReal* approx_f = new FReal[M];
	FReal*        f = new FReal[M];
    {
        for (unsigned int i=0; i<M; ++i) f[i] = FReal(0.);

        const FReal*const positionsX = Y.getSrc()->getPositions()[0];
        const FReal*const positionsY = Y.getSrc()->getPositions()[1];
        const FReal*const positionsZ = Y.getSrc()->getPositions()[2];
        const FReal*const physicalValues = Y.getSrc()->getPhysicalValues();

        for(FSize idxPart = 0 ; idxPart < Y.getSrc()->getNbParticles() ; ++idxPart){
            const FPoint<FReal> y = FPoint<FReal>(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]);
            const FReal        w = physicalValues[idxPart];

            const FReal*const xpositionsX = X.getSrc()->getPositions()[0];
            const FReal*const xpositionsY = X.getSrc()->getPositions()[1];
            const FReal*const xpositionsZ = X.getSrc()->getPositions()[2];

            for(FSize idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
                const FPoint<FReal> x = FPoint<FReal>(xpositionsX[idxPart],xpositionsY[idxPart],xpositionsZ[idxPart]);
                f[idxPartX] += MatrixKernel.evaluate(x,y) * w;
            }
        }
    }
	time.tac();
	std::cout << "Done in " << time.elapsed() << "sec." << std::endl;


	////////////////////////////////////////////////////////////////////
    const FReal*const potentials = X.getSrc()->getPotentials();
    for(FSize idxPart = 0 ; idxPart < X.getSrc()->getNbParticles() ; ++idxPart){
        approx_f[idxPart] = potentials[idxPart];
	}

	
	//for (unsigned int i=0; i<8; ++i)
	//	std::cout << f[i] << "\t" << approx_f[i] << "\t" << approx_f[i]/f[i] << std::endl;
	

    std::cout << "\nRelative L2 error  = " << FMath::FAccurater<FReal>( f, approx_f, M) << std::endl;
    std::cout << "Relative Lmax error = "  << FMath::FAccurater<FReal>( f, approx_f, M) << "\n" << std::endl;

	// free memory
	delete [] approx_f;
	delete [] f;


	return 0;
}


// [--END--]
