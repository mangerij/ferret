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

#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"


#include "../../Src/Kernels/Chebyshev/FChebInterpolator.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


/**
* In this file we show how to use octree
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Test Chebyshev interpolator.");

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
	FReal width = FReal(3.723);

	////////////////////////////////////////////////////////////////////
	LeafClass X;
    FPoint<FReal> cx(0., 0., 0.);
	const unsigned long M = 20000;
	std::cout << "Fill the leaf X of width " << width
						<< " centered at cx=" << cx << " with M=" << M << " target particles" << std::endl;
    {
		for(unsigned long i=0; i<M; ++i){
            FReal x = (FReal(drand48()) - FReal(.5)) * width + cx.getX();
            FReal y = (FReal(drand48()) - FReal(.5)) * width + cx.getY();
            FReal z = (FReal(drand48()) - FReal(.5)) * width + cx.getZ();
            X.push(FPoint<FReal>(x, y, z), FReal(drand48()));
		}
	}


	////////////////////////////////////////////////////////////////////
	LeafClass Y;
    FPoint<FReal> cy(FReal(2.)*width, 0., 0.);
	const unsigned long N = 20000;
	std::cout << "Fill the leaf Y of width " << width
						<< " centered at cy=" << cy	<< " with N=" << N << " target particles" << std::endl;
    {
		for(unsigned long i=0; i<N; ++i){
            FReal x = (FReal(drand48()) - FReal(.5)) * width + cy.getX();
            FReal y = (FReal(drand48()) - FReal(.5)) * width + cy.getY();
            FReal z = (FReal(drand48()) - FReal(.5)) * width + cy.getZ();
            Y.push(FPoint<FReal>(x, y, z), FReal(drand48()));
		}
	}



	////////////////////////////////////////////////////////////////////
	// approximative computation
	const unsigned int ORDER = 10;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
    typedef FChebInterpolator<FReal,ORDER,MatrixKernelClass> InterpolatorClass;
	InterpolatorClass S;

	std::cout << "\nCompute interactions approximatively, interpolation order = " << ORDER << " ..." << std::endl;

	std::cout << "\nP2M ... " << std::flush;
	time.tic();
	// Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
	FReal W[nnodes]; // multipole expansion
	S.applyP2M(cy, width, W, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M
	std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

	std::cout << "M2L ... " << std::flush;
	time.tic();
	// Multipole to local: F_m = \sum_n^L K(\bar x_m, \bar y_n) * W_n
    FPoint<FReal> rootsX[nnodes], rootsY[nnodes];
	FChebTensor<FReal,ORDER>::setRoots(cx, width, rootsX);
	FChebTensor<FReal,ORDER>::setRoots(cy, width, rootsY);

	FReal F[nnodes]; // local expansion
	for (unsigned int i=0; i<nnodes; ++i) {
		F[i] = FReal(0.);
		for (unsigned int j=0; j<nnodes; ++j)
			F[i] += MatrixKernel.evaluate(rootsX[i], rootsY[j]) * W[j];
	}
	std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

	std::cout << "L2P (potential) ... " << std::flush;
	time.tic();
	// Interpolate p_i = \sum_m^L S(x_i,\bar x_m) * F_m
	S.applyL2P(cx, width, F, X.getTargets());
	std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

	std::cout << "L2P (forces) ... " << std::flush;
	time.tic();
	// Interpolate f_i = \sum_m^L P(x_i,\bar x_m) * F_m
	S.applyL2PGradient(cx, width, F, X.getTargets());
	std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

	
	////////////////////////////////////////////////////////////////////
	// direct computation
	std::cout << "Compute interactions directly ..." << std::endl;
	time.tic();

	FReal* approx_f = new FReal [M * 3];
	FReal*        f = new FReal [M * 3];
	FBlas::setzero(M*3, f);

	FReal* approx_p = new FReal[M];
	FReal*        p = new FReal[M];
	FBlas::setzero(M, p);

	{ // start direct computation
		unsigned int counter = 0;
		
        for(FSize idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
            const FPoint<FReal> x = FPoint<FReal>(X.getSrc()->getPositions()[0][idxPartX],
                                     X.getSrc()->getPositions()[1][idxPartX],
                                     X.getSrc()->getPositions()[2][idxPartX]);
            const FReal  wx = X.getSrc()->getPhysicalValues()[idxPartX];
			
            for(FSize idxPartY = 0 ; idxPartY < Y.getSrc()->getNbParticles() ; ++idxPartY){
                const FPoint<FReal> y = FPoint<FReal>(Y.getSrc()->getPositions()[0][idxPartY],
                                         Y.getSrc()->getPositions()[1][idxPartY],
                                         Y.getSrc()->getPositions()[2][idxPartY]);
                const FReal  wy = Y.getSrc()->getPhysicalValues()[idxPartY];

				const FReal one_over_r = MatrixKernel.evaluate(x, y);
				// potential
				p[counter] += one_over_r * wy;
				// force
                FPoint<FReal> force(y - x);
				force *= one_over_r*one_over_r*one_over_r;
				f[counter*3 + 0] += force.getX() * wx * wy;
				f[counter*3 + 1] += force.getY() * wx * wy;
				f[counter*3 + 2] += force.getZ() * wx * wy;
			}
			
            counter++;
		}
	} // end direct computation


	time.tac();
	std::cout << "Done in " << time.elapsed() << "sec." << std::endl;


    ////////////////////////////////////////////////////////////////////
	unsigned int counter = 0;
    for(FSize idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
        approx_p[counter] = X.getSrc()->getPotentials()[idxPartX];
        const FPoint<FReal> force = FPoint<FReal>(X.getSrc()->getForcesX()[idxPartX],
                                    X.getSrc()->getForcesY()[idxPartX],
                                    X.getSrc()->getForcesZ()[idxPartX]);
		approx_f[counter*3 + 0] = force.getX();
		approx_f[counter*3 + 1] = force.getY();
		approx_f[counter*3 + 2] = force.getZ();

        counter++;
	}

	std::cout << "\nPotential error:" << std::endl;
    std::cout << "Relative error   = " << FMath::FAccurater<FReal>( p, approx_p, M) << std::endl;

	std::cout << "\nForce error:" << std::endl;
    std::cout << "Relative L2 error   = " << FMath::FAccurater<FReal>( f, approx_f, M*3) << std::endl;
	std::cout << std::endl;

	// free memory
	delete [] approx_p;
	delete [] p;
	delete [] approx_f;
	delete [] f;
	

	return 0;
}



