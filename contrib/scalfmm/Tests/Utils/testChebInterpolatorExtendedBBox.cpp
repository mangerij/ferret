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

#include "Utils/FTic.hpp"
#include "Utils/FMath.hpp"

#include "Containers/FVector.hpp"

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"


#include "Kernels/Chebyshev/FChebInterpolator.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"
#include "Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


/**
 * Here we compute the interactions between particles contained 
 * in 2 given cluster using the Chebyshev interpolation scheme.
 * This test also illustrates the effect of inserting particles
 * in the extended bounding box of each clusters. 
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
	std::cout << " direct computation.\n";
  std::cout << " Interpolation is performed on the extended bounding boxes since particles\n"; 
  std::cout << " are duplicated and might be located outside the leafs.\n\n";
	//////////////////////////////////////////////////////////////

  MatrixKernelClass MatrixKernel;
	FTic time;

  /// Duplication rule
  // distance between original particule and duplicata
  FReal dradius = .25; // segment size / 2.
    
  // Leaf size 
  FReal alpha = FReal(2.); // width / dradius
  FReal width = FReal(alpha*(2.*dradius)); // alpha * segment size

  // relative position of clusters (X centered at origin)
  // 2/alpha < beta < 2: well separated
  // 1 < beta <= 2/alpha: extended bboxes overlap
  // beta <=1: leaves overlap
  FReal beta = FReal(2.); // cy[0]/width

  // BBox extension
  FReal ExtendedCellWidth = width + FReal(2.)*dradius;

	////////////////////////////////////////////////////////////////////
	LeafClass X;
    FPoint<FReal> cx(0., 0., 0.);
	const unsigned long M = 20000;
	std::cout << "Fill the leaf X of width " << width
						<< " centered at cx=" << cx 
            << " with M=" << M << " target particles";
  int countPartOutX(0);
  {
		for(unsigned long i=0; i<M/2; ++i){
      // generate and insert original particle
      FReal x = (FReal(drand48()) - FReal(.5)) * width + cx.getX();
      FReal y = (FReal(drand48()) - FReal(.5)) * width + cx.getY();
      FReal z = (FReal(drand48()) - FReal(.5)) * width + cx.getZ();
      X.push(FPoint<FReal>(x, y, z), FReal(drand48()));
      // generate and insert duplicata (in 1 random direction)
      FReal ddir[3] = {FReal(drand48()),FReal(drand48()),FReal(drand48())};
      FReal norm_ddir = sqrt(ddir[0]*ddir[0]+ddir[1]*ddir[1]+ddir[2]*ddir[2]);
      FReal xd = x + dradius*ddir[0]/norm_ddir;
      FReal yd = y + dradius*ddir[1]/norm_ddir;
      FReal zd = z + dradius*ddir[2]/norm_ddir;
      X.push(FPoint<FReal>(xd, yd, zd), FReal(drand48()));
      // count duplicata that are outside original leaf
      bool isinX = ((xd<cx.getX() + 0.5*width) && (xd>cx.getX() - 0.5*width));
      bool isinY = ((yd<cx.getY() + 0.5*width) && (yd>cx.getY() - 0.5*width));
      bool isinZ = ((zd<cx.getZ() + 0.5*width) && (zd>cx.getZ() - 0.5*width));
      if(!isinX && !isinY && !isinZ) countPartOutX++;




		}
	}
  std::cout << " ( " << countPartOutX << "/"<< M/2 << " are out)." << std::endl;

	////////////////////////////////////////////////////////////////////
	LeafClass Y;
    FPoint<FReal> cy(beta*width, 0., 0.);
	const unsigned long N = 20000;
	std::cout << "Fill the leaf Y of width " << width
						<< " centered at cy=" << cy	
            << " with N=" << N << " target particles";
  int countPartOutY(0);
  {
		for(unsigned long i=0; i<N/2; ++i){
      // generate and insert original particle
      FReal x = (FReal(drand48()) - FReal(.5)) * width + cy.getX();
      FReal y = (FReal(drand48()) - FReal(.5)) * width + cy.getY();
      FReal z = (FReal(drand48()) - FReal(.5)) * width + cy.getZ();
      Y.push(FPoint<FReal>(x, y, z), FReal(drand48()));
      // generate and insert duplicata (in 1 random direction)
      FReal ddir[3] = {FReal(drand48()),FReal(drand48()),FReal(drand48())};
      FReal norm_ddir = sqrt(ddir[0]*ddir[0]+ddir[1]*ddir[1]+ddir[2]*ddir[2]);
      FReal xd = x + dradius*ddir[0]/norm_ddir;
      FReal yd = y + dradius*ddir[1]/norm_ddir;
      FReal zd = z + dradius*ddir[2]/norm_ddir;
      Y.push(FPoint<FReal>(xd, yd, zd), FReal(drand48()));
      // count duplicata that are outside original leaf
      bool isinX = ((xd<cy.getX() + 0.5*width) && (xd>cy.getX() - 0.5*width));
      bool isinY = ((yd<cy.getY() + 0.5*width) && (yd>cy.getY() - 0.5*width));
      bool isinZ = ((zd<cy.getZ() + 0.5*width) && (zd>cy.getZ() - 0.5*width));
      if(!isinX && !isinY && !isinZ) countPartOutY++;

		}
	}
  std::cout << " ( " << countPartOutY << "/"<< N/2 << " are out)." << std::endl;



	////////////////////////////////////////////////////////////////////
	// approximative computation
	const unsigned int ORDER = 10;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
	typedef FChebInterpolator<FReal,ORDER,MatrixKernelClass> InterpolatorClass;
	InterpolatorClass S; // default ctor is used since no M2M/L2L op required

	std::cout << "\nCompute interactions approximatively, interpolation order = " << ORDER << " ..." << std::endl;

	std::cout << "\nP2M ... " << std::flush;
	time.tic();
	// Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
	FReal W[nnodes]; // multipole expansion
	S.applyP2M(cy, ExtendedCellWidth, W, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M
	std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

	std::cout << "M2L ... " << std::flush;
	time.tic();
	// Multipole to local: F_m = \sum_n^L K(\bar x_m, \bar y_n) * W_n
    FPoint<FReal> rootsX[nnodes], rootsY[nnodes];
	FChebTensor<FReal,ORDER>::setRoots(cx, ExtendedCellWidth, rootsX);
	FChebTensor<FReal,ORDER>::setRoots(cy, ExtendedCellWidth, rootsY);

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
	S.applyL2P(cx, ExtendedCellWidth, F, X.getTargets());
	std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

	std::cout << "L2P (forces) ... " << std::flush;
	time.tic();
	// Interpolate f_i = \sum_m^L P(x_i,\bar x_m) * F_m
	S.applyL2PGradient(cx, ExtendedCellWidth, F, X.getTargets());
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



