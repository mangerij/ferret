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
#include <stdlib.h>

#include "../../Src/Utils/FBlas.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


template <class FReal>
FReal FRandom() { return (FReal(drand48())); }

/**
 * Test functionality of C - interfaced BLAS functions
 */

int main(int argc, char** argv)
{
    FHelpDescribeAndExit(argc, argv, "Simply ensure that blas are compuling and running.");

    typedef double FReal;
	const unsigned int m = 4, n = 4; // to be able to test both, transpose and not transpose operations
	FReal* A = new FReal [m * n]; // matrix: column major ordering
	FReal* x = new FReal [n];
	FReal* y = new FReal [m];
	
	const FReal d = FRandom<FReal>();
	
	for (unsigned int j=0; j<n; ++j) {
		x[j] = FRandom<FReal>();
		for (unsigned int i=0; i<m; ++i) {
			A[j*m + i] = FRandom<FReal>();
		}
	}

//	std::cout << "A = " << std::endl;
//	for (unsigned int i=0; i<m; ++i) {
//		for (unsigned int j=0; j<n; ++j) 
//			std::cout << A[j*m + i] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;	
//	
//	std::cout << "x = " << std::endl;
//	for (unsigned int i=0; i<m; ++i)
//		std::cout << x[i] << std::endl;
//	std::cout << std::endl;	


	
	FReal* z = new FReal [m];
	

	// y = d Ax ////////////////////////////////////
	// cblas
	FBlas::gemv(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		z[i] = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			z[i] += A[j*m + i] * x[j];
		z[i] *= d;
	}
	// compare
	std::cout << "\ny = d Ax (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;


	// y = d A^Tx ////////////////////////////////////
	// cblas
	FBlas::gemtv(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		z[i] = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			z[i] += A[i*m + j] * x[j];
		z[i] *= d;
	}
	// compare
	std::cout << "\ny = d A^Tx (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;


	// y += d Ax ////////////////////////////////////
	// cblas
	FBlas::gemva(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		FReal _z = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			_z += A[j*m + i] * x[j];
		z[i] += _z * d;
	}
	// compare
	std::cout << "\ny += d Ax (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;


	// y += d A^Tx ////////////////////////////////////
	// cblas
	FBlas::gemtva(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		FReal _z = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			_z += A[i*m + j] * x[j];
		z[i] += _z * d;
	}
	// compare
	std::cout << "\ny += d A^Tx (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;



	delete [] A;
	delete [] x;
	delete [] y;
	delete [] z;

	
	return 0;
}

