// ===================================================================================
// Copyright ScalFmm 2011 INRIA
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


// system includes
#include <iostream>
#include <stdexcept>


#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FTic.hpp"


#include "../../Src/Kernels/Chebyshev/FChebTensor.hpp"
#include "../../Src/Kernels/Chebyshev/FChebM2LHandler.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/Utils/FParameterNames.hpp"




int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv, "Test Chebyshev precomputation.",
                         FParameterDefinitions::Epsilon);

    // start timer /////////////////////////////////
    FTic time;

    // define set matrix kernel
    typedef double FReal;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    MatrixKernelClass MatrixKernel;

    // constants
    const FReal epsilon     = FParameters::getValue(argc, argv, FParameterDefinitions::Epsilon.options, FReal(0.1));
    const unsigned int order = 9;

    // number of interpolation points per cell
    const unsigned int nnodes = TensorTraits<order>::nnodes;

    // interpolation points of source (Y) and target (X) cell
    FPoint<FReal> X[nnodes], Y[nnodes];
    // set roots of target cell (X)
    FChebTensor<FReal,order>::setRoots(FPoint<FReal>(0.,0.,0.), FReal(2.), X);


    /*
    // allocate memory
    FReal *Qu, *C, *Qb;
    Qu = Qb = C = nullptr;
    unsigned int ninteractions = 0;

    ////////////////////////////////////////////////
    std::cout << "\nAssembly of 316 times "
                        << nnodes << "x" << nnodes << " M2L operators";
    time.tic();
    // compute 316 m2l operators
    ninteractions = 316;
    C = new FReal [nnodes*nnodes * ninteractions];
    unsigned int counter = 0;
    for (int i=-3; i<=3; ++i) {
        for (int j=-3; j<=3; ++j) {
            for (int k=-3; k<=3; ++k) {
                if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
                    // set roots of source cell (Y)
                    const FPoint<FReal> cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
                    FChebTensor<FReal,order>::setRoots(cy, FReal(2.), Y);
                    // evaluate m2l operator
                    for (unsigned int n=0; n<nnodes; ++n)
                        for (unsigned int m=0; m<nnodes; ++m)
                            C[counter*nnodes*nnodes + n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);
                    // increment interaction counter
                    counter++;
                }
            }
        }
    }
    if (counter != 316)
        std::runtime_error("Number of interactions must correspond to 316");
    std::cout << " took " << time.tacAndElapsed() << " sec." << std::endl;
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    std::cout << "\nSVD compression ";
    time.tic();
    const unsigned int rank = Compress<order>(epsilon, ninteractions, Qu, C, Qb);
    std::cout << "to low rank = " << rank << " (eps = " << epsilon
                        << ") took " << time.tacAndElapsed() << " sec." << std::endl;
    ////////////////////////////////////////////////

    // free memory
    if (C  != nullptr) delete [] C;
    if (Qu != nullptr) delete [] Qu;
    if (Qb != nullptr) delete [] Qb;
    */

    ////////////////////////////////////////////////
    // allocate memory
    FReal *Qu1, *C1, *Qb1;
    Qu1 = Qb1 = C1 = nullptr;
    ////////////////////////////////////////////////
    std::cout << "\nAssembly of an " << nnodes << "x" << nnodes << " M2L operator";
    time.tic();
    // compute 316 m2l operators
    C1 = new FReal [nnodes*nnodes];
    const unsigned int i = 3;
    const unsigned int j = 3;
    const unsigned int k = 3;
    const FPoint<FReal> cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
    FChebTensor<FReal,order>::setRoots(cy, FReal(2.), Y);
    // evaluate m2l operator
    for (unsigned int n=0; n<nnodes; ++n)
        for (unsigned int m=0; m<nnodes; ++m)
            C1[n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);
    std::cout << " took " << time.tacAndElapsed() << " sec." << std::endl;
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // get a copy C2 of the M2L operator C1
    FReal *Qu2, *C2, *Qb2;
    Qu2 = Qb2 = C2 = nullptr;
    C2 = new FReal [nnodes * nnodes];
    FBlas::copy(nnodes*nnodes, C1, C2);
    // Omega_x^{1/2} C2 Omega_y^{1/2}
    FReal weights[nnodes];
    FChebTensor<FReal,order>::setRootOfWeights(weights);
    for (unsigned int n=0; n<nnodes; ++n) {
        FBlas::scal(nnodes, weights[n], C2+n, nnodes); // scale rows
        FBlas::scal(nnodes, weights[n], C2+n*nnodes);  // scale cols
    }
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    std::cout << "\nSVD compression of K ";
    time.tic();
    const unsigned int rank1 = Compress<FReal,order>(epsilon, 1, Qu1, C1, Qb1);
    std::cout << "to low rank = " << rank1 << " (eps = " << epsilon
              << ") took " << time.tacAndElapsed() << " sec." << std::endl;
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    std::cout << "SVD compression of Omega_x^{1/2} K Omega_y^{1/2} ";
    time.tic();
    const unsigned int rank2 = Compress<FReal,order>(epsilon, 1, Qu2, C2, Qb2);
    std::cout << "to low rank = " << rank2 << " (eps = " << epsilon
              << ") took " << time.tacAndElapsed() << " sec." << std::endl;
    ////////////////////////////////////////////////

    // free memory
    if (C1  != nullptr) delete [] C1;
    if (Qu1 != nullptr) delete [] Qu1;
    if (Qb1 != nullptr) delete [] Qb1;

    if (C2  != nullptr) delete [] C2;
    if (Qu2 != nullptr) delete [] Qu2;
    if (Qb2 != nullptr) delete [] Qb2;


    return 0;
} 	
