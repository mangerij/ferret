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


// system includes
#include <iostream>
#include <stdexcept>
#include <climits>



#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FBlas.hpp"

#include "../../Src/Kernels/Chebyshev/FChebTensor.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Interpolation/FInterpSymmetries.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

template < class FReal, int ORDER>
void permuteMatrix(const unsigned int perm[ORDER*ORDER*ORDER], FReal *const Matrix)
{
    const unsigned int nnodes = ORDER*ORDER*ORDER;
    // allocate temporary memory for matrix
    FReal temp[nnodes*nnodes];
    // permute rows
    for (unsigned int r=0; r<nnodes; ++r)
        FBlas::copy(nnodes, Matrix+r, nnodes, temp+perm[r], nnodes);
    // permute columns
    for (unsigned int c=0; c<nnodes; ++c)
        FBlas::copy(nnodes, temp + c * nnodes, Matrix + perm[c] * nnodes);
}





int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv, "Test Chebyshev symmetries.");

    // start timer /////////////////////////////////
    FTic time;

    // define set matrix kernel
    typedef double FReal;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    MatrixKernelClass MatrixKernel;

    // constants
    //const FReal epsilon     = FReal(atof(argv[1]));
    const unsigned int order = 5;

    // number of interpolation points per cell
    const unsigned int nnodes = TensorTraits<order>::nnodes;

    // interpolation points of source (Y) and target (X) cell
    FPoint<FReal> X[nnodes], Y[nnodes];
    // set roots of target cell (X)
    FChebTensor<FReal,order>::setRoots(FPoint<FReal>(0.,0.,0.), FReal(2.), X);

    // allocate 343 pointers to K, but only 16 are actually filled
    FReal** K = new FReal* [343];
    for (unsigned int t=0; t<343; ++t) K[t] = nullptr;

    {
        unsigned int counter = 0;

        for (int i=2; i<=3; ++i) {
            for (int j=0; j<=i; ++j) {
                for (int k=0; k<=j; ++k) {

                    const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
                    K[idx] = new FReal [nnodes*nnodes];

                    //std::cout << i << "," << j << "," << k << "\t" << idx << std::endl;

                    const FPoint<FReal> cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
                    FChebTensor<FReal,order>::setRoots(cy, FReal(2.), Y);
                    for (unsigned int n=0; n<nnodes; ++n)
                        for (unsigned int m=0; m<nnodes; ++m)
                            K[idx][n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);

                    counter++;
                }
            }
        }

        std::cout << "num interactions = " << counter << std::endl;
    }

    // max difference
    FReal maxdiff(0.);

    // permuter
    FInterpSymmetries<order> permuter;

    // permutation vector
    unsigned int perm[nnodes];

    FReal* K0 = new FReal [nnodes*nnodes];

    unsigned int counter = 0;

    for (int i=-3; i<=3; ++i) {
        for (int j=-3; j<=3; ++j) {
            for (int k=-3; k<=3; ++k) {
                if (abs(i)>1 || abs(j)>1 || abs(k)>1) {

                    const FPoint<FReal> cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
                    FChebTensor<FReal,order>::setRoots(cy, FReal(2.), Y);
                    for (unsigned int n=0; n<nnodes; ++n)
                        for (unsigned int m=0; m<nnodes; ++m)
                            K0[n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);

                    // permute
                    const unsigned int pidx = permuter.getPermutationArrayAndIndex(i, j, k, perm);
                    permuteMatrix<FReal,order>(perm, K0);

                    if (K[pidx]==NULL) std::cout << " - not existing index " << pidx << std::endl;

                    FReal mdiff(0.);
                    for (unsigned int n=0; n<nnodes*nnodes; ++n) {
                        FReal diff(K0[n] - K[pidx][n]);
                        if (FMath::Abs(diff)>mdiff) mdiff = diff;
                    }
                    if (FMath::Abs(mdiff)>maxdiff) maxdiff = FMath::Abs(mdiff);

                    if (mdiff > 1e-15) exit(-1);
                    counter++;

                }
            }
        }
    }

    std::cout << "Max error = " << maxdiff << " of counter = " << counter << std::endl;

    for (unsigned int t=0; t<343; ++t) if (K[t]!=NULL) delete [] K[t];
    delete [] K;
    delete [] K0;

    return 0;
} 	



