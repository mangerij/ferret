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
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
 * Test functionality of C - interfaced LAPACK functions
 */

int main(int argc, char ** argv)
{
    FHelpDescribeAndExit(argc, argv, "Test the lapack compilation and linking (only the code is interesting).");

    /*
   * List of tested functions:
   * Cholesky decomposition: FBlas::potrf()
   * TODO SVD: FBlas::gesvd()
   * TODO QR decomposition: FBlas::geqrf()
   */

    typedef double FReal;
    const unsigned int m = 4, n = 4;
    FReal* A = new FReal [m * n]; // matrix: column major ordering

    // A= LL^T ////////////////////////////////////
    // define symmetric definite positive matrix A
    A[0]=5; A[10]=4; A[15]=7;
    A[1]=A[3]=A[4]=A[12]=2;
    A[6]=A[7]=A[9]=A[13]=1;
    A[2]=A[5]=A[8]=3;
    A[11]=A[14]=-1;

    // copy A in C
    FReal* C = new FReal [m * n]; // matrix: column major ordering
    for (unsigned int ii=0; ii<m; ++ii)
        for (unsigned int jj=0; jj<n; ++jj)
            C[ii*m + jj]=A[ii*m + jj];

    std::cout<<"\nA=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << A[i*n+j] << " ";
        std::cout<< std::endl;
    }
    std::cout<<"]"<<std::endl;

    // perform Cholesky decomposition
    std::cout<<"\nCholesky decomposition ";
    int INF = FBlas::potrf(m, A, n);
    if(INF==0) {std::cout<<"succeeded!"<<std::endl;}
    else {std::cout<<"failed!"<<std::endl;}

    std::cout<<"\nA_out=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << A[i*n+j] << " ";
        std::cout<<std::endl;
    }
    std::cout<<"]"<<std::endl;

    // build lower matrix
    FReal* L = new FReal [m * n]; // matrix: column major ordering
    for (unsigned int ii=0; ii<m; ++ii)
        for (unsigned int jj=0; jj<n; ++jj){
            if(ii<=jj)
                L[ii*m + jj]=A[ii*m + jj];
            else
                L[ii*m + jj]=0.;
        }

    std::cout<<"\nL=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << L[i*n+j] << " ";
        std::cout<< std::endl;
    }
    std::cout<<"]"<<std::endl;

    // verify result by computing B=LL^T
    FReal* B = new FReal [m * n]; // matrix: column major ordering
    for (unsigned int ii=0; ii<m; ++ii)
        for (unsigned int jj=0; jj<n; ++jj){
            B[ii*m + jj]=0.;
            for (unsigned int j=0; j<n; ++j)
                B[ii*m + jj]+=L[j*m + ii]*L[j*m + jj];
        }

    std::cout<<"\nA-LL^T=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << B[i*n+j]-C[i*n+j] << " ";
        std::cout<< std::endl;
    }
    std::cout<<"]"<<std::endl;

    delete [] A;
    delete [] B;
    delete [] C;

    return 0;
}
