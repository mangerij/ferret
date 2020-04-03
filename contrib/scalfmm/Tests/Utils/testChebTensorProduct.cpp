// ===================================================================================
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

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FBlas.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Chebyshev/FChebTensor.hpp"
#include "../../Src/Kernels/Chebyshev/FChebInterpolator.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

template <class FReal>
void applyM2M(FReal *const S,	FReal *const w, const unsigned int n,	FReal *const W, const unsigned int N)
{ FBlas::gemtva(n, N, FReal(1.), S,	w, W); }

template <class FReal>
void applym2m(FReal *const S,	FReal *const w, const unsigned int n,	FReal *const W, const unsigned int N)
{ FBlas::gemtm(n, n, n*n, FReal(1.), S, n, w, n, W, n); }

template <class FReal>
void applyL2L(FReal *const S,	FReal *const F, const unsigned int n,	FReal *const f, const unsigned int N)
{ FBlas::gemva(N, n, FReal(1.), S, F, f);	}

template <class FReal>
void applyl2l(FReal *const S,	FReal *const F, const unsigned int n,	FReal *const f, const unsigned int N)
{ FBlas::gemm(n, n, n*n, FReal(1.), S, n, F, n, f, n); }


/**
 * In this file we show how to use octree
 */

int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv, "Test Chebyshev tensor product.");

    typedef double FReal;
	const unsigned int ORDER = 10;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
    FPoint<FReal> X[nnodes];
	FChebInterpolator<FReal,ORDER> Interpolator;

	{
        FChebTensor<FReal,ORDER>::setRoots(FPoint<FReal>(0.,0.,0.), FReal(2.), X);
        FPoint<FReal> x[nnodes];
	
        //	const FPoint<FReal> cx(-.5, -.5, -.5);
        //	const FPoint<FReal> cx(-.5, -.5,  .5);
        //	const FPoint<FReal> cx(-.5,  .5, -.5);
        //	const FPoint<FReal> cx(-.5,  .5,  .5);
        //	const FPoint<FReal> cx( .5, -.5, -.5);
        //	const FPoint<FReal> cx( .5, -.5,  .5);
        //	const FPoint<FReal> cx( .5,  .5, -.5);
        const FPoint<FReal> cx( .5,  .5,  .5);
		const FReal  wx(1.);
		FChebTensor<FReal,ORDER>::setRoots(cx, wx, x);
	
		FReal w[nnodes], f[nnodes];
		for (unsigned int n=0; n<nnodes; ++n) {
			w[n] = f[n] = FReal(n);
			//		std::cout << w[n] << "\t" << X[n] << "\t" << x[n] << std::endl;
		}



		FReal coords[3][ORDER];
		FChebTensor<FReal,ORDER>::setPolynomialsRoots(cx, wx, coords);
	
		//	for (unsigned int n=0; n<ORDER; ++n) {
		//		std::cout << coords[0][n] << "\t"
		//							<< coords[1][n] << "\t"
		//							<< coords[2][n] << std::endl;
		//	}

		FReal S[3][ORDER*ORDER];
		Interpolator.assembleInterpolator(ORDER, coords[0], S[0]);
		Interpolator.assembleInterpolator(ORDER, coords[1], S[1]);
		Interpolator.assembleInterpolator(ORDER, coords[2], S[2]);

		FReal Skron[nnodes * nnodes];
		Interpolator.assembleInterpolator(nnodes, x, Skron);



		FReal W0[nnodes];
		for (unsigned int i=0; i<nnodes; ++i) W0[i] = FReal(0.);
		applyM2M(Skron, w, nnodes, W0, nnodes);


		FReal F0[nnodes];
		for (unsigned int i=0; i<nnodes; ++i) F0[i] = FReal(0.);
		applyL2L(Skron, f, nnodes, F0, nnodes);




		unsigned int perm[3][nnodes];
		for (unsigned int i=0; i<ORDER; ++i) {
			for (unsigned int j=0; j<ORDER; ++j) {
				for (unsigned int k=0; k<ORDER; ++k) {
					const unsigned int index = k*ORDER*ORDER + j*ORDER + i;
					perm[0][index] = k*ORDER*ORDER + j*ORDER + i;
					perm[1][index] = i*ORDER*ORDER + k*ORDER + j;
					perm[2][index] = j*ORDER*ORDER + i*ORDER + k;
				}
			}
		}

		//	for (unsigned int n=0; n<nnodes; ++n)
		//		std::cout << perm[0][n] << "\t" << perm[1][n] << "\t" << perm[2][n] << std::endl;


		FReal W[nnodes];
		for (unsigned int i=0; i<nnodes; ++i) W[i] = FReal(0.);
		applym2m(S[0], w, ORDER, W, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	w[n] = W[perm[1][n]];
		applym2m(S[2], w, ORDER, W, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	w[perm[1][n]] = W[perm[2][n]];
		applym2m(S[1], w, ORDER, W, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	w[perm[2][n]] = W[n];
		FReal m2m_error(0.);
		for (unsigned int n=0; n<nnodes; ++n) {
			//std::cout << n << "\t" << w[n] << " - " << W0[n] << " = " << w[n]-W0[n] << std::endl;
			m2m_error += w[n] - W0[n];
		}

		std::cout << "------------------------------------------"
							<< "\n - M2M: ERROR = " << m2m_error << std::endl;



		FReal F[nnodes];
		for (unsigned int i=0; i<nnodes; ++i) F[i] = FReal(0.);
		applyl2l(S[0], f, ORDER, F, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	f[n] = F[perm[1][n]];
		applyl2l(S[2], f, ORDER, F, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	f[perm[1][n]] = F[perm[2][n]];
		applyl2l(S[1], f, ORDER, F, ORDER);
		for (unsigned int n=0; n<nnodes; ++n)	f[perm[2][n]] = F[n];

		FReal l2l_error(0.);
		for (unsigned int n=0; n<nnodes; ++n) {
			//std::cout << n << "\t" << f[n] << " - " << F0[n] << " = " << f[n]-F0[n] << std::endl;
			l2l_error += f[n] - F0[n];
		}
		
		std::cout << "------------------------------------------"
							<< "\n - L2L: ERROR = " << l2l_error << std::endl;
		
	}

	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	// P2M /////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
		

	const unsigned int M = 10;
	FReal points[3][M];
	FReal weights[M];
    FPoint<FReal> lp[M];
	FReal equivW[nnodes];

    { ////////////////////////////////////////////////////////
		for(unsigned int p=0; p<M; ++p){
            points[0][p] = (FReal(drand48()) - FReal(.5)) * FReal(2.);
            points[1][p] = (FReal(drand48()) - FReal(.5)) * FReal(2.);
            points[2][p] = (FReal(drand48()) - FReal(.5)) * FReal(2.);
            weights[p] = FReal(drand48());
			lp[p].setX(points[0][p]);
			lp[p].setY(points[1][p]);
			lp[p].setZ(points[2][p]);
		}
	} ////////////////////////////////////////////////////////

	
	{ // compute exact equivalent source values
		for(unsigned int i=0; i<nnodes; ++i) equivW[i] = FReal(0.);
		FReal Snorm[M * nnodes];
		Interpolator.assembleInterpolator(M, lp, Snorm);
		applyM2M(Snorm, weights, M, equivW, nnodes);
	}

	FReal W1;
	FReal W2[3][ ORDER-1];
	FReal W4[3][(ORDER-1)*(ORDER-1)];
	FReal W8[   (ORDER-1)*(ORDER-1)*(ORDER-1)];
	
	{ ////////////////////////////////////////////////////////
		W1 = FReal(0.);
		for(unsigned int i=0; i<ORDER-1; ++i)
			W2[0][i] = W2[1][i] = W2[2][i] = FReal(0.);
		for(unsigned int i=0; i<(ORDER-1)*(ORDER-1); ++i)
			W4[0][i] = W4[1][i] = W4[2][i] = FReal(0.);
		for(unsigned int i=0; i<(ORDER-1)*(ORDER-1)*(ORDER-1); ++i)
			W8[i] = FReal(0.);
		
		// loop over source particles
		for (unsigned int p=0; p<M; ++p) {
			FReal T_of_x[3][ORDER];
			T_of_x[0][0] = FReal(1.); T_of_x[0][1] = points[0][p];
			T_of_x[1][0] = FReal(1.); T_of_x[1][1] = points[1][p];
			T_of_x[2][0] = FReal(1.); T_of_x[2][1] = points[2][p];
			const FReal x2 = FReal(2.) * T_of_x[0][1]; // 1 flop
			const FReal y2 = FReal(2.) * T_of_x[1][1]; // 1 flop
			const FReal z2 = FReal(2.) * T_of_x[2][1]; // 1 flop
			for (unsigned int j=2; j<ORDER; ++j) {
				T_of_x[0][j] = x2 * T_of_x[0][j-1] - T_of_x[0][j-2]; // 2 flops
				T_of_x[1][j] = y2 * T_of_x[1][j-1] - T_of_x[1][j-2]; // 2 flops
				T_of_x[2][j] = z2 * T_of_x[2][j-1] - T_of_x[2][j-2]; // 2 flops
			}
			
			W1 += weights[p]; // 1 flop
			for (unsigned int i=1; i<ORDER; ++i) {
				const FReal wx = weights[p] * T_of_x[0][i]; // 1 flop
				const FReal wy = weights[p] * T_of_x[1][i]; // 1 flop
				const FReal wz = weights[p] * T_of_x[2][i]; // 1 flop
				W2[0][i-1] += wx; // 1 flop
				W2[1][i-1] += wy; // 1 flop
				W2[2][i-1] += wz; // 1 flop
				for (unsigned int j=1; j<ORDER; ++j) {
					const FReal wxy = wx * T_of_x[1][j]; // 1 flop
					const FReal wxz = wx * T_of_x[2][j]; // 1 flop
					const FReal wyz = wy * T_of_x[2][j]; // 1 flop
					W4[0][(j-1)*(ORDER-1) + (i-1)] += wxy; // 1 flop
					W4[1][(j-1)*(ORDER-1) + (i-1)] += wxz; // 1 flop
					W4[2][(j-1)*(ORDER-1) + (i-1)] += wyz; // 1 flop
					for (unsigned int k=1; k<ORDER; ++k) {
						const FReal wxyz = wxy * T_of_x[2][k]; // 1 flop
						W8[(k-1)*(ORDER-1)*(ORDER-1) + (j-1)*(ORDER-1) + (i-1)] += wxyz; // 1 flop
					} // flops: (ORDER-1) * 2
				} // flops: (ORDER-1) * (6 + (ORDER-1) * 2) 
			} // flops: (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1) * 2))
			
		} // flops: M * (3 + (ORDER-2) * 6 + (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1) * 2)))

	} ////////////////////////////////////////////////////////


	FReal F2[3][ORDER];
	FReal F4[3][ORDER*ORDER];
	FReal F8[   ORDER*ORDER*ORDER];

	{ ////////////////////////////////////////////////////////
		for(unsigned int i=0; i<ORDER; ++i)
			F2[0][i] = F2[1][i] = F2[2][i] = FReal(0.);
		for(unsigned int i=0; i<ORDER*ORDER; ++i)
			F4[0][i] = F4[1][i] = F4[2][i] = FReal(0.);
		for(unsigned int i=0; i<ORDER*ORDER*ORDER; ++i)
			F8[i] = FReal(0.);

		FReal T_of_y[ORDER * (ORDER-1)];
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        T_of_y[(o-1)*ORDER + j] = FReal(FChebRoots<FReal,ORDER>::T(o, FReal(FChebRoots<FReal,ORDER>::roots[j])));

		for (unsigned int l=0; l<ORDER-1; ++l)
			for (unsigned int i=0; i<ORDER; ++i) {
				F2[0][i] += T_of_y[l*ORDER + i] * W2[0][l];
				F2[1][i] += T_of_y[l*ORDER + i] * W2[1][l];
				F2[2][i] += T_of_y[l*ORDER + i] * W2[2][l];
				
				for (unsigned int m=0; m<ORDER-1; ++m)
					for (unsigned int j=0; j<ORDER; ++j) {
						F4[0][j*ORDER + i] += T_of_y[l*ORDER + i] * T_of_y[m*ORDER + j] * W4[0][m*(ORDER-1) + l];
						F4[1][j*ORDER + i] += T_of_y[l*ORDER + i] * T_of_y[m*ORDER + j] * W4[1][m*(ORDER-1) + l];
						F4[2][j*ORDER + i] += T_of_y[l*ORDER + i] * T_of_y[m*ORDER + j] * W4[2][m*(ORDER-1) + l];
							
						for (unsigned int n=0; n<ORDER-1; ++n)
							for (unsigned int k=0; k<ORDER; ++k)
								F8[k*ORDER*ORDER + j*ORDER + i] +=
									T_of_y[l*ORDER + i] * T_of_y[m*ORDER + j] * T_of_y[n*ORDER + k]	*
									W8[n*(ORDER-1)*(ORDER-1) + m*(ORDER-1) + l];
								}
			}

	} ////////////////////////////////////////////////////////


	FReal W[nnodes];
	{
		//for (unsigned int i=0; i<nnodes; ++i) W[i] = FReal(0.);
		for (unsigned int i=0; i<ORDER; ++i) {
			for (unsigned int j=0; j<ORDER; ++j) {
				for (unsigned int k=0; k<ORDER; ++k) {
					const unsigned int idx = k*ORDER*ORDER + j*ORDER + i;
					W[idx] = (W1 +
										FReal(2.) * (F2[0][i] + F2[1][j] + F2[2][k]) +
										FReal(4.) * (F4[0][j*ORDER+i] + F4[1][k*ORDER+i] + F4[2][k*ORDER+j]) +
										FReal(8.) *  F8[idx]) / (ORDER*ORDER*ORDER);
				}
			}
		}
	}

	FReal p2m_error(0.);
	for (unsigned int i=0; i<nnodes; ++i) {
		p2m_error += W[i] - equivW[i];
		//std::cout << W[i] - equivW[i] << std::endl;
	}

	std::cout << "------------------------------------------"
						<< "\n - P2M: ERROR = " << p2m_error << std::endl;



	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	// L2P /////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////


	FReal exactf[M];
	FReal f[M];
	for(unsigned int i=0; i<M; ++i) exactf[i] = f[i] = FReal(0.);


	{ // compute exact target values
		FReal Snorm[M * nnodes];
		Interpolator.assembleInterpolator(M, lp, Snorm);
		applyL2L(Snorm, W, nnodes, exactf, M);
		//for (unsigned int i=0; i<M; ++i)
		//	std::cout << exactf[i] << std::endl;
	}

	FReal f1;
	{ // sum over interpolation points
		FReal T_of_y[ORDER * (ORDER-1)];
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        T_of_y[(o-1)*ORDER + j] = FReal(FChebRoots<FReal,ORDER>::T(o, FReal(FChebRoots<FReal,ORDER>::roots[j])));
		
		// set everything to zero
		f1 = FReal(0.);
		for(unsigned int i=0; i<ORDER-1; ++i)
			W2[0][i] = W2[1][i] = W2[2][i] = FReal(0.);
		for(unsigned int i=0; i<(ORDER-1)*(ORDER-1); ++i)
			W4[0][i] = W4[1][i] = W4[2][i] = FReal(0.);
		for(unsigned int i=0; i<(ORDER-1)*(ORDER-1)*(ORDER-1); ++i)
			W8[i] = FReal(0.);

		{
			unsigned int nids[nnodes][3];
			FChebTensor<FReal,ORDER>::setNodeIds(nids);

			for (unsigned int idx=0; idx<nnodes; ++idx) {
				f1 += W[idx];
				const unsigned int i = nids[idx][0];
				const unsigned int j = nids[idx][1];
				const unsigned int k = nids[idx][2];
				
				//std::cout << i << "\t" << j << "\t" << k << std::endl;

				for (unsigned int l=0; l<ORDER-1; ++l) {
					W2[0][l] += T_of_y[l*ORDER+i] * W[idx];
					W2[1][l] += T_of_y[l*ORDER+j] * W[idx];
					W2[2][l] += T_of_y[l*ORDER+k] * W[idx];
					for (unsigned int m=0; m<ORDER-1; ++m) {
						W4[0][m*(ORDER-1)+l] += T_of_y[l*ORDER+i] * T_of_y[m*ORDER+j] * W[idx];
						W4[1][m*(ORDER-1)+l] += T_of_y[l*ORDER+i] * T_of_y[m*ORDER+k] * W[idx];
						W4[2][m*(ORDER-1)+l] += T_of_y[l*ORDER+j] * T_of_y[m*ORDER+k] * W[idx];
						for (unsigned int n=0; n<ORDER-1; ++n)
							W8[n*(ORDER-1)*(ORDER-1) + m*(ORDER-1) + l]
								+= T_of_y[l*ORDER+i]*T_of_y[m*ORDER+j]*T_of_y[n*ORDER+k] * W[idx];
					}
				}
			}

		}

	}

	{ // sum over targets
		for (unsigned int p=0; p<M; ++p) {

			FReal T_of_x[3][ORDER];
			{
				T_of_x[0][0] = FReal(1.); T_of_x[0][1] = points[0][p];
				T_of_x[1][0] = FReal(1.); T_of_x[1][1] = points[1][p];
				T_of_x[2][0] = FReal(1.); T_of_x[2][1] = points[2][p];
				const FReal x2 = FReal(2.) * T_of_x[0][1]; // 1 flop
				const FReal y2 = FReal(2.) * T_of_x[1][1]; // 1 flop
				const FReal z2 = FReal(2.) * T_of_x[2][1]; // 1 flop
				for (unsigned int j=2; j<ORDER; ++j) {
					T_of_x[0][j] = x2 * T_of_x[0][j-1] - T_of_x[0][j-2]; // 2 flops
					T_of_x[1][j] = y2 * T_of_x[1][j-1] - T_of_x[1][j-2]; // 2 flops
					T_of_x[2][j] = z2 * T_of_x[2][j-1] - T_of_x[2][j-2]; // 2 flops
				}
			}
			
			FReal f2, f4, f8;
			{
				f2 = f4 = f8 = FReal(0.);
				for (unsigned int l=1; l<ORDER; ++l) {
					f2 +=
						T_of_x[0][l] * W2[0][l-1] +
						T_of_x[1][l] * W2[1][l-1] +
						T_of_x[2][l] * W2[2][l-1];
					for (unsigned int m=1; m<ORDER; ++m) {
						f4 +=
							T_of_x[0][l] * T_of_x[1][m] * W4[0][(m-1)*(ORDER-1)+(l-1)] +
							T_of_x[0][l] * T_of_x[2][m] * W4[1][(m-1)*(ORDER-1)+(l-1)] +
							T_of_x[1][l] * T_of_x[2][m] * W4[2][(m-1)*(ORDER-1)+(l-1)];
						for (unsigned int n=1; n<ORDER; ++n)
							f8 +=
								T_of_x[0][l]*T_of_x[1][m]*T_of_x[2][n] *
								W8[(n-1)*(ORDER-1)*(ORDER-1) + (m-1)*(ORDER-1) + (l-1)];
					}
				}
			}
			f[p] = (f1 + FReal(2.)*f2 + FReal(4.)*f4 + FReal(8.)*f8) / (ORDER*ORDER*ORDER);
			
		}
	}


	FReal l2p_error(0.);
	for (unsigned int i=0; i<M; ++i) {
		l2p_error += f[i] - exactf[i];
		//std::cout << exactf[i] << "\t" << f[i] << std::endl;
	}

	std::cout << "------------------------------------------"
						<< "\n - L2P: ERROR = " << l2p_error << std::endl;




	return 0;
}


// [--END--]

		//for (unsigned int l=0; l<ORDER-1; ++l)
		//	for (unsigned int m=0; m<ORDER-1; ++m)
		//		for (unsigned int n=0; n<ORDER-1; ++n)
		//			
		//			for (unsigned int i=0; i<ORDER; ++i)
		//				for (unsigned int j=0; j<ORDER; ++j)
		//					for (unsigned int k=0; k<ORDER; ++k) {
		//
		//						const unsigned int idx = k*ORDER*ORDER + j*ORDER + i;
		//
		//						W2[0][l] += T_of_y[l*ORDER+i] * W[idx];
		//						W2[1][m] += T_of_y[m*ORDER+j] * W[idx];
		//						W2[2][n] += T_of_y[n*ORDER+k] * W[idx];
		//
		//						W4[0][m*(ORDER-1) + l] += T_of_y[l*ORDER+i] * T_of_y[m*ORDER+j] * W[idx];
		//						W4[1][n*(ORDER-1) + l] += T_of_y[l*ORDER+i] * T_of_y[n*ORDER+k] * W[idx];
		//						W4[2][n*(ORDER-1) + m] += T_of_y[m*ORDER+j] * T_of_y[n*ORDER+k] * W[idx];
		//
		//						W8[n*(ORDER-1)*(ORDER-1) + m*(ORDER-1) + l]
		//							+= T_of_y[l*ORDER+i]*T_of_y[m*ORDER+j]*T_of_y[n*ORDER+k] * W[idx];
		//					}
