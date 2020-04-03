// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#ifndef FBLAS_HPP
#define FBLAS_HPP

#include "FGlobal.hpp"

#ifndef SCALFMM_USE_BLAS
#error The BLAS header is included while SCALFMM_USE_BLAS is turned OFF
#endif

// This file interfaces the blas functions
// to enable a generic use.
// If no blas has been enabled in the cmake,
// the function will be empty


// for real
const double D_ZERO =  0.0;
const double D_ONE  =  1.0;
const double D_MONE = -1.0;
const float  S_ZERO =  0.0;
const float  S_ONE  =  1.0;
const float  S_MONE = -1.0;
// for complex
const double Z_ZERO[2] =  {0.0,0.0};
const double Z_ONE[2]  =  {1.0,0.0};
const double Z_MONE[2] =  {-1.0,0.0};
const float  C_zero[2] =  {0.0,0.0};
const float  C_one[2]  =  {1.0,0.0};
const float  C_MONE[2] =  {-1.0,0.0};

//const double D_PREC = 1e-16;

const unsigned N_ONE = 1;
const int N_MONE = -1;
const char JOB_STR[] = "NTOSVULCR";


extern "C"
{
	// double //////////////////////////////////////////////////////////
	// blas 1
	double ddot_(const unsigned*, const double*, const unsigned*, const double*, const unsigned*);
	void dscal_(const unsigned*, const double*, const double*, const unsigned*);
	void dcopy_(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
	void daxpy_(const unsigned*, const double*, const double*, const unsigned*, double*, const unsigned*);
	// blas 2
	void dgemv_(const char*, const unsigned*, const unsigned*, const double*,
							const double*, const unsigned*, const double*, const unsigned*,
							const double*, double*, const unsigned*);
	// blas 3
	void dgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const double*, double*, const unsigned*,
							double*, const unsigned*, const double*, double*,	const unsigned*);
	// lapack
	void dgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 double*, const unsigned*, double*, double*, const unsigned*,
							 double*, const unsigned*, double*, const unsigned*, int*);
	void dgeqrf_(const unsigned*, const unsigned*, double*, const unsigned*,
							 double*, double*, const unsigned*, int*);
	void dorgqr_(const unsigned*, const unsigned*, const unsigned*,
							 double*, const unsigned*, double*, double*, const unsigned*, int*);
	void dormqr_(const char*, const char*, 
               const unsigned*, const unsigned*, const unsigned*,
							 const double*, const unsigned*, 
               double*, double*, const unsigned*, 
               double*, const unsigned*, int*);
    void dpotrf_(const char*, const unsigned*, double*, const unsigned*, int*);

#ifdef SCALFMM_USE_MKL_AS_BLAS
  // mkl: hadamard product is not implemented in mkl_blas
  void vdmul_(const unsigned* n, const double*, const double*, double*);
  void vsmul_(const unsigned* n, const float*, const float*, float*);
  void vzmul_(const unsigned* n, const double*, const double*, double*);
  void vcmul_(const unsigned* n, const float*, const float*, float*);
#else
  // TODO create interface for hadamard product in case an external Blas is used
#endif

	// single //////////////////////////////////////////////////////////
	// blas 1
	float sdot_(const unsigned*, const float*, const unsigned*,	const float*, const unsigned*);
	void sscal_(const unsigned*, const float*, const float*, const unsigned*);
	void scopy_(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);
	void saxpy_(const unsigned*, const float*, const float*, const unsigned*, float*, const unsigned*);
	// blas 2
	void sgemv_(const char*, const unsigned*, const unsigned*, const float*,
							const float*, const unsigned*, const float*, const unsigned*,
							const float*, float*, const unsigned*);
	// blas 3
	void sgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const float*, float*, const unsigned*,
							float*, const unsigned*, const float*, float*, const unsigned*);
	// lapack
	void sgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 float*, const unsigned*, float*, float*, const unsigned*,
							 float*, const unsigned*, float*, const unsigned*, int*);
	void sgeqrf_(const unsigned*, const unsigned*, float*, const unsigned*,
							 float*, float*, const unsigned*, int*);
	void sorgqr_(const unsigned*, const unsigned*, const unsigned*,
							 float*, const unsigned*, float*, float*, const unsigned*, int*);
	void sormqr_(const char*, const char*, 
               const unsigned*, const unsigned*, const unsigned*,
							 const float*, const unsigned*, 
               float*, float*, const unsigned*, 
               float*, const unsigned*, int*);
    void spotrf_(const char*, const unsigned*, float*, const unsigned*, int*);

	// double complex //////////////////////////////////////////////////
	// blas 1
	void zscal_(const unsigned*, const double*, const double*, const unsigned*);
	void zcopy_(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
	void zaxpy_(const unsigned*, const double*, const double*, const unsigned*, double*, const unsigned*);
	// blas 2
	void zgemv_(const char*, const unsigned*, const unsigned*, const double*,
							const double*, const unsigned*, const double*, const unsigned*,
							const double*, double*, const unsigned*);
	// blas 3
	void zgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const double*, double*, const unsigned*,
							double*, const unsigned*, const double*, double*, const unsigned*);
	void zgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 double*, const unsigned*, double*, double*, const unsigned*,
							 double*, const unsigned*, double*,   int*,  double*,   int*);

	void zgeqrf_(const unsigned*, const unsigned*, double*, const unsigned*,
							 double*, double*, const unsigned*, int*);
    void zpotrf_(const char*, const unsigned*, double*, const unsigned*, int*);

	// single complex //////////////////////////////////////////////////
	// blas 1
	void cscal_(const unsigned*, const float*, const float*, const unsigned*);
	void ccopy_(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);
	void caxpy_(const unsigned*, const float*, const float*, const unsigned*, float*, const unsigned*);
	// blas 2
	void cgemv_(const char*, const unsigned*, const unsigned*, const float*,
							const float*, const unsigned*, const float*, const unsigned*,
							const float*, float*, const unsigned*);
	// blas 3
	void cgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const float*, float*, const unsigned*,
							float*, const unsigned*, const float*, float*, const unsigned*);
	void cgeqrf_(const unsigned*, const unsigned*, float*, const unsigned*,
							 float*, float*, const unsigned*, int*);
    void cpotrf_(const char*, const unsigned*, float*, const unsigned*, int*);


}


// Hadamard (i.e. entrywise) product: 
// NB: The following optimized routines are currently not used 
// since they have not proved their efficiency in comparison 
// with a naive application of the entrywise product.
#ifdef SCALFMM_USE_MKL_AS_BLAS

//#include "mkl_vml.h" 

namespace FMkl{

  // Hadamard product: dest[i]=a[i]*b[i]
  inline void had(const unsigned n, const double* const a, const double* const b, double* const dest)
  { vdmul_(&n, a, b, dest); }
  inline void had(const unsigned n, const float* const a, const float* const b, float* const dest)
  { vsmul_(&n, a, b, dest); }
  inline void c_had(const unsigned n, const double* const a, const double* const b, double* const dest)
  { vzmul_(&n, a, b, dest); }
  inline void c_had(const unsigned n, const float* const a, const float* const b, float* const dest)
  { vcmul_(&n, a, b, dest); }

}

#else

namespace FBlas{

  // TODO create interface for Hadamard product in case an external Blas is used

}

#endif
// end Hadamard product



namespace FBlas {

	// copy
	inline void copy(const unsigned n, double* orig, double* dest)
	{	dcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void copy(const unsigned n, const double* orig, double* dest)
	{	dcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void copy(const unsigned n, float* orig, float* dest)
	{	scopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
    inline void copy(const unsigned n, const float* orig, float* dest)
    {   scopy_(&n, orig, &N_ONE, dest, &N_ONE); }
	inline void c_copy(const unsigned n, double* orig, double* dest)
	{	zcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void c_copy(const unsigned n, const double* orig, double* dest)
	{	zcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void c_copy(const unsigned n, float* orig, float* dest)
	{	ccopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
    inline void c_copy(const unsigned n, const float* orig, float* dest)
    {   ccopy_(&n, orig, &N_ONE, dest, &N_ONE); }

	// copy (variable increment)
	inline void copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
	{	dcopy_(&n, orig, &inco, dest, &incd);	}
	inline void copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
	{	scopy_(&n, orig, &inco, dest, &incd);	}
	inline void c_copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
	{	zcopy_(&n, orig, &inco, dest, &incd);	}
	inline void c_copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
	{	ccopy_(&n, orig, &inco, dest, &incd);	}

	// scale
	inline void scal(const unsigned n, const double d, double* const x)
	{	dscal_(&n, &d, x, &N_ONE); }
	inline void scal(const unsigned n, const float d, float* const x)
	{	sscal_(&n, &d, x, &N_ONE); }
	inline void c_scal(const unsigned n, const double d, double* const x)
	{	zscal_(&n, &d, x, &N_ONE); }
	inline void c_scal(const unsigned n, const float d, float* const x)
	{	cscal_(&n, &d, x, &N_ONE); }

	// scale (variable increment)
	inline void scal(const unsigned n, const double d, double* const x, const unsigned incd)
	{	dscal_(&n, &d, x, &incd); }
	inline void scal(const unsigned n, const float d, float* const x, const unsigned incd)
	{	sscal_(&n, &d, x, &incd); }
	inline void c_scal(const unsigned n, const double d, double* const x, const unsigned incd)
	{	zscal_(&n, &d, x, &incd); }
	inline void c_scal(const unsigned n, const float d, float* const x, const unsigned incd)
	{	cscal_(&n, &d, x, &incd); }

	// set zero
	inline void setzero(const unsigned n, double* const x)
	{	for (unsigned i=0; i<n; ++i) x[i] = 0.0; }
	inline void setzero(const unsigned n, float* const x)
	{	for (unsigned i=0; i<n; ++i) x[i] = 0.0f; }
	inline void c_setzero(const unsigned n, double* const x)
	{	for (unsigned i=0; i<n; ++i) x[i*2] = x[i*2+1] = 0.0; }
	inline void c_setzero(const unsigned n, float* const x)
	{	for (unsigned i=0; i<n; ++i) x[i*2] = x[i*2+1] = 0.0f; }

	// y += x
	inline void add(const unsigned n, double* const x, double* const y)
	{	daxpy_(&n, &D_ONE, x, &N_ONE, y, &N_ONE);	}
	inline void add(const unsigned n, float* const x, float* const y)
	{	saxpy_(&n, &S_ONE, x, &N_ONE, y, &N_ONE);	}
	inline void c_add(const unsigned n, float* const x, float* const y)
	{	caxpy_(&n, C_one, x, &N_ONE, y, &N_ONE);	}
	inline void c_add(const unsigned n, double* const x,double* const y)
	{	zaxpy_(&n, Z_ONE, x, &N_ONE, y, &N_ONE);	}

	// y += d x
	inline void axpy(const unsigned n, const double d, const double* const x, double* const y)
	{	daxpy_(&n, &d, x, &N_ONE, y, &N_ONE);	}
	inline void axpy(const unsigned n, const float d, const float* const x, float* const y)
	{	saxpy_(&n, &d, x, &N_ONE, y, &N_ONE);	}
	inline void c_axpy(const unsigned n, const float* d, const float* const x, float* const y)
	{	caxpy_(&n, d, x, &N_ONE, y, &N_ONE);	}
	inline void c_axpy(const unsigned n, const double* d, const double* const x, double* const y)
	{	zaxpy_(&n, d, x, &N_ONE, y, &N_ONE);	}



	//	// y = d Ax
	//	inline void gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	//	{	cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, D_ZERO, y, N_ONE); }
	//	inline void gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	//	{	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, S_ZERO, y, N_ONE); }
	// y = d Ax
	inline void gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE); }
	inline void gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE); }
	inline void c_gemv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, C_zero, y, &N_ONE); }
	inline void c_gemv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, Z_ZERO, y, &N_ONE); }

	//	// y += d Ax
	//	inline void gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	//	{	cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, D_ONE, y, N_ONE); }
	//	inline void gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	//	{	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, S_ONE, y, N_ONE); }
	// y += d Ax
	inline void gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);	}
	inline void gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);	}
	inline void c_gemva(const unsigned m, const unsigned n, const float* d, const float* A, const float *x, float *y)
	{	cgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, C_one, y, &N_ONE);	}
	inline void c_gemva(const unsigned m, const unsigned n, const double* d, const double* A, const double *x, double *y)
	{	zgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, Z_ONE, y, &N_ONE);	}

	//	// y = d A^T x
	//	inline void gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	//	{ cblas_dgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, D_ZERO, y, N_ONE); }
	//	inline void gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	//	{	cblas_sgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, S_ZERO, y, N_ONE); }
	// y = d A^T x
	inline void gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE); }
	inline void gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE); }
	inline void c_gemtv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, C_zero, y, &N_ONE); }
	inline void c_gemtv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, Z_ZERO, y, &N_ONE); }
	inline void c_gemhv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, C_zero, y, &N_ONE); } // hermitian transposed
	inline void c_gemhv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, Z_ZERO, y, &N_ONE); } // hermitian transposed

	//	// y += d A^T x
	//	inline void gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	//	{	cblas_dgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, D_ONE, y, N_ONE); }
	//	inline void gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	//	{	cblas_sgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, S_ONE, y, N_ONE); }
	// y += d A^T x
	inline void gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);	}
	inline void gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);	}
	inline void c_gemtva(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, C_one, y, &N_ONE);	}
	inline void c_gemtva(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, Z_ONE, y, &N_ONE); }
	inline void c_gemhva(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, C_one, y, &N_ONE);	} // hermitian transposed
	inline void c_gemhva(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, Z_ONE, y, &N_ONE);	} // hermitian transposed




	// C = d A B, A is m x p, B is p x n
	inline void gemm(unsigned m, unsigned p, unsigned n, double d,
									 double* A, unsigned ldA, double* B, unsigned ldB, double* C, unsigned ldC)
	{	dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &D_ZERO, C, &ldC);	}
	inline void gemm(unsigned m, unsigned p, unsigned n, float d,
									 float* A, unsigned ldA, float* B, unsigned ldB, float* C, unsigned ldC)
	{	sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &S_ZERO, C, &ldC);	}
	inline void c_gemm(const unsigned m, const unsigned p, const unsigned n, const float* d,
										 float* A, const unsigned ldA, float* B, const unsigned ldB, float* C, const unsigned ldC)
	{
		cgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, C_zero, C, &ldC);	}
	inline void c_gemm(const unsigned m, const unsigned p, const unsigned n, const double* d,
										 double* A, const unsigned ldA, double* B, const unsigned ldB, double* C, const unsigned ldC)
	{
		zgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}

	// C += d A B, A is m x p, B is p x n
	inline void gemma(unsigned m, unsigned p, unsigned n, double d,
										double* A, unsigned ldA, double* B, unsigned ldB,	double* C, unsigned ldC)
	{	dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &D_ONE, C, &ldC); }
	inline void gemma(unsigned m, unsigned p, unsigned n, float d,
										float* A, unsigned ldA, float* B, unsigned ldB,	float* C, unsigned ldC)
	{	sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &S_ONE, C, &ldC); }
	inline void c_gemma(unsigned m, unsigned p, unsigned n, float* d,
											float* A, unsigned ldA, float* B, unsigned ldB,	float* C, unsigned ldC)
	{	cgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, C_one, C, &ldC); }
	inline void c_gemma(unsigned m, unsigned p, unsigned n, double* d,
											double* A, unsigned ldA, double* B, unsigned ldB,	double* C, unsigned ldC)
	{	zgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }

	// C = d A^T B, A is m x p, B is m x n
	inline void gemtm(unsigned m, unsigned p, unsigned n, double d,
										double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
	{	dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &D_ZERO, C, &ldC);	}
	inline void gemtm(unsigned m, unsigned p, unsigned n, float d,
										float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
	{	sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &S_ZERO, C, &ldC);	}
	inline void c_gemtm(unsigned m, unsigned p, unsigned n, float* d,
											float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
	{	cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_zero, C, &ldC);	}
	inline void c_gemtm(unsigned m, unsigned p, unsigned n, double* d,
											double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
	{	zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}
	inline void c_gemhm(unsigned m, unsigned p, unsigned n, float* d, // hermitialn transposed
											float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
	{	cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_zero, C, &ldC);	}
	inline void c_gemhm(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
											double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
	{	zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}

	// C += d A^T B, A is m x p, B is m x n
	inline void gemtma(unsigned m, unsigned p, unsigned n, double d,
										 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &D_ONE, C, &ldC); }
	inline void gemtma(unsigned m, unsigned p, unsigned n, float d,
										 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &S_ONE, C, &ldC); }
	inline void c_gemtma(unsigned m, unsigned p, unsigned n, float* d,
											 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_one, C, &ldC); }
	inline void c_gemtma(unsigned m, unsigned p, unsigned n, double* d,
											 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }
	inline void c_gemhma(unsigned m, unsigned p, unsigned n, float* d, // hermitian transposed
											 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_one, C, &ldC); }
	inline void c_gemhma(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
											 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }
	
	
	// C = d A B^T, A is m x p, B is n x p
	inline void gemmt(unsigned m, unsigned p, unsigned n, double d,
										double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &D_ZERO, C, &ldC);	}
	inline void gemmt(unsigned m, unsigned p, unsigned n, float d,
										float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
	{	sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &S_ZERO, C, &ldC);	}
	inline void c_gemmt(unsigned m, unsigned p, unsigned n, float d,
											float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	cgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, C_zero, C, &ldC); }
	inline void c_gemmt(unsigned m, unsigned p, unsigned n, double d,
											double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
	{	zgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}
	inline void c_gemmh(unsigned m, unsigned p, unsigned n, float d, // hermitian transposed
											float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB, C_zero, C, &ldC); }
	inline void c_gemmh(unsigned m, unsigned p, unsigned n, double d, // hermitian transposed
											double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
	{	zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}

	// C += d A B^T, A is m x p, B is n x p
	inline void gemmta(unsigned m, unsigned p, unsigned n, double d,
										 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &D_ONE, C, &ldC); }
	inline void gemmta(unsigned m, unsigned p, unsigned n, float d,
										 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &S_ONE, C, &ldC); }
	inline void c_gemmta(unsigned m, unsigned p, unsigned n, float* d,
											 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	cgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, d, A, &ldA, B, &ldB, C_one, C, &ldC); }
	inline void c_gemmta(unsigned m, unsigned p, unsigned n, double* d,
											 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	zgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }
	inline void c_gemmha(unsigned m, unsigned p, unsigned n, float* d, // hermitian transposed
											 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, d, A, &ldA, B, &ldB, C_one, C, &ldC); }
	inline void c_gemmha(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
											 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }


	// singular value decomposition
    //
	inline int gesvd(unsigned m, unsigned n, double* A, double* S, double* VT, unsigned ldVT,
									 unsigned nwk, double* wk)
	{
		int INF;
		dgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
		return INF;
	}
	//
	//    A = U * SIGMA * conjugate-transpose(V)
	// JOB_STR+2 = 'O':  the first min(m,n) columns of U (the left singular vectors) are overwritten on the array A;
	inline int c_gesvd(unsigned m, unsigned n, double* A, double* S, double* VT, unsigned ldVT,
			int& nwk, double* wk,double* rwk)
	{
		int INF;
		zgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, rwk,&INF);
		return INF;
	}

	inline int gesvd(unsigned m, unsigned n, float* A, float* S, float* VT, unsigned ldVT,
									 unsigned nwk, float* wk)
	{
		int INF;
		sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
		return INF;
	}

    // singular value decomposition (SO)
    inline int gesvdSO(unsigned m, unsigned n, double* A, double* S, double* U, unsigned ldU,
                       unsigned nwk, double* wk)
    {
        int INF;
        dgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
        return INF;
    }
    inline int gesvdSO(unsigned m, unsigned n, float* A, float* S, float* U, unsigned ldU,
                       unsigned nwk, float* wk)
    {
        int INF;
        sgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
        return INF;
    }

    // singular value decomposition (AA)
    inline int gesvdAA(unsigned m, unsigned n, double* A, double* S, double* U, unsigned ldU,
                       unsigned nwk, double* wk)
    {
        int INF;
        dgesvd_("A", "A", &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
        return INF;
    }
    inline int gesvdAA(unsigned m, unsigned n, float* A, float* S, float* U, unsigned ldU,
                       unsigned nwk, float* wk)
    {
        int INF;
        sgesvd_("A", "A", &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
        return INF;
    }

	// Scalar product v1'*v2
	inline double scpr(const unsigned n, const double* const v1, const double* const v2)
	{	return ddot_(&n, v1, &N_ONE, v2, &N_ONE); }
	inline float scpr(const unsigned n, const float* const v1, const float* const v2)
	{	return sdot_(&n, v1, &N_ONE, v2, &N_ONE);	}



	// QR factorisation
	inline int geqrf(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
	{
		int INF;
		dgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}
	inline int geqrf(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
	{
		int INF;
		sgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}
	
	inline int c_geqrf(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
	{
		int INF;
		cgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}
	
	inline int c_geqrf(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
	{
		int INF;
		zgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}


	// return full of Q-Matrix (QR factorization) in A
	inline int orgqr_full(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
	{
		int INF;
		dorgqr_(&m, &m, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}
	inline int orgqr_full(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
	{
		int INF;
		sorgqr_(&m, &m, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}
	// return the leading n columns of Q-Matrix (QR factorization) in A
	inline int orgqr(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
	{
		int INF;
		dorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}
	inline int orgqr(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
	{
		int INF;
		sorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
		return INF;
	}



    // apply Q-Matrix (from QR factorization) to C
    // LEFT: Q(^T)C
	inline int left_ormqr(const char* TRANS, const unsigned m, const unsigned n, const double* A, double* tau, double* C, unsigned nwk, double* wk)
	{
		int INF;
		dormqr_("L", TRANS, &m, &n, &m, A, &m, tau, C, &m, wk, &nwk, &INF);
		return INF;
	}
	inline int left_ormqr(const char* TRANS, const unsigned m, const unsigned n, const float* A, float* tau, float* C, unsigned nwk, float* wk)
	{
		int INF;
		sormqr_("L", TRANS, &m, &n, &m, A, &m, tau, C, &m, wk, &nwk, &INF);
		return INF;
	}
    // RIGHT: CQ(^T)
	inline int right_ormqr(const char* TRANS, const unsigned m, const unsigned n, const double* A, double* tau, double* C, unsigned nwk, double* wk)
	{
		int INF;
		dormqr_("R", TRANS, &m, &n, &n, A, &n, tau, C, &m, wk, &nwk, &INF);
		return INF;
	}
	inline int right_ormqr(const char* TRANS, const unsigned m, const unsigned n, const float* A, float* tau, float* C, unsigned nwk, float* wk)
	{
		int INF;
		sormqr_("R", TRANS, &m, &n, &n, A, &n, tau, C, &m, wk, &nwk, &INF);
		return INF;
	}

    // Cholesky decomposition: A=LL^T (if A is symmetric definite positive)
    inline int potrf(const unsigned m, double* A, const unsigned n)
    { 
		int INF;  
        dpotrf_("L", &m, A, &n, &INF); 
        return INF;
    }
    inline int potrf(const unsigned m, float* A, const unsigned n)
    { 
		int INF;  
        spotrf_("L", &m, A, &n, &INF); 
        return INF;
    }

} // end namespace FCBlas

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//#else
//enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
//enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};
//
//template <typename T>
//void cblas_gemv(const CBLAS_ORDER order ,
//								const CBLAS_TRANSPOSE TransA , const int M , const int N ,
//								const void *alpha , const void *A , const int lda ,
//								const void *X , const int incX , const void *beta ,
//								void *Y , const int incY){
//}
//template <typename T>
//void cblas_dotu_sub( const int N , const void *X , const int incX ,
//										 const void *Y , const int incY , void *dotu){
//}


#endif //FBLAS_HPP

