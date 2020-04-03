
#include <x86intrin.h>
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2
#include <pmmintrin.h> // SSE3
#ifdef __SSSE3__
#include <tmmintrin.h>  //SSSE3
#endif
#ifdef __SSSE4_1__
#include <smmintrin.h> // SSE4
#endif
int main(){
	__m128d res0d, res1d;
	res0d = _mm_hadd_pd(res0d, res1d);

	__m128 res0, res1;
	res0 = _mm_hadd_ps(res0, res1);

	return 0;
}
