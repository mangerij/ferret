#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  //SSE2
#include <pmmintrin.h> //SSE3
#ifdef __SSSE3__
#include <tmmintrin.h>  //SSSE3
#endif
#ifdef __SSSE4_1__
#include <smmintrin.h> // SSE4
#endif

int main() {
	__m128d tx, ty ;
	tx += ty ;
  return 0;
}
