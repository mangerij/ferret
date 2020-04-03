
#include "immintrin.h"


int main() {
#ifdef __MIC__
	__m512 tx, ty ;
	tx += ty ;
#endif
  return 0;
}

