int main(){

    int i ;
#ifdef __INTEL_COMPILER

     i = 0;

#else

#error 'Not Intel Compiler "

#endif
}
