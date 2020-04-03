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
// @FUSE_FFT
// ================

#include <iostream>
#include <stdlib.h>

#include <fftw3.h>

#include "Utils/FGlobal.hpp"
#include "Utils/FComplex.hpp"

#include "Utils/FTic.hpp"

#include "Utils/FParameterNames.hpp"

#include "Utils/FDft.hpp"


int main(int argc, char** argv)
{
    FHelpDescribeAndExit(argc, argv, "Test the FFT wrapper.");

    typedef float FReal;
    const FReal FRandMax = FReal(RAND_MAX);

    FTic time;

    //////////////////////////////////////////////////////////////////////////////
    // INITIALIZATION

    // dimension
    static const int dim = 3;

    // size (pick a power of 2 for better performance of the FFT algorithm)
    const int size = 50;
    int total_size = 1;
    for(int d = 0; d<dim;++d) total_size *= size;

    // input/output types
    typedef FReal           FftR2CInputType;
    typedef FComplex<FReal> FftR2COutputType;
    typedef FComplex<FReal> FftC2CInputType;
    typedef FComplex<FReal> FftC2COutputType;

    // fftw arrays
    FftR2CInputType*  FftR2CInput  = new FftR2CInputType[total_size];
    FftR2CInputType*  FftR2CInvOut = new FftR2CInputType[total_size];
    FftR2COutputType* FftR2COutput = new FftR2COutputType[total_size];

    FftC2CInputType*  FftC2CInput  = new FftC2CInputType[total_size];
    FftC2CInputType*  FftC2CInvOut = new FftC2CInputType[total_size];
    FftC2COutputType* FftC2COutput = new FftC2COutputType[total_size];


    // fftw wrappers
    typedef FFftw<FReal          ,FComplex<FReal>,dim> FftwR2CClass;
    typedef FFftw<FComplex<FReal>,FComplex<FReal>,dim> FftwC2CClass;

    std::cout<< "Init FFTW wrappers: ";
    time.tic();

    FftwR2CClass FftwR2C(size);
    FftwC2CClass FftwC2C(size);

    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    //////////////////////////////////////////////////////////////////////////////
    // EXECUTION
    // generate random physical data and set output to zero
    for( int s=0; s<total_size; ++s){

        FftR2CInput[s]  = FftR2CInputType(FReal(rand())/FRandMax);
        FftR2CInvOut[s] = FftR2CInputType(0.);
        FftR2COutput[s] = FftR2COutputType(0.,0.);

        FftC2CInput[s]  = FftC2CInputType(FReal(rand())/FRandMax,FReal(rand())/FRandMax);
        FftC2CInvOut[s] = FftC2CInputType(0.,0.);
        FftC2COutput[s] = FftC2COutputType(0.,0.);

    }

//    // display data in  physical space
//    std::cout<< "Physical data (R2C): "<<std::endl;
//    for( int s=0; s<total_size; ++s)
//      std::cout<< FftR2CInput[s] << ", ";
//    std::cout<<std::endl;
//
//    std::cout<< "Physical data (C2C): "<<std::endl;
//    for( int s=0; s<total_size; ++s)
//      std::cout<< FftC2CInput[s] << ", ";
//    std::cout<<std::endl;

    // perform fft
    std::cout<< "Perform Forward FFT: ";
    time.tic();

    FftwR2C.applyDFT(FftR2CInput,FftR2COutput);
    FftwC2C.applyDFT(FftC2CInput,FftC2COutput);

    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

//    // display transform in Fourier space
//    // beware the real data FFT stores only N/2+1 complex output values
//    std::cout<< "Transformed data (R2C): "<<std::endl;
//    for( int s=0; s<total_size; ++s)
//      std::cout<< FftR2COutput[s] << ", ";
//    std::cout<<std::endl;
//
//    std::cout<< "Transformed data (C2C): "<<std::endl;
//    for( int s=0; s<total_size; ++s)
//      std::cout<< FftC2COutput[s] << ", ";
//    std::cout<<std::endl;

    // perform ifft of transformed data (in order to get physical values back)
    std::cout<< "Perform Normalized Backward FFT: ";
    time.tic();

    FftwR2C.applyIDFTNorm(FftR2COutput,FftR2CInvOut);
    FftwC2C.applyIDFTNorm(FftC2COutput,FftC2CInvOut);

    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

//    // display data in physical space
//    std::cout<< "Physical data from 1/N*IFFT(FFT(Physical data)) (R2C) : "<<std::endl;
//    for( int s=0; s<total_size; ++s)
//      std::cout<< FftR2CInvOut[s] << ", ";
//    std::cout<<std::endl;
//
//    std::cout<< "Physical data from 1/N*IFFT(FFT(Physical data)) (C2C) : "<<std::endl;
//    for( int s=0; s<total_size; ++s)
//      std::cout<< FftC2CInvOut[s] << ", ";
//    std::cout<<std::endl;

    //////////////////////////////////////////////////////////////////////////////
    // VALIDATION
    bool isWorking = true;
    FReal threshold = (sizeof(FReal)==8 ? FReal(1.e-12) : FReal(1.e-6));

    FReal* FftC2CInvOutReal = new FReal[total_size]; FReal* FftC2CInputReal = new FReal[total_size];
    FReal* FftC2CInvOutImag = new FReal[total_size]; FReal* FftC2CInputImag = new FReal[total_size];
    for( int s=0; s<total_size; ++s){
        FftC2CInvOutReal[s] = FftC2CInvOut[s].getReal(); FftC2CInputReal[s] = FftC2CInput[s].getReal();
        FftC2CInvOutImag[s] = FftC2CInvOut[s].getImag(); FftC2CInputImag[s] = FftC2CInput[s].getImag();
    }

    if(FMath::FAccurater<FReal>(FftR2CInput,     FftR2CInvOut,     total_size).getRelativeInfNorm() > threshold) isWorking=false;
    if(FMath::FAccurater<FReal>(FftR2CInput,     FftR2CInvOut,     total_size).getRelativeL2Norm() > threshold) isWorking=false;
    if(FMath::FAccurater<FReal>(FftC2CInputReal, FftC2CInvOutReal, total_size).getRelativeInfNorm() > threshold) isWorking=false;
    if(FMath::FAccurater<FReal>(FftC2CInputImag, FftC2CInvOutImag, total_size).getRelativeInfNorm() > threshold) isWorking=false;
    if(FMath::FAccurater<FReal>(FftC2CInputReal, FftC2CInvOutReal, total_size).getRelativeL2Norm() > threshold) isWorking=false;
    if(FMath::FAccurater<FReal>(FftC2CInputImag, FftC2CInvOutImag, total_size).getRelativeL2Norm() > threshold) isWorking=false;

    if(isWorking) std::cout << "FFT Wrapper is working (errors < threshold = " << threshold << ")." << std::endl;
    else std::cout << "FFT Wrapper is NOT working (errors > threshold = " << threshold << ")." << std::endl;

    //free memory
    delete[] FftR2CInput;
    delete[] FftR2CInvOut;
    delete[] FftR2COutput;
    delete[] FftC2CInput;
    delete[] FftC2CInvOut;
    delete[] FftC2COutput;
    delete[] FftC2CInvOutReal;
    delete[] FftC2CInputReal;
    delete[] FftC2CInvOutImag;
    delete[] FftC2CInputImag;

}// end test
