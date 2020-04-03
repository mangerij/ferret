// ===================================================================================
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

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include  "ScalFmmConfig.h"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Kernels/P2P/FP2P.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Utils/FParameterNames.hpp"
//
/// \file  DirectComputation.cpp
//!
//! \brief DirectComputation: Driver to compute direct interaction between N particles for 1/r kernel.
//!
//! DirectComputation: Driver to compute direct interaction between N particles for 1/r kernel.
//! the particles are read from file given by -fin argument and potential, forces are stored in FMA format.
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary).
//!                             Only our FMA (.bma, .bfma) is allowed "
//!     \param    -fout filenameOUT   output file  with extension (default output.bfma)
//!      \param   -verbose : print index x y z Q V fx fy fz
//!

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         ">> This executable has to be used to compute  interaction either for periodic or non periodic system.\n"
                         ">> Example -fin filenameIN.{fma or bfma)     -fout filenameOUT{fma or bfma) \n"
                         ">> Default input file : Data/unitCubeXYZQ20k.fma\n"
                         ">> Only our FMA (.bma, .bfma) is allowed as input.\n"
                         ">> Output file  with extension (default output.bfma).",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::EnabledVerbose);

    //////////////////////////////////////////////////////////////
    typedef double FReal;
    const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/unitCubeXYZQ20k.fma");
    const std::string filenameIn(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options,  defaultFile.c_str()));
    const std::string filenameOut(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "output.bfma"));
    //
    FTic counter;

    // -----------------------------------------------------
    //  LOADER
    //  -----------------------------------------------------
    // ---------------------------------------------------------------------------------
    // Read  particles in the Octree
    // ---------------------------------------------------------------------------------
    std::cout << "Opening : " << filenameIn << "\n";
    //
    FFmaGenericLoader<FReal> loader(filenameIn);
    //
    FSize nbParticles = static_cast<int>(loader.getNumberOfParticles());
    std::cout << "Read " << nbParticles << " particles ..." << std::endl;
    double BoxWith=loader.getBoxWidth();
    FPoint<FReal> Centre(loader.getCenterOfBox().getX(), loader.getCenterOfBox().getY() , loader.getCenterOfBox().getZ());
    std::cout << "\tWidth : " <<BoxWith << " \t center x : " << loader.getCenterOfBox().getX()
              << " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;

    counter.tic();
    //
    FmaRWParticle<FReal, 4,8> *  particles = new FmaRWParticle<FReal, 4,8>[nbParticles];
    memset(particles, 0, sizeof(FmaRWParticle<FReal, 4,8>) * nbParticles) ;
    //
    double totalCharge = 0.0;
    //
    //	int nbDataToRead = particles[0].getReadDataNumber();
    for(int idx = 0 ; idx<nbParticles ; ++idx){
        //
        loader.fillParticle(particles[idx].getPtrFirstData(), particles[idx].getReadDataNumber());
        //	loader.fillParticle(particles[idx].getPtrFirstData(), nbDataToRead);    // OK
        //  loader.fillParticle(particles[idx]); // OK
        //    std::cout << idx <<"  "<<  particles[idx].getPosition() << " "<<particles[idx].getPhysicalValue() << " "<<particles[idx].getPotential()
        //			<<"  " << particles[idx].getForces()[0]<<"  " <<particles[idx].getForces()[1]<<"  " <<particles[idx].getForces()[2]<<"  " <<std::endl;
        //
        totalCharge += particles[idx].getPhysicalValue() ;
    }

    counter.tac();

    std::cout << std::endl;
    std::cout << "Total Charge         = "<< totalCharge <<std::endl;
    std::cout << std::endl;

    std::cout << "Done  " << "(@ reading Particles  " << counter.elapsed() << " s)." << std::endl;
    //
    // ----------------------------------------------------------------------------------------------------------
    //                                   COMPUTATION
    // ----------------------------------------------------------------------------------------------------------
    // interaction kernel evaluator
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    const MatrixKernelClass MatrixKernel;
    FReal denergy = 0.0;
    //
    //  computation
    //
    {
        printf("Compute :\n");
        counter.tic();
#pragma omp parallel shared(nbParticles, particles,denergy)
        {
#pragma omp for
            for(int idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
                //
                // compute with all other except itself
                //
                // Compute force and potential between  particles[idxTarget] and particles inside the box
                //
                for(int idxOther = 0; idxOther < nbParticles ; ++idxOther){
                    if( idxOther != idxTarget ){
                        FP2P::NonMutualParticles(
                                    particles[idxOther].getPosition().getX(), particles[idxOther].getPosition().getY(),
                                    particles[idxOther].getPosition().getZ(),particles[idxOther].getPhysicalValue(),
                                    particles[idxTarget].getPosition().getX(), particles[idxTarget].getPosition().getY(),
                                    particles[idxTarget].getPosition().getZ(),particles[idxTarget].getPhysicalValue(),
                                    &particles[idxTarget].setForces()[0],&particles[idxTarget].setForces()[1],
                                &particles[idxTarget].setForces()[2],particles[idxTarget].setPotential(),&MatrixKernel);
                    }
                }
            } // end for
            // Compute the energy
#pragma omp  for reduction(+:denergy)
            for(int idx = 0 ; idx < nbParticles ; ++idx){
                denergy +=  particles[idx].getPotential()*(particles[idx].getPhysicalValue()) ;
            }
        } // end pragma parallel
        //
        denergy *= 0.5 ;
        counter.tac();
        //
        printf("Energy =   %.14e\n",denergy);
        //
        std::cout << "Done  " << "(@ Direct computation done = " << counter.elapsed() << " s)." << std::endl;
        std::cout << "\n"<< "END  "
                  << "-------------------------------------------------------------------------"
                  << std::endl << std::endl ;
    } // END

    //
    // ----------------------------------------------------------------
    //  Save  computation in binary format
    //
    //

    std::cout << "Generate " << filenameOut <<"  for output file" << std::endl;
    //
    std::cout << " nbParticles: " << nbParticles <<"  " << sizeof(nbParticles) <<std::endl;
    std::cout << " denergy: " << denergy <<"  " << sizeof(denergy) <<std::endl;
    std::cout << " Box size: " << loader.getBoxWidth() << "  " << sizeof(loader.getBoxWidth())<<std::endl;
    //
    FFmaGenericWriter<FReal> writer(filenameOut) ;
    writer.writeHeader(Centre,BoxWith, nbParticles,*particles) ;
    writer.writeArrayOfParticles(particles, nbParticles);
    //
    // end generate
    // -----------------------------------------------------
    //
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::EnabledVerbose.options)){
        denergy = 0 ;
        for(int idx = 0 ; idx < nbParticles ; ++idx){
            std::cout << ">> index " << idx << std::endl;
            std::cout << " x   " << particles[idx].getPosition().getX() << " y  " << particles[idx].getPosition().getY() << " z  " << particles[idx].getPosition().getZ() << std::endl;
            std::cout << " Q   " << particles[idx].getPhysicalValue()   << " V  " << particles[idx].getPotential() << std::endl;
            std::cout << " fx  " << particles[idx].getForces()[0]       << " fy " << particles[idx].getForces()[1]       << " fz " << particles[idx].getForces()[2] << std::endl;
            std::cout << "\n";
            denergy +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
        }
    }
    std::cout << " ENERGY " << denergy << std::endl;
    //
    delete[] particles;
    return 0;
}



