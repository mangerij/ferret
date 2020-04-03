/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
//
#include "Files/FFmaGenericLoader.hpp"
#include "Files/FDlpolyLoader.hpp"
//
#include "Utils/FGlobal.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FParameters.hpp"
#include "Files/FGenerateDistribution.hpp"
#include "Files/FExportWriter.hpp"

#include "Utils/FParameterNames.hpp"

//
/// \file  changeFmaFormat.cpp
//!
//! \brief changeFormat: Driver to transform a FMA format and to build a visualization file
//!
//!  Driver to transform a FMA format and/or to build a visualization file<br>
//! For a description of the FMA format see FFmaGenericLoader<br>
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary)
//!     \param  -fdlpoly name  file coming from a DLpoly simulation
//!
//!     \param   -fout name: generic name for files (without extension) and save data
//!                  with following format in name.fma or name.bfma if -bin is set"
//!      \param   -bin save output in binary mode (name file  name.bfma
//!
//!
//! \b examples
//!
//!  Transform an ascii file in a binary file
//!
//!   changeFormat -fin unitCubeXYZQ100.fma  -fout unitCubeXYZQ100 -bin

//    \param  -visufmt format for the visu file (vtk, vtp, cvs or cosmo). vtp is the default

int main(int argc, char ** argv){

    FHelpDescribeAndExit(argc, argv,
                         "Driver to change the format of the input file. "
                         "fdlpoly is not supported for now.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::OutputVisuFile);
    typedef double FReal;
    const std::string filename(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options,   "data.fma"));

    FFmaGenericLoader<FReal> loader(filename);

    const FSize NbPoints  = loader.getNumberOfParticles();
    const unsigned int nbData   = loader.getNbRecordPerline() ;
    const FSize arraySize =nbData*NbPoints;

    FReal * particles = new FReal[arraySize] ;
    std::memset(particles,0,arraySize*sizeof(FReal));

    FSize j = 0 ;
    for(FSize idxPart = 0 ; idxPart < NbPoints ;++idxPart, j+=nbData){
        loader.fillParticle(&particles[j],nbData);
    }


    /////////////////////////////////////////////////////////////////////////
    //                                           Save data
    /////////////////////////////////////////////////////////////////////////
    //
    //  Generate file for ScalFMM FMAGenericLoader
    //
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
        std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
          FFmaGenericWriter<FReal> writer(name) ;
        writer.writeHeader( loader.getCenterOfBox(), loader.getBoxWidth() , NbPoints, sizeof(FReal), nbData) ;
        writer.writeArrayOfReal(particles, nbData, NbPoints);
    }
    //
    //   Generate file for visualization purpose
    //
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputVisuFile.options)){
        std::string outfilename(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.vtp"));
        driverExportData(outfilename, particles , NbPoints,loader.getNbRecordPerline() );
    }
    //
    delete[] particles ;
    return 0;
}
