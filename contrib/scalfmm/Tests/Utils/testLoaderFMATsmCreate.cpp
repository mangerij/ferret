// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, B��renger Bramas, Matthias Messner
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
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

// This file can generate basic particles files in the FMA format

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Create a file for the TSM (target source model)",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::PhysicalValue);

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable can create a FMA-like particles files (this is not a fma file!)";
    std::cout << ">> You can pass a filename in parameter else the program will use\n";
    std::cout << ">> a default filename.\n";
    std::cout << ">> The format of the file is : \n";
    std::cout << ">> [number of particles] \n";
    std::cout << ">> [boxe width] [boxe x center] [boxe y center] [boxe z center]\n";
    std::cout << ">> [x] [y] [z] [physical value] [1 if target 0 if source]...\n";
    //////////////////////////////////////////////////////////////
    // Nb of particles
    typedef double FReal;
    const long NbParticles = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, long(20000));
    const FReal physicalValue = FParameters::getValue(argc,argv,FParameterDefinitions::PhysicalValue.options, FReal(0.1));

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char* const Output = FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "../Data/test20k.tsm.fma");
    std::cout << "Creating : " << Output << "\n";

    // Create file
    std::ofstream myfile;
    myfile.open (Output);

    if(!myfile.is_open()){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }

    std::cout << "Generating " << NbParticles << " in " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticles << "\n";
    myfile << BoxWidth << "\t" << XCenter << "\t" << YCenter << "\t" << ZCenter;
	srand48( static_cast<long>(time(nullptr))) ;

    // Generate particles
    for( long idx = 0 ; idx < NbParticles ; ++idx ){
        const FReal px = ((FReal(drand48())) * BoxWidth * FReal(2)) + XCenter - BoxWidth;
        const FReal py = ((FReal(drand48())) * BoxWidth * FReal(2)) + YCenter - BoxWidth;
        const FReal pz = ((FReal(drand48())) * BoxWidth * FReal(2)) + ZCenter - BoxWidth;

        const int isTarget = drand48() > 0.5 ? 1 : 0;

        myfile << "\n" << px << "\t" << py << "\t" <<  pz << "\t" << physicalValue << "\t" << isTarget;
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}



