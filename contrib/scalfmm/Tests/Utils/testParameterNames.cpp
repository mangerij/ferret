#include "../../Src/Utils/FParameterNames.hpp"

#include <iostream>

int main(int argc, char** argv){

    //FHelpAndExit(argc, argv, FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeSubHeight );
    FHelpDescribeAndExit(argc, argv, "It is made to test the parameters common to all binaries.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeSubHeight);

    std::cout << "Pass -help to the executable to know how to use it\n";

    return 0;
}
