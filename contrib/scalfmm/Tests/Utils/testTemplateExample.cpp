// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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


/** This file gives an example of the template helper usage.
  *
  **/

#include "../../Src/Utils/FTemplate.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include <iostream>

struct RunContainer{
    template <const int P>
    static void For(double mydouble, int myint){
        std::cout << "For >> P = " << P << " mydouble " << mydouble << " myint " << myint << "\n";
    }

    template <const int P>
    static void Run(double mydouble, int myint){
        std::cout << "Run >> P = " << P << " mydouble " << mydouble << " myint " << myint << "\n";
    }

    template <const int P>
    static void Run(){
        std::cout << "Run without param >> P = " << P << "\n";
    }
};

struct RunClass {
    template <const int P>
    void For(double mydouble, int myint){
        std::cout << "RunClass::For >> P = " << P << " mydouble " << mydouble << " myint " << myint << "\n";
    }

    template <const int P>
    void Run(double mydouble, int myint){
        std::cout << "RunClass::Run >> P = " << P << " mydouble " << mydouble << " myint " << myint << "\n";
    }
};


// Compile with g++ -std=c++11 main.cpp -o test.exe
int main(int argc, char** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Only the code is interesting in order to understand the multi-template compilation.");

    // This will call the For method (maybe several times)
    FForAll::For<int, 0, 20, 1, RunContainer>(45.4, 55);

    // This will call the For method (maybe several times)
    RunClass runner;
    FForAllThis::For<int, 0, 20, 1, RunClass>(&runner, 45.4, 55);

    // This will call the Run method
    FRunIf::Run<int, 0, 20, 1, RunContainer>(3, 45.4, 55);

    // This will call the Run method
    FRunIfThis::Run<int, 0, 20, 1, RunClass>(&runner, 3, 45.4, 55);

    // Should not run, because 21 is out of the compiled interval
    FRunIf::Run<int, 0, 20, 1, RunContainer>(21, 45.4, 55);

    // Call the Run method with no param
    FRunIf::Run<int, 0, 20, 1, RunContainer>(18);

    // Call the Run method with no param for odd P
    FRunIfFunctional::Run<int, 0, 20, 1, RunContainer>([&](int p){ return p&1;});

    // Call the Run method with params for even P
    FRunIfFunctional::Run<int, 0, 20, 1, RunContainer>([&](int p){ return !(p&1);}, 45.4, 55);

    return 0;
}
