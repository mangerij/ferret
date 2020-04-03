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
#include <iostream>
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameterNames.hpp"

#include <cstdlib>
#include <unistd.h>

/**
* Here we show an example of using FTic
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Only the code is interesting in order to understand the use of timers.");
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use FTic time counter.\n";
    //////////////////////////////////////////////////////////////
    {
	FTic counter;	
	counter.tic();
	usleep(1500000);
	//Sleep(1500); //on windows
	counter.tac();
	std::cout << counter.elapsed() << " (s)\n";
    }
    {
        FTic counter;
        usleep(1500000);
        //Sleep(1500); //on windows
        std::cout << counter.tacAndElapsed() << " (s)\n";
    }
    {
        FTic counter;
        usleep(1500000);
        //Sleep(1500); //on windows
        counter.tac();
        counter.tic();
        usleep(1500000);
        //Sleep(1500); //on windows
        std::cout << counter.tacAndElapsed() << " (s)\n";
        std::cout << counter.cumulated() << " (s)\n";
    }
    return 0;
}

