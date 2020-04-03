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
#include "FUTester.hpp"

#include "Containers/FTreeCoordinate.hpp"

// compile by g++ utestMorton.cpp -o utestMorton.exe

/**
* This file is a unit test for the FTreeCoordinate class
*/


/** this class test the list container */
class TestMorton : public FUTester<TestMorton> {
        void Morton(){
            {
                FTreeCoordinate pos(5,1,7);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                uassert(pos == cp);
                uassert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
            {
                FTreeCoordinate pos(2,8,3);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                uassert(pos == cp);
                uassert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
            {
                FTreeCoordinate pos(51,11,47);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                uassert(pos == cp);
                uassert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
	}

        void Position(){
            {
                FTreeCoordinate pos(0,0,0);
                uassert(pos.getMortonIndex(1) == 0);
            }
            {
                FTreeCoordinate pos(1,1,1);
                uassert(pos.getMortonIndex(1) == 7);
            }
            {
                FTreeCoordinate pos(0,1,1);
                uassert(pos.getMortonIndex(1) == 3);
            }
            {
                FTreeCoordinate pos(2,2,2);
                uassert(pos.getMortonIndex(2) == (7 << 3) );
            }
            {
                FTreeCoordinate pos(1,2,4);
                uassert(pos.getMortonIndex(3) == 84 );// 001 010 100 =>> 001010100 => 84d
            }
	}
		
	// set test
	void SetTests(){
            AddTest(&TestMorton::Morton,"Test Morton");
            AddTest(&TestMorton::Position,"Test Position");
	}
};

// You must do this
TestClass(TestMorton)


