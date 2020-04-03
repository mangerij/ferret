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

#include "Containers/FBoolArray.hpp"


/**
* This file is a unit test for the FBoolArray class
*/


/** this class test the bool array container */
class TestArray : public FUTester<TestArray> {

    void TestGetSet(){
        FBoolArray array(500);
        for(int idx = 0 ; idx < 500 ; ++idx){
            uassert(!array.get(idx));
        }

        for(int idx = 0 ; idx < 500 ; ++idx){
            array.set(idx, true);
            uassert(array.get(idx));
            array.set(idx, false);
            uassert(!array.get(idx));
        }

        for(int idx = 0 ; idx < 500 ; ++idx){
            array.set(idx, true);
        }
        array.setToZeros();
        for(int idx = 0 ; idx < 500 ; ++idx){
            uassert(!array.get(idx));
        }
        array.setToOnes();
        for(int idx = 0 ; idx < 500 ; ++idx){
            uassert(array.get(idx));
        }
    }

    void TestGetSet2(){
        FBoolArray array(100);

        for(int idx = 0 ; idx < 100 ; ++idx){
            if(idx%3){
                array.set(idx, true);
                uassert(array.get(idx));
            }
            else{
                uassert(!array.get(idx));
            }
        }
    }

    void TestEqual(){
        FBoolArray array1(10);
        FBoolArray array2(10);


        uassert(array1 == array2);

        array1.set(1, true);
        uassert(array1 != array2);

        array2.set(1, true);
        uassert(array1 == array2);

        array1.set(5, true);
        array2 = array1;
        uassert(array1 == array2);
    }

    // set test
    void SetTests(){
        AddTest(&TestArray::TestGetSet,"Test Get & Set");
        AddTest(&TestArray::TestGetSet2,"Test Get & Set 2");
        AddTest(&TestArray::TestEqual,"Test Equal");
    }
};

// You must do this
TestClass(TestArray)


