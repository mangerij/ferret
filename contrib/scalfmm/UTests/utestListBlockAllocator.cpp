
// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Brenger Bramas, Matthias Messner
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

#include "Containers/FBlockAllocator.hpp"

#include <cstring>

/**
* This file is a unit test for the List block allocator
*/

/**
* This class is simply used to count alloc dealloc
*/
static const int SizeArray = 50;
class TestObject{
public:
    static int counter;
    static int dealloced;

    int array[SizeArray];

    TestObject(){
        memset(array, 0, SizeArray * sizeof(int));
        ++counter;
    }
    TestObject(const TestObject&){
        ++counter;
    }
    ~TestObject(){
        ++dealloced;
    }
};

int TestObject::counter(0);
int TestObject::dealloced(0);


/** this class test the list container */
class TestBlock : public FUTester<TestBlock> {
    // Called before each test : simply set counter to 0
    void PreTest(){
        TestObject::counter = 0;
        TestObject::dealloced = 0;
    }

    // test copy
    void TestBlockFunction(){
        FListBlockAllocator<TestObject, 10> alloc;

        const int NbAlloc = 2;
        TestObject* ptr[NbAlloc];
        for(int idx = 0 ; idx < NbAlloc ; ++idx){
            TestObject* dl1 = alloc.newObject();
            TestObject* dl2 = alloc.newObject();
            alloc.deleteObject(dl1);
            ptr[idx] = alloc.newObject();
            alloc.deleteObject(dl2);
        }

        // Just put on volatile to avoid loop unrolling (which leads to a bug with my gcc!)
        volatile int volSizeArray = SizeArray;
        volatile int volNbAlloc = NbAlloc;

        for(int idx = 0 ; idx < volNbAlloc ; ++idx){
            for(int idxval = 0 ; idxval < volSizeArray ; ++idxval){
                ptr[idx]->array[idxval] += (idxval * idx);
            }
        }

        for(int idx = 0 ; idx < NbAlloc ; ++idx){
            for(int idxval = 0 ; idxval < SizeArray ; ++idxval){
                uassert(ptr[idx]->array[idxval] == (idxval * idx));
            }
        }

        for(int idx = 0 ; idx < NbAlloc ; ++idx){
            alloc.deleteObject(ptr[idx]);
        }

        uassert(TestObject::counter == (3*NbAlloc));
        uassert(TestObject::counter == TestObject::dealloced);
    }

    // set test
    void SetTests(){
            AddTest(&TestBlock::TestBlockFunction,"Test Allocate Deallocate");
    }
};

// You must do this
TestClass(TestBlock)


