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

#include "Containers/FVector.hpp"


/**
* This file is a unit test for the FVector class
*/

/**
* This class is simply used to count alloc dealloc
*/
class TestObject{
public:
	static int counter;
	static int dealloced;

	TestObject(){
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


/** this class test the vector container */
class TestVector : public FUTester<TestVector> {
	// Called before each test : simply set counter to 0
	void PreTest(){
		TestObject::counter = 0;
		TestObject::dealloced = 0;
	}

	// test size
	void TestSize(){
                FVector<TestObject> vector;
                vector.push(TestObject());
                vector.push(TestObject());
                vector.push(TestObject());
                uassert(vector.getSize() == 3);
		
                uassert((TestObject::counter - TestObject::dealloced) == vector.getSize());

                vector.clear();
                uassert(vector.getSize() == 0);

                uassert(TestObject::counter == TestObject::dealloced);
	}
	
	// test copy
	void TestCopy(){
                FVector<TestObject> vector;
                vector.push(TestObject());
                vector.push(TestObject());
                vector.push(TestObject());

                {
                    FVector<TestObject> vector2(vector);
                    uassert(vector.getSize() == vector2.getSize());
                    uassert((TestObject::counter - TestObject::dealloced) == (vector.getSize() + vector2.getSize()));
                }
                {
                    FVector<TestObject> vector2(vector.getSize()/2);
                    vector2 = vector;
                    uassert(vector.getSize() == vector2.getSize());
                    uassert((TestObject::counter - TestObject::dealloced) == (vector.getSize() + vector2.getSize()));
                }
	}

	// test iter
	void TestIter(){		
                FVector<TestObject> vector;
		{
                        FVector<TestObject>::ConstBasicIterator iter(vector);
                        uassert(!iter.hasNotFinished());
		}
		{
                        vector.push(TestObject());
                        vector.push(TestObject());
                        vector.push(TestObject());

                        FVector<TestObject>::ConstBasicIterator iter(vector);
                        uassert(iter.hasNotFinished());

			int counter = 0;
			while(iter.hasNotFinished()){ iter.gotoNext(); ++counter; }
                        uassert(!iter.hasNotFinished());
                        uassert(counter == vector.getSize());
		}
	}

        // test remove
        void TestRemove(){
                FVector<TestObject> vector;
                {
                        FVector<TestObject>::BasicIterator iter(vector);
                        uassert(!iter.hasNotFinished());
                }
                {
                        vector.push(TestObject());
                        vector.push(TestObject());
                        vector.push(TestObject());

                        FVector<TestObject>::BasicIterator iter(vector);
                        uassert(iter.hasNotFinished());

                        iter.gotoNext();
                        iter.remove();
                        uassert(iter.hasNotFinished());

                        iter.gotoNext();
                        uassert(!iter.hasNotFinished());

                        uassert(2 == vector.getSize());
                }
        }
		
	// set test
	void SetTests(){
            AddTest(&TestVector::TestSize,"Test Size");
            AddTest(&TestVector::TestCopy,"Test Copy");
            AddTest(&TestVector::TestIter,"Test Iter");
            AddTest(&TestVector::TestRemove,"Test Remove");
	}
};

// You must do this
TestClass(TestVector)


