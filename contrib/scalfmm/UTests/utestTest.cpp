// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#include "FUTester.hpp"

// compile by g++ utestTest.cpp -o utestTest.exe

/** this class show a simple example of unit test */
class MyTest : public FUTester<MyTest> {
	void Before(){
		Print("Before running the test");
	}

	void TestOne(){
                uassert(true);
                //or uassert(false); make an error
                uassert(1 == 1);
	}
	
	void TestTwo(){
		equal(1 , 1);
                different(1 , 2);
	}
	
	void After(){
		Print("After running the test");
	}
	
	void PreTest(){
		Print("Before each test");
	}
	
	void PostTest(){
		Print("After each test");
	}
	
	// You must implement it
	void SetTests(){
            AddTest(&MyTest::TestOne);
            AddTest(&MyTest::TestTwo,"My Second Test");
	}
};

// You must do this
TestClass(MyTest)


