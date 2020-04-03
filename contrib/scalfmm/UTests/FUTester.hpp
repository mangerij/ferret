// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#ifndef UTESTER_HPP
#define UTESTER_HPP

#include "ScalFmmConfig.h"

#include <iostream>
#include <list>
#include <string>
#include <cstdio>


#define TestClass(X)\
    int main(void){\
    X Controller;\
    return Controller.Run();\
    }\


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* Please read the license
*
* This class is used to make simple unit test cases
*
* Please refer to testUTest.cpp to see an example
* @warning Create a derived class that implement SetTests() and use TestClass() macro
*
* We recommend to have a look to a unit test to better understand how it works,
* as for example @example TestList
*/
template <class TestClass>
class FUTester{
    // Test function pointer
    typedef void (TestClass::*TestFunc)(void);

    /** Test descriptor */
    struct TestFuncDescriptor{
        TestFunc func;		//< Test adress
        std::string name;	//< Test name
    };


    std::list<TestFuncDescriptor> tests;	//< all tests

    int totalTests;				//< number of tests

    int currentTest;			//< current processing test in the run
    int currentStep;			//< current processing step in the run

    int failedSteps;			//< number of failed step in the current test
    int failedTests;			//< number of failed tests

protected:
    /** Constructor */
    FUTester(){
        totalTests = 0;
    }

    /** Callback before processing test */
    virtual void Before(){}

    /** Callback after processing test */
    virtual void After(){}

    /** Callback before each unit test */
    virtual void PreTest(){}

    /** Callback after each unit test */
    virtual void PostTest(){}

    /**
    * This function has to add tests
        * <code> AddTest(&MyTest::TestOne); </code>
    */
    virtual void SetTests() = 0;

    /**
    * Add a test without giving a name
    * @param inFunc test function address
    */
    void AddTest(TestFunc inFunc){
        char buff[256];
        sprintf(buff,"Unnamed Test number %d",totalTests+1);
        AddTest(inFunc,buff);
    }

    /**
    * Add a test with a name
    * @param inFunc test function address
    * @param inFuncName function name
    */
    void AddTest(TestFunc inFunc, const std::string& inFuncName){
        ++totalTests;
        TestFuncDescriptor desc;
        desc.func = inFunc;
        desc.name = inFuncName;
        tests.push_back(desc);
    }

    /**
    * To print a message manually in the test
    * @param value a object that ostream can work on
    */
    template <class Output>
    void Print(const Output& value){
        std::cout<< "--- Output from program : " << value << "\n";
    }

    /**
    * To test
    * @param result the test result
    * if result is false test failed
    */
    void uassert(const bool result){
        ++currentStep;
        if(!result){
            std::cout << ">> Step " << currentStep << " Failed\n";
            ++failedSteps;
        }
    }

    /**
    * To test equality
    * @param v1 value one
    * @param v2 value 2
    * if v1 is not equal v2 test failed
    */
    template <class T>
    void equal(const T& v1, const T& v2){
        uassert(v1 == v2);
    }

    /**
    * To test equality
    * @param v1 value one
    * @param v2 value 2
    * if v1 is equal v2 test failed
    */
    template <class T>
    void different(const T& v1, const T& v2){
        uassert(v1 != v2);
    }

public :
    /**
    * Processing the test
        * return application exit code (= nb of errors)
    */
    int Run(){
        tests.clear();
        // register tests
        SetTests();

        TestClass* const toTest = static_cast<TestClass*>(this);
        currentTest = 0;
        failedTests = 0;

        Before();

        // for each tests
        const typename std::list<TestFuncDescriptor>::const_iterator end = tests.end();
        for(typename std::list<TestFuncDescriptor>::iterator iter = tests.begin() ; iter != end ; ++iter){
            currentStep = 0;
            failedSteps = 0;

            std::cout << "[Start] " << (*iter).name << "\n";

            PreTest();
            TestFunc ff = (*iter).func;
            (toTest->*ff)();
            PostTest();

            if(failedSteps){
                std::cout << "[Finished] FAILED (" << failedSteps << "/" << currentStep<< " steps failed)\n";
                ++failedTests;
            }
            else{
                std::cout << "[Finished] PASSED (" << currentStep << " steps)\n";
            }

            ++currentTest;
        }


        After();

        std::cout <<"Test is over, " << (totalTests-failedTests) << " Passed, " << failedTests << " Failed\n";

        return failedTests;
    }

};

#ifdef SCALFMM_USE_MPI

#include "Utils/FMpi.hpp"

#define TestClassMpi(X)						\
    int main(int argc, char** argv){				\
    X Controller(argc,argv);					\
    return Controller.Run();					\
    }								\

template <class TestClass>
class FUTesterMpi : public FUTester<TestClass>{
protected:
    FMpi app;

    //Constructor with params to initialize FMpi
    FUTesterMpi(int argc, char ** argv) : app(argc,argv){
    }

    /**
   * To print a message manually in the test
   * @param value a object that ostream can work on
   */
    template <class Output>
    void Print(const Output& value){
        if(app.global().processId()==0){
            std::cout<< "--- Output from program : " << value << "\n";
        }
    }


};

#endif  
#endif
