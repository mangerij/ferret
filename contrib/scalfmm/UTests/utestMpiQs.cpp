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
#ifndef UTESTMPIQS_HPP
#define UTESTMPIQS_HPP

#include "Utils/FGlobal.hpp"
#include "FUTester.hpp"

#include "Utils/FMpi.hpp"
#include "Utils/FQuickSortMpi.hpp"

#include <memory>
#include <limits>
#include <cstdlib>

// ==== CMAKE =====
// @FUSE_MPI
// ================

/** this class test the mpi quick sort */
class TestMpiQs : public FUTesterMpi<TestMpiQs> {
    ////////////////////////////////////////////////////////////
    /// Check function
    ////////////////////////////////////////////////////////////

    /** To test if an array is sorted */
    template <class ValueType, class IndexType>
    void CheckIfSorted(const ValueType array[], const IndexType size){
        for(int idx = 1 ; idx < size ; ++idx){
            uassert(array[idx-1] <= array[idx]);
        }
    }
    /** To test if the global sorting is correct */
    template <class ValueType, class IndexType>
    void CheckBorder(const ValueType array[], const IndexType size){
        const int myRank = app.global().processId();
        const int nbProcess = app.global().processCount();

        ValueType atMyLeft;     // The value on my left [myRank-1][size-1]
        ValueType atMyRight;    // The value on my right[myRank+1][0]
        ValueType myLeft;       // My left value [myRank][0]
        ValueType myRight;      // My right value [myRank][size-1]

        // If there is someone on my left
        if(myRank != 0){
            // If I do not have value I shoud replace it
            if(size == 0){
                // If there is no-one on my right I take the max
                if(myRank == nbProcess-1){
                    myLeft = std::numeric_limits<ValueType>::max();
                }
                // [A] Else I receive the value from the right and say it is my left
                else{
                    FMpi::Assert(MPI_Recv(&atMyRight, sizeof(atMyRight), MPI_BYTE, myRank+1, 0, app.global().getComm(), MPI_STATUS_IGNORE) , __LINE__);
                    myLeft = atMyRight;
                }
            }
            else{
                // Take my left value
                myLeft = array[0];
            }
            // Send it to my left neighbors
            FMpi::Assert(MPI_Send(&myLeft, sizeof(myLeft), MPI_BYTE , myRank-1, 0, app.global().getComm()) , __LINE__);
        }
        // If there is someone on my right
        if(myRank != nbProcess-1){
            // I should receive the value (if not already done in [A])
            if(size != 0){
                FMpi::Assert(MPI_Recv(&atMyRight, sizeof(atMyRight), MPI_BYTE, myRank+1, 0, app.global().getComm(), MPI_STATUS_IGNORE) , __LINE__);
            }
            // If I do not have value I shoud replace it
            if(size == 0){
                // If there is no-one on my left I take the min
                if(myRank == 0){
                    myRight = std::numeric_limits<ValueType>::min();
                }
                // [B] Else I receive the value from the left and say it is my right
                else{
                    FMpi::Assert(MPI_Recv(&atMyLeft, sizeof(atMyLeft), MPI_BYTE, myRank-1, 0, app.global().getComm(), MPI_STATUS_IGNORE) , __LINE__);
                    myRight = atMyLeft;
                }
            }
            else{
                // Take my right value
                myRight = array[size-1];
            }
            // Send it to my right neighbors
            FMpi::Assert(MPI_Send(&myRight, sizeof(myRight), MPI_BYTE , myRank+1, 0, app.global().getComm()) , __LINE__);
        }
        // If there is someone on my left
        if(myRank != 0){
            // If not already receive in [B]
            if(size != 0){
                FMpi::Assert(MPI_Recv(&atMyLeft, sizeof(atMyLeft), MPI_BYTE, myRank-1, 0, app.global().getComm(), MPI_STATUS_IGNORE) , __LINE__);
            }
        }
        // Test only if I hold data and if someone on my left
        if(myRank != 0 && size != 0){
            uassert(atMyLeft <= myLeft);
        }
        // Test only if I hold data and if someone on my right
        if(myRank != nbProcess-1 && size != 0){
            uassert(myRight <= atMyRight);
        }
    }

    ////////////////////////////////////////////////////////////
    /// The tests
    ////////////////////////////////////////////////////////////

    void TestSmallSort(){
        //const int myRank = app.global().processId();
        const int nbProcess = app.global().processCount();

        const int nbElements = nbProcess;
        std::unique_ptr<long[]> elements(new long[nbElements]);

        for(int idx = 0 ; idx < nbElements ; ++idx){
            elements[idx] = idx;
        }

        const int nbElementsInTest = app.global().reduceSum(nbElements);

        int outSize;
        long* outElements = nullptr;
        FQuickSortMpi<long, long, int>::QsMpi(elements.get(), nbElements, &outElements, &outSize, app.global());

        CheckIfSorted(outElements, outSize);
        CheckBorder(outElements, outSize);

        uassert(nbElementsInTest == app.global().reduceSum(outSize));
        delete[] outElements;
    }

    void TestTinySort(){
        const int myRank = app.global().processId();
        const int nbProcess = app.global().processCount();

        const int nbElements = (myRank == 0 ? nbProcess : 0);
        std::unique_ptr<long[]> elements(new long[nbElements]);

        for(int idx = 0 ; idx < nbElements ; ++idx){
            elements[idx] = idx;
        }

        const int nbElementsInTest = app.global().reduceSum(nbElements);

        int outSize;
        long* outElements = nullptr;
        FQuickSortMpi<long, long, int>::QsMpi(elements.get(), nbElements, &outElements, &outSize, app.global());

        CheckIfSorted(outElements, outSize);
        CheckBorder(outElements, outSize);

        uassert(nbElementsInTest == app.global().reduceSum(outSize));
        delete[] outElements;
    }

    void TestSameSort(){
        //const int myRank = app.global().processId();
        const int nbProcess = app.global().processCount();

        const int nbElements = nbProcess * 100;
        std::unique_ptr<long[]> elements(new long[nbElements]);

        for(int idx = 0 ; idx < nbElements ; ++idx){
            elements[idx] = nbProcess;
        }

        const int nbElementsInTest = app.global().reduceSum(nbElements);

        int outSize;
        long* outElements = nullptr;
        FQuickSortMpi<long, long, int>::QsMpi(elements.get(), nbElements, &outElements, &outSize, app.global());

        CheckIfSorted(outElements, outSize);
        CheckBorder(outElements, outSize);

        uassert(nbElementsInTest == app.global().reduceSum(outSize));
        delete[] outElements;
    }


    void TestUniqueSort(){
        const int myRank = app.global().processId();
        const int nbProcess = app.global().processCount();

        const int nbElements = nbProcess * 100;
        std::unique_ptr<long[]> elements(new long[nbElements]);

        for(int idx = 0 ; idx < nbElements ; ++idx){
            elements[idx] = myRank;
        }

        const int nbElementsInTest = app.global().reduceSum(nbElements);

        int outSize;
        long* outElements = nullptr;
        FQuickSortMpi<long, long, int>::QsMpi(elements.get(), nbElements, &outElements, &outSize, app.global());

        CheckIfSorted(outElements, outSize);
        CheckBorder(outElements, outSize);

        uassert(nbElementsInTest == app.global().reduceSum(outSize));
        delete[] outElements;
    }

    void TestBigSort(){
        //const int myRank = app.global().processId();
        //const int nbProcess = app.global().processCount();

        const int nbElements = 500;
        std::unique_ptr<long[]> elements(new long[nbElements]);

        for(int idx = 0 ; idx < nbElements ; ++idx){
            elements[idx] = long(drand48() * 10000);
        }

        const int nbElementsInTest = app.global().reduceSum(nbElements);

        int outSize;
        long* outElements = nullptr;
        FQuickSortMpi<long, long, int>::QsMpi(elements.get(), nbElements, &outElements, &outSize, app.global());

        CheckIfSorted(outElements, outSize);
        CheckBorder(outElements, outSize);

        uassert(nbElementsInTest == app.global().reduceSum(outSize));
        delete[] outElements;
    }



    void TestPivotSort(){
        //const int myRank = app.global().processId();
        //const int nbProcess = app.global().processCount();

        const int nbElements = 500;
        std::unique_ptr<long[]> elements(new long[nbElements]);

        for(int idx = 0 ; idx < nbElements ; ++idx){
            elements[idx] = 10000;
        }
        elements[0] = 9999;

        const int nbElementsInTest = app.global().reduceSum(nbElements);

        int outSize;
        long* outElements = nullptr;
        FQuickSortMpi<long, long, int>::QsMpi(elements.get(), nbElements, &outElements, &outSize, app.global());

        CheckIfSorted(outElements, outSize);
        CheckBorder(outElements, outSize);

        uassert(nbElementsInTest == app.global().reduceSum(outSize));
        delete[] outElements;
    }

    ////////////////////////////////////////////////////////////
    /// Starter functions
    ////////////////////////////////////////////////////////////

    // set test
    void SetTests(){
        AddTest(&TestMpiQs::TestTinySort,"Test tiny sort with values only on p0 at starting");
        AddTest(&TestMpiQs::TestSmallSort,"Test small sort");
        AddTest(&TestMpiQs::TestSameSort,"Test with all the same value");
        AddTest(&TestMpiQs::TestUniqueSort,"Test with unique value");
        AddTest(&TestMpiQs::TestBigSort,"Test with random values");
        AddTest(&TestMpiQs::TestPivotSort,"Test with big pivot");
    }
public:
    TestMpiQs(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }
};

// You must do this
TestClassMpi(TestMpiQs)

#endif // UTESTMPIQS_HPP
