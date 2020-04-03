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

#include "../Src/Utils/FAlignedMemory.hpp"
#include "../Src/Utils/FTemplate.hpp"

/**
* This file is a unit test for the alignement class
*/



template <std::size_t Value>
struct IncClass{
    static const std::size_t NextValue = (Value<<1);
};

struct CheckClass {
    static int nbObjects;
    int idObject;
    int flag;
    CheckClass(){
        idObject = nbObjects;
        nbObjects += 1;
        flag = -1;
    }
    ~CheckClass(){
        nbObjects -= 1;
    }
};
int CheckClass::nbObjects = 0;


/** this class test the bool array container */
class TestAlignement : public FUTester<TestAlignement> {

    void TestSimple(){
        for(std::size_t idxSize = 1 ; idxSize < 1000 ; idxSize *= 10){
            {
                void* ptr = FAlignedMemory::AllocateBytes<2>(idxSize);
                uassert((reinterpret_cast<std::size_t>(ptr) & (2-1)) == 0);
                FAlignedMemory::DeallocBytes(ptr);
            }
            {
                void* ptr = FAlignedMemory::AllocateBytes<8>(idxSize);
                uassert((reinterpret_cast<std::size_t>(ptr) & (8-1)) == 0);
                FAlignedMemory::DeallocBytes(ptr);
            }
            {
                void* ptr = FAlignedMemory::AllocateBytes<16>(idxSize);
                uassert((reinterpret_cast<std::size_t>(ptr) & (16-1)) == 0);
                FAlignedMemory::DeallocBytes(ptr);
            }
            {
                void* ptr = FAlignedMemory::AllocateBytes<32>(idxSize);
                uassert((reinterpret_cast<std::size_t>(ptr) & (22-1)) == 0);
                FAlignedMemory::DeallocBytes(ptr);
            }
            {
                void* ptr = FAlignedMemory::AllocateBytes<64>(idxSize);
                uassert((reinterpret_cast<std::size_t>(ptr) & (64-1)) == 0);
                FAlignedMemory::DeallocBytes(ptr);
            }
        }
    }

public:
    template <std::size_t Allign>
    void For(){
        for(std::size_t idxSize = 1 ; idxSize < 1000 ; idxSize *= 10){
            void* ptr = FAlignedMemory::AllocateBytes<Allign>(idxSize);
            //uassert((reinterpret_cast<std::size_t>(ptr) & (Allign-1)) == 0);
            FAlignedMemory::DeallocBytes(ptr);
        }
    }

    void TestAll(){
        FForAllThisWithInc::For<std::size_t, 1, 128+1, IncClass, TestAlignement>(this);
    }


    void TestArray(){
        CheckClass::nbObjects = 0;

        const int nb = 5;

        CheckClass* array = FAlignedMemory::AllocateArray<8, CheckClass>( nb );
        uassert(CheckClass::nbObjects == nb);
        for(int idx = 0 ; idx < nb ; ++idx){
            uassert(array[idx].flag == -1);
            uassert(array[idx].idObject == idx);
        }

        FAlignedMemory::DeallocArray( array );
        uassert(CheckClass::nbObjects == 0);
    }


    // set test
    void SetTests(){
        AddTest(&TestAlignement::TestSimple,"Test if aligned");
        AddTest(&TestAlignement::TestAll,"Test if all aligned");
        AddTest(&TestAlignement::TestArray,"Test if array are correct");
    }
};

// You must do this
TestClass(TestAlignement)


