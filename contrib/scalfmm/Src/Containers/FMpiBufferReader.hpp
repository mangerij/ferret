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
#ifndef FMPIBUFFERREADER_HPP
#define FMPIBUFFERREADER_HPP

#include <memory>
#include "../Utils/FMpi.hpp"
#include "FAbstractBuffer.hpp"
#include "../Utils/FAssert.hpp"

/** @author Cyrille Piacibello
 * This class provide the same features as FBufferWriter using MPI_Pack system
 *
 * Put some data
 * then insert back if needed
 * finally use data pointer as you like
 */
class FMpiBufferReader : public FAbstractBufferReader {
    MPI_Comm comm;            //< Communicator needed by MPI_Pack functions
    FSize arrayCapacity;        //< Allocated space
    std::unique_ptr<char[]> array;  //< Allocated Array
    FSize currentIndex;

public :
    /*Constructor with a default arrayCapacity of 512 bytes */
    explicit FMpiBufferReader(const MPI_Comm inComm = MPI_COMM_WORLD, const FSize inDefaultCapacity = 512):
        comm(inComm),
        arrayCapacity(inDefaultCapacity),
        array(new char[inDefaultCapacity]),
        currentIndex(0){
        FAssertLF(array, "Cannot allocate array");
    }

    /** Change the comm (or to set it later) */
    void setComm(const MPI_Comm inComm){
        comm = inComm;
    }

    /** To change the capacity (but reset the head to 0) */
    void cleanAndResize(const FSize newCapacity){
        if(newCapacity != arrayCapacity){
            arrayCapacity = newCapacity;
            array.reset(new char[newCapacity]);
        }
        currentIndex = 0;
    }

    /** Destructor
   */
    virtual ~FMpiBufferReader(){
    }

    /** Get allocated memory pointer */
    char* data() override {
        return array.get();
    }

    /** Get allocated memory pointer */
    const char* data() const override  {
        return array.get();
    }

    /** get the filled space */
    FSize getSize() const override {
        return currentIndex;
    }

    /** Size of the memory initialized */
    FSize getCapacity() const{
        return arrayCapacity;
    }

    /** Move the read index to a position */
    void seek(const FSize inIndex) override {
        FAssertLF(inIndex <= arrayCapacity, "FMpiBufferReader :: Aborting :: Can't move index because buffer isn't long enough ",inIndex," ",arrayCapacity);
        currentIndex = inIndex;
    }

    /** Get the read position */
    FSize tell() const override  {
        return currentIndex;
    }

    /** Get a value with memory cast */
    template <class ClassType>
    ClassType getValue(){
        FAssertLF(arrayCapacity < std::numeric_limits<int>::max());
        FAssertLF(currentIndex < std::numeric_limits<int>::max());
        int previousIndex = int(currentIndex);
        ClassType value;
        FMpi::Assert(MPI_Unpack(array.get(),int(arrayCapacity),&previousIndex,&value,1,FMpi::GetType(value),comm), __LINE__);
        seek(FSize(sizeof(value) + currentIndex));
        FAssertLF(previousIndex == currentIndex);
        return value;
    }

    /** Get a value with memory cast at a specified index */
    template <class ClassType>
    ClassType getValue(const FSize ind){
        ClassType value;
        FAssertLF(arrayCapacity < std::numeric_limits<int>::max());
        FAssertLF(ind < std::numeric_limits<int>::max());
        int previousIndex = int(ind);
        FMpi::Assert(MPI_Unpack(array.get(),int(arrayCapacity),&previousIndex,&value,1,FMpi::GetType(value),comm), __LINE__);
        seek(FSize(sizeof(value)+ind));
        FAssertLF(previousIndex == currentIndex);
        return value;
    }

    /** Fill a value with memory cast */
    template <class ClassType>
    void fillValue(ClassType* const inValue){
        FAssertLF(arrayCapacity < std::numeric_limits<int>::max());
        FAssertLF(currentIndex < std::numeric_limits<int>::max());
        int previousIndex = int(currentIndex);
        FMpi::Assert(MPI_Unpack(array.get(),int(arrayCapacity),&previousIndex,inValue,1,FMpi::GetType(*inValue),comm), __LINE__);
        seek(FSize(sizeof(ClassType) + currentIndex));
        FAssertLF(previousIndex == currentIndex);
    }

    /** Fill one/many value(s) with memcpy */
    template <class ClassType>
    void fillArray(ClassType* const inArray, const FSize inSize){
        FAssertLF(arrayCapacity < std::numeric_limits<int>::max());
        FAssertLF(currentIndex < std::numeric_limits<int>::max());
        FAssertLF(inSize < std::numeric_limits<int>::max());
        int previousIndex = int(currentIndex);
        FMpi::Assert(MPI_Unpack(array.get(),int(arrayCapacity),&previousIndex,inArray,int(inSize),FMpi::GetType(*inArray),comm), __LINE__);
        seek(FSize(sizeof(ClassType) * inSize + currentIndex));
        FAssertLF(previousIndex == currentIndex);
    }

    /** Same as fillValue */
    template <class ClassType>
    FMpiBufferReader& operator>>(ClassType& object){
        fillValue(&object);
        return *this;
    }

};
#endif

