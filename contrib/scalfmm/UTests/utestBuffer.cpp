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
#
#include <cstddef>

#include "FUTester.hpp"

#include "Containers/FBufferReader.hpp"
#include "Containers/FBufferWriter.hpp"



/** this class test the buffers container */
class TestBuffer : public FUTester<TestBuffer> {

        // test size
        void TestWriteRead(){
            FBufferWriter writer;

            const int BytesTested = static_cast<int>(sizeof(int)+sizeof(char)+sizeof(double)+sizeof(float));
            const int NbTest = 5;
            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                writer << idxWrite << char(idxWrite) << double(idxWrite) << float(idxWrite);
            }

            uassert(writer.getSize() == (NbTest*BytesTested));

            FBufferReader reader(writer.getSize());
            uassert(reader.getSize() == writer.getSize());

            memcpy(reader.data(), writer.data(), writer.getSize());
            for(int idxRead = 0 ; idxRead < NbTest ; ++idxRead){
                int intval;
                char charval;
                double doubleval;
                float floatval;
                reader >> intval >> charval >> doubleval >> floatval;

                uassert(intval == idxRead);
                uassert(charval == char(idxRead));
                uassert(doubleval == double(idxRead));
                uassert(floatval == float(idxRead));

                uassert(reader.tell() == (BytesTested * (idxRead+1)));
            }

            uassert(reader.tell() == reader.getSize());
            reader.seek(0);
            uassert(reader.tell() == 0);

            for(int idxRead = 0 ; idxRead < NbTest ; ++idxRead){
                uassert(reader.FBufferReader::getValue<int>() == idxRead);
                uassert(reader.FBufferReader::getValue<char>() == char(idxRead));
                uassert(reader.FBufferReader::getValue<double>() == double(idxRead));
                uassert(reader.FBufferReader::getValue<float>() == float(idxRead));

                uassert(reader.tell() == (BytesTested * (idxRead+1)));
            }

            uassert(reader.tell() == reader.getSize());
        }


        void TestWriteAt(){
            FBufferWriter writer;

            const int SizeOfInt = int(sizeof(int));
            const int NbTest = 5;
            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                const FSize position = writer.getSize();

                uassert(position == (NbTest * SizeOfInt * idxWrite) + (idxWrite * SizeOfInt));

                writer.FBufferWriter::write<int>(0);

                uassert(writer.getSize() == (NbTest * SizeOfInt * idxWrite) + (idxWrite * SizeOfInt) + SizeOfInt);

                for(int count = 1 ; count <= NbTest ; ++count) writer << count;
                writer.writeAt(position, idxWrite);
            }

            FBufferReader reader(writer.getSize());
            memcpy(reader.data(), writer.data(), writer.getSize());

            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                uassert(reader.FBufferReader::getValue<int>() == idxWrite);
                for(int count = 1 ; count <= NbTest ; ++count){
                    uassert(reader.FBufferReader::getValue<int>() == count);
                }
            }
        }


        // set test
        void SetTests(){
            AddTest(&TestBuffer::TestWriteRead,"Test Write then Read");
            AddTest(&TestBuffer::TestWriteAt,"Test Write at then Read");
        }
};

// You must do this
TestClass(TestBuffer)
