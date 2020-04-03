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
#include "Containers/FNeighborIndexes.hpp"

/**
* This file is a unit test for the FNeigborIndexes classes
*/


/** this class test the list container */
class TestIndexes : public FUTester<TestIndexes> {
    void MortonLimite(){
        for(int idxLevel = 1 ; idxLevel < MaxTreeHeight ; ++idxLevel){
            const int limit = FMath::pow2(idxLevel);
            for(int x = 0 ; x < 3 ; ++x){
                const int xbox = (x*(limit-1))/2;
                for(int y = 0 ; y < 3 ; ++y){
                    const int ybox = (y*(limit-1))/2;
                    for(int z = 0 ; z < 3 ; ++z){
                        const int zbox = (z*(limit-1))/2;

                        const MortonIndex mindex = FTreeCoordinate(xbox, ybox, zbox).getMortonIndex(idxLevel);

                        FCoordinateNeighborIndex coordindexes(mindex, idxLevel);
                        FBitsNeighborIndex bitsindex(mindex, idxLevel);

                        uassert(coordindexes.minX() == bitsindex.minX());
                        uassert(coordindexes.minY() == bitsindex.minY());
                        uassert(coordindexes.minZ() == bitsindex.minZ());

                        uassert(coordindexes.maxX() == bitsindex.maxX());
                        uassert(coordindexes.maxY() == bitsindex.maxY());
                        uassert(coordindexes.maxZ() == bitsindex.maxZ());

                        for(int idxX = coordindexes.minX() ; idxX <= coordindexes.maxX() ; ++idxX){
                            for(int idxY = coordindexes.minY() ; idxY <= coordindexes.maxY() ; ++idxY){
                                for(int idxZ = coordindexes.minZ() ; idxZ <= coordindexes.maxZ() ; ++idxZ){
                                    if(idxX || idxY || idxZ){
                                        uassert(coordindexes.getIndex(idxX, idxY, idxZ)
                                                == bitsindex.getIndex(idxX, idxY, idxZ));
                                        const MortonIndex neigh = bitsindex.getIndex(idxX, idxY, idxZ);
                                        const FTreeCoordinate neighCoord(neigh, idxLevel);
                                        uassert(xbox + idxX == neighCoord.getX());
                                        uassert(ybox + idxY == neighCoord.getY());
                                        uassert(zbox + idxZ == neighCoord.getZ());

                                        uassert(bitsindex.areNeighbors(mindex, neigh));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void Morton(){
        srand(0);
        for(int idxLevel = 1 ; idxLevel < MaxTreeHeight ; ++idxLevel){
            const int limit = FMath::pow2(idxLevel);
            for(int idxTest = 0 ; idxTest < 100 ; ++idxTest){
                const int xbox = int(drand48()*double(limit));
                const int ybox = int(drand48()*double(limit));
                const int zbox = int(drand48()*double(limit));

                const MortonIndex mindex = FTreeCoordinate(xbox, ybox, zbox).getMortonIndex(idxLevel);

                FCoordinateNeighborIndex coordindexes(mindex, idxLevel);
                FBitsNeighborIndex bitsindex(mindex, idxLevel);

                uassert(coordindexes.minX() == bitsindex.minX());
                uassert(coordindexes.minY() == bitsindex.minY());
                uassert(coordindexes.minZ() == bitsindex.minZ());

                uassert(coordindexes.maxX() == bitsindex.maxX());
                uassert(coordindexes.maxY() == bitsindex.maxY());
                uassert(coordindexes.maxZ() == bitsindex.maxZ());

                for(int idxX = coordindexes.minX() ; idxX <= coordindexes.maxX() ; ++idxX){
                    for(int idxY = coordindexes.minY() ; idxY <= coordindexes.maxY() ; ++idxY){
                        for(int idxZ = coordindexes.minZ() ; idxZ <= coordindexes.maxZ() ; ++idxZ){
                            if(idxX || idxY || idxZ){
                                uassert(coordindexes.getIndex(idxX, idxY, idxZ)
                                        == bitsindex.getIndex(idxX, idxY, idxZ));
                                const MortonIndex neigh = bitsindex.getIndex(idxX, idxY, idxZ);
                                const FTreeCoordinate neighCoord(neigh, idxLevel);
                                uassert(xbox + idxX == neighCoord.getX());
                                uassert(ybox + idxY == neighCoord.getY());
                                uassert(zbox + idxZ == neighCoord.getZ());

                                uassert(bitsindex.areNeighbors(mindex, neigh));
                                if(bitsindex.areNeighbors(mindex, neigh) == false){
                                    printf("Stop here\n");// TODO remove
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void MortonAll(){
        {
            const int idxLevel = 5;
            const int limit = FMath::pow2(idxLevel);
            for(int xbox = 0 ; xbox < limit ; ++xbox){
                for(int ybox = 0 ; ybox < limit ; ++ybox){
                    for(int zbox = 0 ; zbox < limit ; ++zbox){

                        const MortonIndex mindex = FTreeCoordinate(xbox, ybox, zbox).getMortonIndex(idxLevel);

                        FCoordinateNeighborIndex coordindexes(mindex, idxLevel);
                        FBitsNeighborIndex bitsindex(mindex, idxLevel);

                        uassert(coordindexes.minX() == bitsindex.minX());
                        uassert(coordindexes.minY() == bitsindex.minY());
                        uassert(coordindexes.minZ() == bitsindex.minZ());

                        uassert(coordindexes.maxX() == bitsindex.maxX());
                        uassert(coordindexes.maxY() == bitsindex.maxY());
                        uassert(coordindexes.maxZ() == bitsindex.maxZ());

                        for(int idxX = coordindexes.minX() ; idxX <= coordindexes.maxX() ; ++idxX){
                            for(int idxY = coordindexes.minY() ; idxY <= coordindexes.maxY() ; ++idxY){
                                for(int idxZ = coordindexes.minZ() ; idxZ <= coordindexes.maxZ() ; ++idxZ){
                                    if(idxX || idxY || idxZ){
                                        uassert(coordindexes.getIndex(idxX, idxY, idxZ)
                                                == bitsindex.getIndex(idxX, idxY, idxZ));
                                        const MortonIndex neigh = bitsindex.getIndex(idxX, idxY, idxZ);
                                        const FTreeCoordinate neighCoord(neigh, idxLevel);
                                        uassert(xbox + idxX == neighCoord.getX());
                                        uassert(ybox + idxY == neighCoord.getY());
                                        uassert(zbox + idxZ == neighCoord.getZ());

                                        uassert(bitsindex.areNeighbors(mindex, neigh));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void Neighbors(){
        {
            const int idxLevel = 5;
            const MortonIndex limit = FMath::pow2(idxLevel)*3L;
            for(MortonIndex idxV1 = 1 ; idxV1 < limit ; ++idxV1){
                FBitsNeighborIndex bitsindex(idxV1, idxLevel);
                const FTreeCoordinate coord1(idxV1, idxLevel);

                for(MortonIndex idxV2 = 0 ; idxV2 < limit ; ++idxV2){
                    const bool isneig = bitsindex.areNeighbors(idxV1, idxV2);

                    const FTreeCoordinate coord2(idxV2, idxLevel);
                    const bool isreallyneig = (FMath::Abs(coord1.getX()-coord2.getX()) <= 1)
                            && (FMath::Abs(coord1.getY()-coord2.getY()) <= 1)
                            && (FMath::Abs(coord1.getZ()-coord2.getZ()) <= 1);

                    uassert(isneig == isreallyneig);
                }
            }
        }
    }


    // set test
    void SetTests(){
        AddTest(&TestIndexes::MortonLimite,"Test NeighborIndexes at limits");
        AddTest(&TestIndexes::Morton,"Test NeighborIndexes");
        AddTest(&TestIndexes::MortonAll,"Test All NeighborIndexes");
        AddTest(&TestIndexes::Neighbors,"Test Neighbors");
    }
};

// You must do this
TestClass(TestIndexes)



