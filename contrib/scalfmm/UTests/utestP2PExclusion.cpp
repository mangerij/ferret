
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

#include "Core/FP2PExclusion.hpp"
#include "Utils/FMath.hpp"

#include <memory>

/**
* This file is a unit test for the FNeigborIndexes classes
*/


/** this class test the list container */
class TestExclusion : public FUTester<TestExclusion> {
    const int Size = 100;

    void Exclusion2(){
        const int Width = 2;
        std::unique_ptr<int[]> grid(new int[Size*Size*Size]);
        for(int idxShape = 0 ; idxShape < FP2PExclusion<Width>::SizeShape ; ++idxShape){
            memset(grid.get(), 0, sizeof(int)*Size*Size*Size);

            for(int idxX = 0 ; idxX < Size ; ++idxX){
                for(int idxY = 0 ; idxY < Size ; ++idxY){
                    for(int idxZ = 0 ; idxZ < Size ; ++idxZ){
                        if(FP2PExclusion<Width>::GetShapeIdx(idxX,idxY,idxZ) == idxShape){
                            for(int idxX_neig = FMath::Max(0,idxX-1) ; idxX_neig < FMath::Min(Size,idxX+2) ; ++idxX_neig){
                                for(int idxY_neig = FMath::Max(0,idxY-1) ; idxY_neig < FMath::Min(Size,idxY+2) ; ++idxY_neig){
                                    for(int idxZ_neig = FMath::Max(0,idxZ-1) ; idxZ_neig < FMath::Min(Size,idxZ+2) ; ++idxZ_neig){
                                        uassert(grid[(idxX_neig*Size + idxY_neig)*Size + idxZ_neig] == 0);
                                        grid[(idxX_neig*Size + idxY_neig)*Size + idxZ_neig] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void Exclusion1(){
        const int Width = 1;
        std::unique_ptr<int[]> grid(new int[Size*Size*Size]);
        for(int idxShape = 0 ; idxShape < FP2PExclusion<Width>::SizeShape ; ++idxShape){
            memset(grid.get(), 0, sizeof(int)*Size*Size*Size);

            for(int idxX = 0 ; idxX < Size ; ++idxX){
                for(int idxY = 0 ; idxY < Size ; ++idxY){
                    for(int idxZ = 0 ; idxZ < Size ; ++idxZ){
                        if(FP2PExclusion<Width>::GetShapeIdx(idxX,idxY,idxZ) == idxShape){
                            for(int idxX_neig = FMath::Max(0,idxX-1) ; idxX_neig < idxX ; ++idxX_neig){
                                for(int idxY_neig = FMath::Max(0,idxY-1) ; idxY_neig < idxY ; ++idxY_neig){
                                    for(int idxZ_neig = FMath::Max(0,idxZ-1) ; idxZ_neig < idxZ ; ++idxZ_neig){
                                        uassert(grid[(idxX_neig*Size + idxY_neig)*Size + idxZ_neig] == 0);
                                        grid[(idxX_neig*Size + idxY_neig)*Size + idxZ_neig] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void Middle(){
        std::unique_ptr<int[]> grid(new int[Size*Size*Size]);
        for(int idxShape = 0 ; idxShape < FP2PMiddleExclusion::SizeShape ; ++idxShape){
            memset(grid.get(), 0, sizeof(int)*Size*Size*Size);

            for(int idxX = 0 ; idxX < Size ; ++idxX){
                for(int idxY = 0 ; idxY < Size ; ++idxY){
                    for(int idxZ = 0 ; idxZ < Size ; ++idxZ){
                        if(FP2PMiddleExclusion::GetShapeIdx(idxX,idxY,idxZ) == idxShape){
                            for(int idxX_neig = FMath::Max(0,idxX-1) ; idxX_neig < FMath::Min(Size,idxX+2) ; ++idxX_neig){
                                for(int idxY_neig = FMath::Max(0,idxY-1) ; idxY_neig < FMath::Min(Size,idxY+2) ; ++idxY_neig){
                                    for(int idxZ_neig = FMath::Max(0,idxZ-1) ; idxZ_neig < FMath::Min(Size,idxZ+2) ; ++idxZ_neig){
                                        const int diffx = idxX_neig-idxX;
                                        const int diffy = idxY_neig-idxY;
                                        const int diffz = idxZ_neig-idxZ;
                                        const int idx = (diffx+1)*9 + (diffy+1)*3 + (diffz+1);
                                        if(idx <= 14){
                                            uassert(grid[(idxX_neig*Size + idxY_neig)*Size + idxZ_neig] == 0);
                                            grid[(idxX_neig*Size + idxY_neig)*Size + idxZ_neig] = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    // set test
    void SetTests(){
        AddTest(&TestExclusion::Exclusion2,"Test 2 exclustion");
        AddTest(&TestExclusion::Exclusion1,"Test 1 exclustion");
        AddTest(&TestExclusion::Middle,"Test middle exclustion");
    }
};

// You must do this
TestClass(TestExclusion)



