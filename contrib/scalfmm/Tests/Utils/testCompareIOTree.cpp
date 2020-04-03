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

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FTreeIO.hpp"

#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

// Simply create particles and try the kernels
int main(int argc, char ** argv){

    FHelpDescribeAndExit(argc, argv,
                         "Load octrees that have been saved and compare everything from leaves to cells.\n"
                         "Using the Spherical Harmonics old kernel.",
                         FParameterDefinitions::SHDevelopment, FParameterDefinitions::InputFileOne,
                         FParameterDefinitions::InputFileTwow
                         );

    typedef double FReal;
    typedef FSphericalCell<FReal>                 CellClass;
    typedef FP2PParticleContainer<FReal>         ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to compare two trees.\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 8);

    // -----------------------------------------------------
    CellClass::Init(DevP, true);
    OctreeClass tree1(5, 3, 0, FPoint<FReal>());
    OctreeClass tree2(5, 3, 0, FPoint<FReal>());

    // -----------------------------------------------------
    const char* const filename1 = FParameters::getStr(argc,argv,FParameterDefinitions::InputFileOne.options, "tree.data");
    const char* const filename2 = FParameters::getStr(argc,argv,FParameterDefinitions::InputFileOne.options, "dtree.data");
    std::cout << "Compare tree " << filename1 << " and " << filename2 << std::endl;

    FTreeIO<FReal>::Load<OctreeClass, CellClass, LeafClass, ContainerClass >(filename1, tree1);
    FTreeIO<FReal>::Load<OctreeClass, CellClass, LeafClass, ContainerClass >(filename2, tree2);

    // -----------------------------------------------------
    std::cout << "Check Result\n";
    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator1(&tree1);
        octreeIterator1.gotoBottomLeft();

        OctreeClass::Iterator octreeIterator2(&tree2);
        octreeIterator2.gotoBottomLeft();

        int nbLeaves = 0;

        do{
            if( octreeIterator1.getCurrentGlobalIndex() != octreeIterator2.getCurrentGlobalIndex()){
                std::cout << "Index is different\n";
                break;
            }

            if( octreeIterator1.getCurrentListSrc()->getNbParticles() != octreeIterator2.getCurrentListSrc()->getNbParticles()){
                std::cout << "Number of particles different on leaf " << octreeIterator1.getCurrentGlobalIndex() <<
                             " tree1 " << octreeIterator1.getCurrentListSrc()->getNbParticles() <<
                             " tree2 " << octreeIterator2.getCurrentListSrc()->getNbParticles() << std::endl;
            }

            nbLeaves += 1;

            if( octreeIterator1.moveRight() ){
                if( !octreeIterator2.moveRight() ){
                    std::cout << "Not the same number of leaf, tree2 end before tree1\n";
                    break;
                }
            }
            else {
                if( octreeIterator2.moveRight() ){
                    std::cout << "Not the same number of leaf, tree1 end before tree2\n";
                }
                break;
            }

        } while(true);

        std::cout << "There are " << nbLeaves << " leaves ...\n";
    }
    { // Ceck if there is number of NbPart summed at level 1
        OctreeClass::Iterator octreeIterator1(&tree1);
        octreeIterator1.gotoBottomLeft();

        OctreeClass::Iterator octreeIterator2(&tree2);
        octreeIterator2.gotoBottomLeft();

        for(int idxLevel = tree1.getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
            int nbCells = 0;
            do{
                if( octreeIterator1.getCurrentGlobalIndex() != octreeIterator2.getCurrentGlobalIndex()){
                    std::cout << "Index is different\n";
                    break;
                }

                const CellClass*const cell1 = octreeIterator1.getCurrentCell();
                const CellClass*const cell2 = octreeIterator2.getCurrentCell();

                FReal cumul = 0;
                for(int idx = 0; idx < FSphericalCell<FReal>::GetPoleSize(); ++idx){
                    cumul += FMath::Abs( cell1->getMultipole()[idx].getImag() - cell2->getMultipole()[idx].getImag() );
                    cumul += FMath::Abs( cell1->getMultipole()[idx].getReal() - cell2->getMultipole()[idx].getReal() );
                }
                if( cumul > 0.00001 || FMath::IsNan(cumul)){
                    std::cout << "Pole Data are different. Cumul " << cumul << " at level " << idxLevel
                              << " index is " << octreeIterator1.getCurrentGlobalIndex() << std::endl;
                }
                cumul = 0;
                for(int idx = 0; idx < FSphericalCell<FReal>::GetLocalSize(); ++idx){
                    cumul += FMath::Abs( cell1->getLocal()[idx].getImag() - cell2->getLocal()[idx].getImag() );
                    cumul += FMath::Abs( cell1->getLocal()[idx].getReal() - cell2->getLocal()[idx].getReal() );
                }
                if( cumul > 0.00001 || FMath::IsNan(cumul)){
                    std::cout << "Local Data are different. Cumul " << cumul << " at level " << idxLevel
                              << " index is " << octreeIterator1.getCurrentGlobalIndex() << std::endl;
                }

                nbCells += 1;
                if( octreeIterator1.moveRight() ){
                    if( !octreeIterator2.moveRight() ){
                        std::cout << "Not the same number of leaf, tree2 end before tree1\n";
                        break;
                    }
                }
                else {
                    if( octreeIterator2.moveRight() ){
                        std::cout << "Not the same number of leaf, tree1 end before tree2\n";
                    }
                    break;
                }
            } while(true);

            octreeIterator1.moveUp();
            octreeIterator1.gotoLeft();

            octreeIterator2.moveUp();
            octreeIterator2.gotoLeft();

            std::cout << "There are " << nbCells << " cells at level " << idxLevel << " ...\n";
        }
    }

    std::cout << "Done\n";

    return 0;
}




