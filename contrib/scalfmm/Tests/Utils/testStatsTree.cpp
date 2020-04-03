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

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"


#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Print a lot of information about the tree build from a given file.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;

    typedef FBasicParticleContainer<FReal,0,FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to show some stat about the tree.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    std::cout << "Opening : " << filename << "\n";

    FFmaGenericLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating and Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition );
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;


    // -----------------------------------------------------

    { // get stats
        {
            std::cout << "[STAT] level is " << NbLevels << std::endl;
            std::cout << "[STAT] potentials leafs number is " << (1 << (3* (NbLevels-1) )) << std::endl;

            FReal averageParticles = 0;
            {
                int nbLeafs = 0;
                OctreeClass::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    averageParticles += FReal(octreeIterator.getCurrentListTargets()->getNbParticles());
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                averageParticles /= FReal(nbLeafs);

                std::cout << "[STAT] Nb leafs : " << nbLeafs << std::endl;
            }
            std::cout << "[STAT] Average particles on leafs = " << averageParticles << std::endl;

            FReal varianceParticles = 0;
            {
                int nbLeafs = 0;
                OctreeClass::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    varianceParticles += FReal(octreeIterator.getCurrentListTargets()->getNbParticles() * octreeIterator.getCurrentListTargets()->getNbParticles());
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                varianceParticles /= FReal(nbLeafs);
            }
            std::cout << "[STAT] Variances of particles on leafs is = " << (varianceParticles - (averageParticles*averageParticles)) << std::endl;
        }

        {
            FReal averageNeighbors = 0;
            int nbLeafs = 0;
            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                ContainerClass* neighbors[27];
                // need the current particles and neighbors particles
                averageNeighbors += FReal(tree.getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),NbLevels-1));
                ++nbLeafs;
            } while(octreeIterator.moveRight());
            std::cout << "[STAT] Average neighbors for each leafs = " << (averageNeighbors/FReal(nbLeafs)) << std::endl;
        }

        {
            long long int totalCells = 0;
            long long int totalM2L = 0;
            long long int totalM2ML2L = 0;

            int nbCellsAtTop = 0;
            int nbCellsAtBottom = 0;

            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel){

                int nbCellsAtLevel = 0;
                int nbChildAtLevel = 0;
                int nbNeighborsAtLevel = 0;

                do{
                    ++nbCellsAtLevel;

                    if( idxLevel != NbLevels - 1 ){
                        FBasicCell** child = octreeIterator.getCurrentChild();
                        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                            if(child[idxChild]) ++nbChildAtLevel;
                        }
                    }

                    const FBasicCell* neighbors[343];
                    nbNeighborsAtLevel += tree.getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);

                } while(octreeIterator.moveRight());

                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                std::cout << "[STAT] Level = " << idxLevel << "\n";
                std::cout << "[STAT] >> Nb Cells = " << nbCellsAtLevel << "\n";
                std::cout << "[STAT] >> Nb M2M/L2L interactions = " << nbChildAtLevel << "\n";
                std::cout << "[STAT] >> Average M2M/L2L interactions = " << FReal(nbChildAtLevel)/FReal(nbCellsAtLevel) << "\n";
                std::cout << "[STAT] >> Nb M2L interactions = " << nbNeighborsAtLevel << "\n";
                std::cout << "[STAT] >> Average M2L interactions = " << FReal(nbNeighborsAtLevel)/FReal(nbCellsAtLevel) << "\n";

                totalCells += (long long int)(nbCellsAtLevel);
                totalM2L += (long long int)(nbNeighborsAtLevel);
                totalM2ML2L += (long long int)(nbChildAtLevel);
                nbCellsAtTop = nbCellsAtLevel;
                if( idxLevel == NbLevels - 1 ) nbCellsAtBottom = nbCellsAtLevel;
            }

            std::cout << "[STAT] For all the tree\n";
            std::cout << "[STAT] >> Total Nb Cells = " << totalCells-nbCellsAtTop << "\n";
            std::cout << "[STAT] >> Total Nb M2M/L2L interactions = " << totalM2ML2L << "\n";
            std::cout << "[STAT] >> Total Average M2M/L2L interactions = " << FReal(totalM2ML2L-nbCellsAtTop)/FReal(totalCells-nbCellsAtBottom) << "\n";
            std::cout << "[STAT] >> Total Nb M2L interactions per cell = " << totalM2L << "\n";
            std::cout << "[STAT] >> Total Average M2L interactions per cell = " << FReal(totalM2L)/FReal(totalCells) << "\n";

       }
    }

    // -----------------------------------------------------


    return 0;
}



