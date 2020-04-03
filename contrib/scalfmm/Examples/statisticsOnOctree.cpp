// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>

#include "Utils/FParameters.hpp"
#include "Containers/FOctree.hpp"

#include "Components/FBasicCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Components/FBasicParticleContainer.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FMath.hpp"
#include "Files/FFmaGenericLoader.hpp"

#include "Utils/FParameterNames.hpp"

/**
 * @brief Example of ScallFFM's use.
 * @author B. Bramas, O. Coulaud, Q. Khan
 *
 * This example presents the basics of tree traversal with ScallFFM. Reading the
 * source you can find out how to load a tree from an FMA formated file and how
 * to move around the tree. Some statistics are produced to show how to get the
 * information stored.
 * 
 */


int main(int argc, char ** argv) 
{
    typedef double FReal;
    // Typedefs to shorten code
    typedef FBasicCell                                      CellClass;
    typedef FBasicParticleContainer<FReal, 0, FReal >                    ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                   LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass, LeafClass > OctreeClass;

    namespace FPD = FParameterDefinitions;

    // Constant strings for output
    const std::string lhead("[STAT] ");
    const std::string newline = "\n" + lhead;

    //     Parameter handling
    // -------------------------------------------------------------------------
    // Custom parameter
    const FParameterNames LocalOptionHist = {
        {"-histP", "--histogram-stat"},                        // option enabling flags
        "Build the histogram of the particle number per leaf." // option description
    };

    // ScalFFM generic options processing + program's arguments description
    FHelpDescribeAndExit(argc, argv,
                         "Driver to obtain statistics on the octree.",
                         FPD::InputFile,
                         FPD::OctreeHeight,
                         FPD::OctreeSubHeight,
                         FPD::OutputFile, 
                         LocalOptionHist);

    // Octree, input and output parameters
    const int TreeHeight          = FParameters::getValue(argc,argv,FPD::OctreeHeight.options, 4);
    const int SubTreeHeight       = FParameters::getValue(argc,argv,FPD::OctreeSubHeight.options, 1);
    const std::string inFileName  =  FParameters::getStr(argc, argv, FPD::InputFile.options, "../Data/test20k.fma");
    const std::string outFileName =
        FParameters::getStr(argc, argv, FPD::OutputFile.options, "output.txt");

    std::cout << "Parameters  " << std::endl
              << "\tOctree Depth   : " << TreeHeight    << std::endl
              << "\tSubOctree depth: " << SubTreeHeight << std::endl
              << "\tInput file     : " << inFileName    << std::endl
              << "\tOutput file    : " << outFileName   << std::endl
              << std::endl;


    //     Creating and Inserting particles in the tree
    // -------------------------------------------------------------------------
    // The loader reads FMA formated files. The file header is read
    // automatically : we have basic information on the tree
    // structure. Initially, the tree is empty.
    FFmaGenericLoader<FReal> loader(inFileName.c_str());
    OctreeClass tree(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());

    std::cout << newline
              << "Particle file box." << newline
              << "\t\tCentre: " << loader.getCenterOfBox() << newline
              << "\t\tLength: " << loader.getBoxWidth()    
              << std::endl << std::endl;

    FReal  physicalValue;
    FPoint<FReal> particlePosition,
        minPos(loader.getBoxWidth(),loader.getBoxWidth(),loader.getBoxWidth()),
        maxPos(0., 0., 0.);

    // Insertion of particles in the tree.
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        // read next particle in file
        loader.fillParticle(&particlePosition,&physicalValue);
        // insert particle into tree
        tree.insert(particlePosition);
        // Box bounds stats
        minPos.setX( FMath::Min(minPos.getX(), particlePosition.getX()) ) ;
        minPos.setY( FMath::Min(minPos.getY(), particlePosition.getY()) ) ;
        minPos.setZ( FMath::Min(minPos.getZ(), particlePosition.getZ()) ) ;
        maxPos.setX( FMath::Max(maxPos.getX(), particlePosition.getX()) ) ;
        maxPos.setY( FMath::Max(maxPos.getY(), particlePosition.getY()) ) ;
        maxPos.setZ( FMath::Max(maxPos.getZ(), particlePosition.getZ()) ) ;
    }

    std::cout << lhead + "Particle bounding box " << newline
              << "\t\tMin corner: " << minPos << newline
              << "\t\tMax corner: " << maxPos << std::endl << std::endl;


    //     Statistics computation
    // -------------------------------------------------------------------------
    {   // Moving around the leaf level to get information on particles

        long int allLeaves =  (1 << (3 * (TreeHeight-1) )); // maximum number of leaves
        FReal averageParticles = 0.0,  varianceParticles = 0.0;
        FSize nbLeafs = 0, nbPart = 0, nbTPart = 0; // number of particles, total number of particles
        FSize minParticles = std::numeric_limits<FSize>::max(),
            maxParticles = 0;

        // Compute statistics on particles
        // An octree iterator is used to get through the tree.
        OctreeClass::Iterator octreeIterator(&tree);
        // We move to the leftmost leaf.
        octreeIterator.gotoBottomLeft();
        do {
            nbPart            =  octreeIterator.getCurrentListTargets()->getNbParticles() ;
            minParticles      =  FMath::Min(minParticles,nbPart) ;
            maxParticles      =  FMath::Max(maxParticles,nbPart) ;
            nbTPart           += nbPart;
            varianceParticles += FReal(nbPart*nbPart) ;
            ++ nbLeafs;
        } while(octreeIterator.moveRight()); // advancing through the level

        averageParticles  = FReal(nbTPart)/FReal(nbLeafs);
        varianceParticles = varianceParticles/FReal(nbLeafs) - averageParticles*averageParticles;


        std::cout.precision(4);
        std::cout << lhead
                  << "Leaf level      : " << TreeHeight-1 << newline
                  << "Max leaf count  : " << allLeaves << newline
                  << "Non empty leaves: " << nbLeafs 
                  << " (" << 100 * static_cast<FReal>(nbLeafs) / static_cast<FReal>(allLeaves)
                  << "%)" << newline;
        std::cout << "Particles in leaves:" << newline
                  << "\t\tMin     :" << std::setw(8) << minParticles
                  << "\t\tMax     :" << std::setw(8) << maxParticles     << newline
                  << "\t\tAverage :" << std::setw(8) << averageParticles
                  << "\t\tVariance:" << std::setw(8) << varianceParticles 
                  << std::endl << std::endl;


        //  Histogram of particles per leaf
        if( FParameters::existParameter(argc, argv, LocalOptionHist.options) ) {
            std::vector<FSize> hist(maxParticles + 1, 0);

            octreeIterator.gotoBottomLeft();
            do {
                nbPart = octreeIterator.getCurrentListSrc()->getNbParticles();
                ++hist[nbPart] ;
            } while(octreeIterator.moveRight());

            // write data
            std::ofstream outfile(outFileName);
            if( ! outfile.good() ) {
                std::cerr << "Cannot open file " << outFileName << std::endl;
                exit(EXIT_FAILURE);
            }
            outfile << "# Particle per leaf histogram. " << hist.size() << " chunk" << std::endl;
            for(size_t i = 0 ; i < hist.size() ; ++i){
                outfile << i << "  " << hist[i] << std::endl;
            }
        }
        // Statistics on particles done

        FReal averageNeighbors = 0.0, varianceNeighbors = 0.0 ;
        FSize nbBox, minBox = 100000, maxBox=0;
        ContainerClass*  neighborsP2P[27];

        octreeIterator.gotoBottomLeft();
        do {// P2P Neighbors
            nbBox = tree.getLeafsNeighbors(neighborsP2P, 
                                           octreeIterator.getCurrentGlobalCoordinate(),
                                           TreeHeight-1) ;
            // need the current particles and neighbors particles
            minBox             = FMath::Min(minBox,nbBox) ;
            maxBox             = FMath::Max(maxBox,nbBox) ;
            averageNeighbors  += FReal(nbBox);
            varianceNeighbors += FReal(nbBox*nbBox) ;
        } while(octreeIterator.moveRight());

        averageNeighbors  /= FReal(nbLeafs) ;
        varianceNeighbors  = varianceNeighbors/FReal(nbLeafs)-averageNeighbors*averageNeighbors;


        std::cout << lhead
                  << "P2P Neighbors for each leaf "  << newline
                  << "\t\tMin     :" << std::setw(8) << minBox 
                  << "\t\tMax     :" << std::setw(8) << maxBox << newline
                  << "\t\tAverage :" << std::setw(8) << averageNeighbors
                  << "\t\tVariance:" << std::setw(8) << varianceNeighbors
                  << std::endl << std::endl;
    } // End of leaves statistics


    { // Internal cells statistics
        long long int totalCells  = 0;
        long long int totalM2L    = 0;
        long long int totalM2ML2L = 0;
        int nbCellsAtTop    = 0;
        int nbCellsAtBottom = 0;

        // As for the leaf level, we use an iterator got run through the tree.
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        // For each level, see end of loop to climb a level in the tree.
        for(int idxLevel = TreeHeight - 1 ; idxLevel >= 1 ; --idxLevel){

            int nbCellsAtLevel = 0;
            int nbChildrenAtLevel = 0, adaptiveCell = 0, nbChildrenForMyCell;
            int nbNeighborsAtLevel = 0;

            int   nbM2LNeighbors = 0, minM2L = 500, maxM2L = -1;
            FReal averageM2LNeighbors = 0.0, varianceM2LNeighbors = 0.0;
            const CellClass* neighborsM2L[343];

            do {
                ++ nbCellsAtLevel;
                // Count children
                if( idxLevel != TreeHeight - 1 ) {
                    nbChildrenForMyCell = 0;
                    CellClass** children = octreeIterator.getCurrentChild();
                    for( int idxChild = 0 ; idxChild < 8 ; ++ idxChild ) {
                        // Non existing children are NULL
                        if( children[idxChild] )
                            ++ nbChildrenForMyCell;
                    }
                    nbChildrenAtLevel += nbChildrenForMyCell;
                    if( nbChildrenForMyCell > 1 )
                        ++ adaptiveCell;
                }
                //  M2L Neighbors
                nbM2LNeighbors = 
                    tree.getInteractionNeighbors(neighborsM2L,
                                                 octreeIterator.getCurrentGlobalCoordinate(),
                                                 idxLevel);
                nbNeighborsAtLevel += nbM2LNeighbors;

                minM2L               =  FMath::Min(minM2L,nbM2LNeighbors) ;
                maxM2L               =  FMath::Max(maxM2L,nbM2LNeighbors) ;
                averageM2LNeighbors  += FReal(nbM2LNeighbors) ;
                varianceM2LNeighbors += FReal(nbM2LNeighbors*nbM2LNeighbors) ;
            } while(octreeIterator.moveRight());

            averageM2LNeighbors  /= FReal(nbCellsAtLevel);
            varianceM2LNeighbors =  varianceM2LNeighbors / nbCellsAtLevel 
                - averageM2LNeighbors * averageM2LNeighbors;

            totalCells   += (long long int)(nbCellsAtLevel);
            totalM2L     += (long long int)(nbNeighborsAtLevel);
            totalM2ML2L  += (long long int)(nbChildrenAtLevel);
            nbCellsAtTop =  nbCellsAtLevel;
            if( idxLevel == TreeHeight - 1 )
                nbCellsAtBottom = nbCellsAtLevel;


            std::cout << "[STAT] "
                      << "Level " << idxLevel
                      << newline
                      << "\tCell count : " << nbCellsAtLevel
                      << " (adaptive are cells with > 1 children)" << newline
                      << "\t\tAdaptive:" << std::setw(8) << adaptiveCell
                      << "\t\tNon adap:" << std::setw(8) << nbCellsAtLevel-adaptiveCell
                      << newline
                      << "\tM2M/L2L interact:" << std::setw(8) << nbChildrenAtLevel 
                      <<       "\t\tAverage :" << std::setw(8) 
                      << FReal(nbChildrenAtLevel)/FReal(nbCellsAtLevel)
                      << newline
                      << "\tM2L interactions:" << std::setw(8) << nbNeighborsAtLevel
                      << newline;
            std::cout << "\tM2L Neighbors / leaf: "
                      << newline
                      << "\t\tMin     :" << std::setw(8) << minM2L
                      << "\t\tMax     :" << std::setw(8) << maxM2L
                      << newline
                      << "\t\tAverage :" << std::setw(8) << averageM2LNeighbors
                      << "\t\tVariance:" << std::setw(8) << varianceM2LNeighbors
                      << std::endl << std::endl;


            //  Go to next level
            octreeIterator.moveUp();
            // Move to leftmost cell on new level
            octreeIterator.gotoLeft();
        }

        // Global statistics on the octree
        std::cout 
            << lhead + "Global stats" << newline
            << "\tCells = "  << totalCells-nbCellsAtTop << newline
            << "\tM2L interactions  / cell :" << std::setw(8) 
            << totalM2L << newline
            << "\tAverage M2L inter / cell :" << std::setw(8)
            << FReal(totalM2L)/FReal(totalCells) << newline
            << "\tM2M/L2L interactions     :" << std::setw(8)
            << totalM2ML2L << newline
            << "\tAverage M2M/L2L interact :" << std::setw(8)
            << FReal(totalM2ML2L-nbCellsAtTop) / FReal(totalCells-nbCellsAtBottom)
            << std::endl;

    }

    return 0;
}



