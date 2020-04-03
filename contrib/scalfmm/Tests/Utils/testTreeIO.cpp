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
#include <string>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"


#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Files/FTreeIO.hpp"

#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

void usage() {
	std::cout << "Exemple to store and load a tree" << std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "      -help         to see the parameters    " << std::endl
			<<	  "      -depth       the depth of the octree   "<< std::endl
			<<	  "      -subdepth  specifies the size of the sub octree   " << std::endl
			<<     "      -fin    name    name specifies the file of the particle distribution" << std::endl
			<<     "      -fout   name    file name of the octree  " << std::endl;
}
// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Load and store a tree (only the code is interesting).",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::OutputFile);

    typedef double FReal;

	typedef FSphericalCell<FReal>                 CellClass;
	typedef FP2PParticleContainer<FReal>         ContainerClass;

	typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
	typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

	///////////////////////What we do/////////////////////////////
	std::cout << ">> This executable has to be used to load or retrieve an entier tree.\n";
	//////////////////////////////////////////////////////////////
    const unsigned int TreeHeight       = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);

	FTic counter;
    const std::string filenameIN     = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    const std::string filenameOUT = FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "tmp_tree.data");
	std::cout << "Opening : " << filenameIN << "\n";

	FFmaGenericLoader<FReal> loader(filenameIN);
	//
	// -----------------------------------------------------
	CellClass::Init(5);
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

	// -----------------------------------------------------

	std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
	counter.tic();

	for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
		FReal physicalValue = 0.0;
		loader.fillParticle(&particlePosition,&physicalValue);
		tree.insert(particlePosition, physicalValue );
	}

	counter.tac();
	std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

	// -----------------------------------------------------

	std::cout << "Save tree in binary format ..." << std::endl;

    FTreeIO<FReal>::Save<OctreeClass, CellClass, LeafClass, ContainerClass >(filenameOUT.c_str(), tree);

	// -----------------------------------------------------

	std::cout << "Load tree in binary format  ..." << std::endl;

    FTreeIO<FReal>::Load<OctreeClass, CellClass, LeafClass, ContainerClass >(filenameOUT.c_str(),tree);

	return 0;
}



