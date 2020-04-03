// ===================================================================================
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
#include <iomanip>

#include <cstdio>  //printf
#include <cstdlib>
#include <cstring>  //memset

#include <cmath>
#include <algorithm>
#include <string>

#include  "ScalFmmConfig.h"
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FParameterNames.hpp"
#include "../../Src/Files/FIOVtk.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"


#include "../../Src/Files/FDlpolyLoader.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Kernels/P2P/FP2P.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#ifdef SCALFMM_USE_BLAS
// chebyshev kernel
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"
#endif

/** DLpoly particle is used in the gadget program
 * here we try to make the same simulation
 */
template <class FReal>
struct MDParticle {
    FPoint<FReal> position;
	FReal forces[3];
	FReal physicalValue;
	FReal potential;
	int index;
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         ">> This executable has to be used to compute direct interaction either for periodic or non periodic system.\n"
                         ">> options are -depth H -subdepth SH  [-per perdeep,  -noper] -fin filenameIN (-bin)  -fout filenameOUT \n"
                         ">> Recommended files : ../Data/EwalTest_Periodic.run ../Data/EwalTest_NoPeriodic.run.",
                         FParameterDefinitions::OutputFile, FParameterDefinitions::EnabledVerbose,
                         FParameterDefinitions::PeriodicityDisabled, FParameterDefinitions::PeriodicityNbLevels,
                         FParameterDefinitions::OutputFile, FParameterDefinitions::InputFile,
                         FParameterDefinitions::InputBinFormat);
    typedef double FReal;
    typedef FP2PParticleContainerIndexed<FReal>                    ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FInterpMatrixKernelR<FReal>                                              MatrixKernelClass;
  const MatrixKernelClass MatrixKernel;

#ifdef  SCALFMM_USE_BLAS
	// begin Chebyshev kernel
	// accuracy
	const unsigned int ORDER = 12;
	// typedefs
    typedef FChebCell<FReal,ORDER>                                                  CellClass;
    typedef FOctree<FReal,CellClass,ContainerClass,LeafClass>                       OctreeClass;
    typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;

#else
    typedef FSphericalCell<FReal>                                    CellClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<FReal, CellClass, ContainerClass >     KernelClass;
	const int DevP          = FParameters::getValue(argc,argv,"-P", 9);
#endif

    typedef FFmmAlgorithmPeriodic<FReal, OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
	//	typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass >         FmmClassNoPer;
	//////////////////////////////////////////////////////////////

    const int NbLevels         = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options,   4);
    const int SizeSubLevels    = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options,  2);
    const int PeriodicDeep     = FParameters::getValue(argc,argv,FParameterDefinitions::PeriodicityNbLevels.options, 3);
    const std::string filenameIn(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/forceNacl_128_dlpolyPer.bin"));
    const char* const filenameOut = FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "periodicDirect.out");
	//  file for -saveError option
	std::ofstream errorfile("outputEwaldError.txt",std::ios::out);

	FTic counter;

	// -----------------------------------------------------
	//  LOADER
	//  -----------------------------------------------------
	std::cout << "Opening : " << filenameIn << "\n";
    FDlpolyLoader<FReal>  *loader = nullptr ;
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::InputBinFormat.options)){
        loader  = new FDlpolyBinLoader<FReal>(filenameIn.c_str());
	}
	else {
        loader  = new FDlpolyAsciiLoader<FReal>(filenameIn.c_str());
	}

	if(! loader->isOpen()){
		std::cout << "Loader Error, " << filenameIn << " is missing\n";
		return 1;
	}
	OctreeClass tree(NbLevels, SizeSubLevels, loader->getBoxWidth(), loader->getCenterOfBox());
	// ---------------------------------------------------------------------------------
	// Insert particles in the Octree
	// ---------------------------------------------------------------------------------   std::cout << "Creating & Inserting " << loader->getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tWidth : " << loader->getBoxWidth() << " \t center x : " << loader->getCenterOfBox().getX()
	    																														<< " y : " << loader->getCenterOfBox().getY() << " z : " << loader->getCenterOfBox().getZ() << std::endl;
	std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

	counter.tic();
    FPoint<FReal> electricMoment(0.0,0.0,0.0) ;
	// const --> then shared
    MDParticle<FReal> * const particles = new MDParticle<FReal>[loader->getNumberOfParticles()];
    std::memset(particles, 0, sizeof(MDParticle<FReal>) * loader->getNumberOfParticles()) ;
    MDParticle<FReal>* particlesDirect = nullptr;
    particlesDirect = new MDParticle<FReal>[loader->getNumberOfParticles()];
    std::memset(particlesDirect, 0, sizeof(MDParticle<FReal>) * loader->getNumberOfParticles()) ;
	//
	FSize nbParticles = (loader->getNumberOfParticles());
	double totalCharge = 0.0;
	//
	for(FSize idxPart = 0 ; idxPart < loader->getNumberOfParticles() ; ++idxPart){
		//
		loader->fillParticle(&particles[idxPart].position, particles[idxPart].forces,
				&particles[idxPart].physicalValue,&particles[idxPart].index);
		particlesDirect[idxPart].position      = particles[idxPart].position;
		particlesDirect[idxPart].index         = particles[idxPart].index ;
		particlesDirect[idxPart].physicalValue = particles[idxPart].physicalValue;
		//
		totalCharge += particles[idxPart].physicalValue ;
		//
		electricMoment.incX(particles[idxPart].physicalValue*particles[idxPart].position.getX() );
		electricMoment.incY(particles[idxPart].physicalValue*particles[idxPart].position.getY() );
		electricMoment.incZ(particles[idxPart].physicalValue*particles[idxPart].position.getZ() );
		//
		// reset forces and insert in the tree
		//
		tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
	}

	counter.tac();

	std::cout << std::endl;
	std::cout << "Total Charge         = "<< totalCharge <<std::endl;
	std::cout << "Electric Moment      = "<< electricMoment <<std::endl;
	std::cout << "Electric Moment norm = "<< electricMoment.norm2()  <<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << " s)." << std::endl;
	//
	//
	// ---------------------------------------------------------------------------------
	FTreeCoordinate min, max;
    if( FParameters::existParameter(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options) ){
		FmmClass algo(&tree,PeriodicDeep);
		algo.repetitionsIntervals(&min, &max);
		std::cout << "Simulated box: " << algo.extendedBoxWidth()<<std::endl;

	}
	else{
		std::cout << "Simulated box: " << loader->getBoxWidth()<<std::endl;
	}
	// ----------------------------------------------------------------------------------------------------------
	//                                  DIRECT COMPUTATION
	// ----------------------------------------------------------------------------------------------------------
	FReal denergy = 0.0;
	//
	// direct computation
	//
	{
		printf("Compute direct:\n");
		printf("Box [%d;%d][%d;%d][%d;%d]\n", min.getX(), max.getX(), min.getY(),
				max.getY(), min.getZ(), max.getZ());
		//
		//
		// particles       : Input results
		// particlesDirect : Direct results
		//
		counter.tic();
		denergy = 0.0 ;
//#pragma omp parallel for shared(loader,nbParticles,min, max,particlesDirect) schedule(static,1)
#pragma omp parallel shared(loader,nbParticles,min, max,particlesDirect)
		{
#pragma omp for schedule(static,1)
			for(int idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
				//
				// compute with all other except itself
				//
				// Compute force and potential between  particles[idxTarget] and particles inside the box
				//
				for(int idxOther = 0; idxOther < nbParticles ; ++idxOther){
					if( idxOther != idxTarget ){
						FP2P::NonMutualParticles(
								particles[idxOther].position.getX(), particles[idxOther].position.getY(),
								particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
								particlesDirect[idxTarget].position.getX(), particlesDirect[idxTarget].position.getY(),
								particlesDirect[idxTarget].position.getZ(),particlesDirect[idxTarget].physicalValue,
								&particlesDirect[idxTarget].forces[0],&particlesDirect[idxTarget].forces[1],
								&particlesDirect[idxTarget].forces[2],&particlesDirect[idxTarget].potential, &MatrixKernel);
					}
				}
				//
				// compute with image boxes
				//
				for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
					for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
						for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
							{
								if(idxX == 0 && idxY == 0 && idxZ == 0) continue;

                                const FPoint<FReal> offset(loader->getBoxWidth() * FReal(idxX),
										loader->getBoxWidth() * FReal(idxY),
										loader->getBoxWidth() * FReal(idxZ));
								//
								for(int idxSource = 0 ; idxSource < nbParticles ; ++idxSource){
                                    MDParticle<FReal> source = particles[idxSource];
									source.position += offset;
									FP2P::NonMutualParticles(
											source.position.getX(), source.position.getY(),source.position.getZ(),source.physicalValue,
											particlesDirect[idxTarget].position.getX(), particlesDirect[idxTarget].position.getY(),particlesDirect[idxTarget].position.getZ(),particlesDirect[idxTarget].physicalValue,
											&particlesDirect[idxTarget].forces[0],&particlesDirect[idxTarget].forces[1],&particlesDirect[idxTarget].forces[2],&particlesDirect[idxTarget].potential, &MatrixKernel
									);
								}
							}
						}
					}
				} // End nested loop - compute with image boxes
			} // end for
		// Compute the energy
//#pragma omp  for shared(nbParticles,particlesDirect) reduction(+:denergy)
#pragma omp  for reduction(+:denergy)
		for(int idx = 0 ; idx < nbParticles ; ++idx){
			denergy +=  particlesDirect[idx].potential*particlesDirect[idx].physicalValue ;
		}
} // end pragma parallel

		denergy *= 0.5 ;
		counter.tac();
		//
		printf("Energy DIRECT=   %.14e\n",denergy);
		//
		std::cout << "Done  " << "(@DIRECT Algorithm = " << counter.elapsed() << " s)." << std::endl;
		std::cout << "\n"<< "END DIRECT "
				<< "-------------------------------------------------------------------------"
				<< std::endl << std::endl ;
	} // END DIRECT

	//
	// ----------------------------------------------------------------
	//  Save direct computation in binary format
	// write binary output file
	//
	std::cout << "Generate " << filenameOut <<"  from input file" << std::endl;
	std::fstream fileout(filenameOut,std::ifstream::in|std::ifstream::out| std::ios::binary| std::ios::trunc);
	if(!fileout) {
		std::cout << "Cannot open file."<< std::endl;
		return 1;
	}

	//
	std::cout << " nbParticles: " << nbParticles <<"  " << sizeof(nbParticles) <<std::endl;
	std::cout << " denergy: " << denergy <<"  " << sizeof(denergy) <<std::endl;
	std::cout << " Box size: " << loader->getBoxWidth() << "  " << sizeof(loader->getBoxWidth())<<std::endl;
	double boxsize[3] ; boxsize[0] = boxsize[1]= boxsize[2]=loader->getBoxWidth();
	int PER[4] ;
    if( FParameters::existParameter(argc, argv, FParameterDefinitions::PeriodicityDisabled.options) ){
		PER[0] = 0 ;
		PER[1] = PER[2] = PER[3] = 0 ;
	}
	else {
		PER[0] = 1 ;
		PER[1] = PER[2] = PER[3] = PeriodicDeep ;
	}
	fileout.write((char*)&nbParticles,sizeof(nbParticles));
	fileout.write((char*)&boxsize,sizeof(double)*3);
	fileout.write((char*)&PER,sizeof(int)*4);
	fileout.write((char*)&denergy,sizeof(denergy));
	//
    fileout.write ((char*)&particlesDirect[0], sizeof(MDParticle<FReal>)*nbParticles);
	fileout.flush();
	//
	// end generate
	// -----------------------------------------------------
	//
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::EnabledVerbose.options)){
		denergy = 0 ;
		for(int idx = 0 ; idx < nbParticles ; ++idx){
			std::cout << ">> index " << particlesDirect[idx].index << std::endl;
			std::cout << " x   " << particlesDirect[idx].position.getX() << " y  " << particlesDirect[idx].position.getY() << " z  " << particlesDirect[idx].position.getZ() << std::endl;
			std::cout << " fx  " << particlesDirect[idx].forces[0]       << " fy " << particlesDirect[idx].forces[1]       << " fz " << particlesDirect[idx].forces[2] << std::endl;
			std::cout << " Q   " << particlesDirect[idx].physicalValue   << " V  " << particlesDirect[idx].potential << std::endl;
			std::cout << "\n";
			denergy +=  particlesDirect[idx].potential*particlesDirect[idx].physicalValue ;
		}
	}
    std::cout << " ENERGY " << denergy*0.5 << std::endl;
	delete[] particles;
	delete[] particlesDirect;

	return 0;
}



