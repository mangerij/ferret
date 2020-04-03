// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================

#include <iostream>
#include <cstdio>

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../Src/FmmApi.h"

/** This program show an example of use of the fmm api
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM API\n";
    //////////////////////////////////////////////////////////////

    int NbLevels      = FParameters::getValue(argc,argv,"-h", 7);
    int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const int NbPart  = FParameters::getValue(argc,argv,"-nb", 2000000);
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoader<FReal> loader(NbPart, 1, FPoint<FReal>(0.5,0.5,0.5), 1);

    void* FmmCoreHandle;
    FmmCore_init(&FmmCoreHandle);

    FmmCore_setParameter(FmmCoreHandle, FMMCORE_TREE_HEIGHT, &NbLevels);
    FReal boxWidth = loader.getBoxWidth();
    FmmCore_setParameter(FmmCoreHandle, FMMCORE_ROOT_BOX_WIDTH, &boxWidth);
    FmmCore_setParameter(FmmCoreHandle, FMMCORE_ROOT_BOX_CENTER, loader.getCenterOfBox().getDataValue());

    void* FmmKernelHandle;
    FmmKernel_init(FmmCoreHandle, &FmmKernelHandle);

    FmmCore_setKernelData(FmmCoreHandle, FmmKernelHandle);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    FSize r.getNumberOfParticles();
    FReal* potentials = new FReal[nbPart];
    FReal* positions = new FReal[nbPart*3];

    {
        FPoint<FReal> part;
        const FReal physicalValue = 0.1;

        for(int idx = 0 ; idx < nbPart ; ++idx){
            loader.fillParticle(&part);

            potentials[idx] = physicalValue;
            positions[3*idx] = part.getX();
            positions[3*idx+1] = part.getY();
            positions[3*idx+2] = part.getZ();
        }

        FmmCore_setPositions(FmmCoreHandle, &nbPart, positions);
        FmmCore_setPotentials(FmmCoreHandle, potentials);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    FmmCore_doComputation(FmmCoreHandle);

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Getting results ..." << std::endl;
    counter.tic();

    FReal* fields = new FReal[4*nbPart];

    FmmCore_getField(FmmCoreHandle, fields);

    counter.tac();
    std::cout << "Done  " << "(@retrieve = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Computing direct ..." << std::endl;
    counter.tic();

    FReal* fieldsdirect = new FReal[4*nbPart];
    FMemUtils::setall(fieldsdirect, 0.0, nbPart * 4);

    for(int idxTarget = 0 ; idxTarget < nbPart ; ++idxTarget){

        for(int idxSource = idxTarget + 1 ; idxSource < nbPart ; ++idxSource){
            FReal dx = positions[idxSource*3] - positions[idxTarget*3];
            FReal dy = positions[idxSource*3+1] - positions[idxTarget*3+1];
            FReal dz = positions[idxSource*3+2] - positions[idxTarget*3+2];

            FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
            FReal inv_distance = FMath::Sqrt(inv_square_distance);

            inv_square_distance *= inv_distance;
            inv_square_distance *= potentials[idxTarget] * potentials[idxSource];

            dx *= inv_square_distance;
            dy *= inv_square_distance;
            dz *= inv_square_distance;

            fieldsdirect[4*idxTarget+1] += dx;
            fieldsdirect[4*idxTarget+2] += dy;
            fieldsdirect[4*idxTarget+3] += dz;
            fieldsdirect[4*idxTarget+0] += inv_distance * potentials[idxSource];

            fieldsdirect[4*idxSource+1] += -dx;
            fieldsdirect[4*idxSource+2] += -dy;
            fieldsdirect[4*idxSource+3] += -dz;
            fieldsdirect[4*idxSource+0] += inv_distance * potentials[idxTarget];

        }
    }


    counter.tac();
    std::cout << "Done  " << "(@direct = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////


    std::cout << "Comparing results ..." << std::endl;
    counter.tic();

    FMath::FAccurater<FReal> forces;
    FMath::FAccurater<FReal> potential;

    for(int idx = 0 ; idx < nbPart ; ++idx){
        if(idx < 10){
            std::cout << fieldsdirect[4*idx]  << " fmm " << fields[4*idx] << std::endl;
        }
        forces.add(     fieldsdirect[4*idx]      ,fields[4*idx]);
        forces.add(     fieldsdirect[4*idx+1]    ,fields[4*idx+1]);
        forces.add(     fieldsdirect[4*idx+2]    ,fields[4*idx+2]);
        potential.add(  fieldsdirect[4*idx+3]    ,fields[4*idx+3]);
    }

    counter.tac();
    std::cout << "Done  " << "(@comparing = " << counter.elapsed() << "s)." << std::endl;

    std::cout << "\tForces inf " << forces.getInfNorm() << " normL2 " << forces.getL2Norm() << std::endl;
    std::cout << "\tPotential inf " << potential.getInfNorm() << " normL2 " << potential.getL2Norm() << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    delete[] fieldsdirect;
    delete[] fields;
    delete[] potentials;
    delete[] positions;

    FmmCore_free(FmmCoreHandle);
    FmmKernel_free(FmmKernelHandle);

    return 0;
}



