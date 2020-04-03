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

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMpi.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FTreeMpiCsvSaver.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


template <class FReal>
class VelocityContainer : public FP2PParticleContainer<FReal> {
    typedef FP2PParticleContainer<FReal> Parent;

    FVector<FPoint<FReal>> velocities;

public:
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FPoint<FReal>& velocity, Args... args){
        Parent::push(inParticlePosition, args... );
        velocities.push(velocity);
    }

    const FVector<FPoint<FReal>>& getVelocities() const{
        return velocities;
    }

    FVector<FPoint<FReal>>& getVelocities() {
        return velocities;
    }

    void fillToCsv(const FSize partIdx, FReal values[4]) const {
        values[0] = Parent::getPositions()[0][partIdx];
        values[1] = Parent::getPositions()[1][partIdx];
        values[2] = Parent::getPositions()[2][partIdx];
        values[3] = Parent::getPotentials()[partIdx];
    }
};


template <class FReal>
class GalaxyLoader : public FFmaGenericLoader<FReal> {
public:
    GalaxyLoader(const std::string & filename) : FFmaGenericLoader<FReal>(filename) {
    }

    void fillParticle(FPoint<FReal>* position, FReal* physivalValue, FPoint<FReal>* velocity){
        FReal x,y,z,data, vx, vy, vz;
        (*FFmaGenericLoader<FReal>::file)  >> x >> y >> z >> data >> vx >> vy >> vz;
        position->setPosition(x,y,z);
        *physivalValue = (data);
        velocity->setPosition(vx,vy,vz);
    }
};

template <class FReal>
struct TestParticle{
    FPoint<FReal> position;
    FReal physicalValue;
    FReal potential;
    FReal forces[3];
    FPoint<FReal> velocity;
    const FPoint<FReal>& getPosition(){
        return position;
    }
};

template < class FReal,class ParticleClass>
class Converter {
public:
    template <class ContainerClass>
    static ParticleClass GetParticle(ContainerClass* containers, const int idxExtract){
        const FReal*const positionsX = containers->getPositions()[0];
        const FReal*const positionsY = containers->getPositions()[1];
        const FReal*const positionsZ = containers->getPositions()[2];
        const FReal*const forcesX = containers->getForcesX();
        const FReal*const forcesY = containers->getForcesY();
        const FReal*const forcesZ = containers->getForcesZ();
        const FReal*const physicalValues = containers->getPhysicalValues();
        const FReal*const potentials = containers->getPotentials();
        FVector<FPoint<FReal>> velocites = containers->getVelocities();

        TestParticle<FReal> part;
        part.position.setPosition( positionsX[idxExtract],positionsY[idxExtract],positionsZ[idxExtract]);
        part.physicalValue = physicalValues[idxExtract];
        part.forces[0] = forcesX[idxExtract];
        part.forces[1] = forcesY[idxExtract];
        part.forces[2] = forcesZ[idxExtract];
        part.potential = potentials[idxExtract];
        part.velocity  = velocites[idxExtract];

        return part;
    }

    template <class OctreeClass>
    static void Insert(OctreeClass* tree, const ParticleClass& part){
        tree->insert(part.position , part.velocity, part.physicalValue, part.forces[0],
                part.forces[1],part.forces[2],part.potential);
    }
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Convert the data from a file into a csv file to load into Paraview for example.\n"
                         "It puts the file into the /tmp dir and the code is an example of using FTreeMpiCsvSaver.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::InputFile);

    typedef double FReal;
    typedef FBasicCell              CellClass;
    typedef VelocityContainer<FReal>  ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////
    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 6);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);

    GalaxyLoader<FReal> loader(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/galaxy.fma"));

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    {
        FPoint<FReal> position, velocity;
        FReal physicalValue;

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&position, &physicalValue, &velocity);
            if( (idxPart+1) % (app.global().processId()+1) == 0) tree.insert(position,velocity,physicalValue);
        }
    }

    // -----------------------------------------------------

    {
        FTreeMpiCsvSaver<FReal, OctreeClass, ContainerClass> saver("/tmp/test%d.csv", app.global() , false);
        saver.exportTree(&tree);
    }

    // -----------------------------------------------------

    {
        FTreeMpiCsvSaver<FReal, OctreeClass, ContainerClass> saver("/tmp/htest%d.csv", app.global() , true);
        saver.exportTree(&tree);
    }

    return 0;
}
