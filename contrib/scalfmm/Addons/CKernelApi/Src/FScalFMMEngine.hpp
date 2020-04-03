// ===================================================================================
// Copyright ScalFmm 2014 I
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


/**
 * @file This file contain a class, gathering all the function that
 * can be called by the ScalFMM API. Documentation for each function
 * can be found in the C Header.
 *
 */

#ifndef FSCALFMMENGINE_HPP
#define FSCALFMMENGINE_HPP


#include "Utils/FAssert.hpp"

//For tree
#include "Components/FSimpleLeaf.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

//For interpolation
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

//For lagrange interpolation
// #include "Kernels/Uniform/FUnifCell.hpp"
// #include "Kernels/Uniform/FUnifKernel.hpp"

//For chebyshev Interpolation
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Utils/FAlgorithmTimers.hpp"
#include "Components/FParticleType.hpp"
#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Utils/FTemplate.hpp"
#include "Core/FCoreCommon.hpp"

/**
 * @class FScalFMMEngine
 */
template<class FReal>
class FScalFMMEngine{

protected:
    scalfmm_kernel_type kernelType;

    scalfmm_algorithm Algorithm;
    FVector<bool>* progress;
    int nbPart;
    FAlgorithmTimers * algoTimer;
    FAbstractAlgorithm * abstrct;
public:

    FScalFMMEngine() : Algorithm(multi_thread), progress(nullptr), nbPart(0), algoTimer(nullptr), abstrct(nullptr){
        progress = new FVector<bool>();
    }

    virtual ~FScalFMMEngine() {
        //Do not delete algoTimer because abstract and algoTimer are two pointers on the same thing
        if(abstrct){
            delete abstrct;
            abstrct = nullptr;
            algoTimer = nullptr;
        }
        delete progress;
    }

    //First function displayed there are common function for every
    //kernel
    scalfmm_kernel_type getKernelType(){
        return this->kernelType;
    }

    //To change default algorithm
    void algorithm_config(scalfmm_algorithm config){
        this->Algorithm = config;
    }


    //Functions displayed there are function that are to be redefined
    //by specific Engine

    //Function about the tree
    virtual void build_tree(int TreeHeight,FReal BoxWidth,FReal* BoxCenter,Scalfmm_Cell_Descriptor user_cell_descriptor){
        FAssertLF(0,"Nothing has been done yet, exiting");
    }

    virtual void tree_insert_particles( int NbPositions, FReal * arrayX, FReal * arrayY, FReal * arrayZ, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }

    virtual void tree_insert_particles_xyz( int NbPositions, FReal * XYZ, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }

    virtual void set_physical_values( int nbPhysicalValues, FReal * physicalValues, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }

    virtual void get_physical_values( int nbPhysicalValues, FReal * physicalValues, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }

    virtual void set_physical_values_npart( int nbPhysicalValues,
                                            int* idxOfParticles, FReal * physicalValues, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void get_physical_values_npart( int nbPhysicalValues,
                                            int* idxOfParticles, FReal * physicalValues, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }

    virtual void get_forces_xyz( int nbParts, FReal * forcesToFill, PartType type){
    }

    virtual void get_forces(int nbParts, FReal * fX, FReal* fY, FReal* fZ, PartType type){
    }

    virtual void get_forces_npart(int nbParts, int* idxOfParticles ,FReal * fX, FReal* fY, FReal* fZ, PartType type){
    }

    virtual void get_forces_xyz_npart(int nbParts, int* idxOfParticles, FReal * forcesToFill, PartType type){
    }

    virtual void add_to_positions_xyz(int NbPositions,FReal * updatedXYZ,PartType type){
    }

    virtual void add_to_positions(int NbPositions,FReal * X, FReal * Y , FReal * Z, PartType type){
    }

    virtual void tree_abstract_insert(int NbPartToInsert, int nbAttributeToInsert, int * strideForEachAtt,
                                      FReal* rawDatas){
    }


    //To set initial condition
    virtual void set_forces_xyz( int nbParts, FReal * forcesToFill, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_forces_xyz_npart( int nbParts, int* idxOfParticles, FReal * forcesToFill, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_forces( int nbParts, FReal * fX, FReal* fY, FReal* fZ, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_forces_npart( int nbParts, int* idxOfParticles, FReal * fX, FReal* fY, FReal* fZ, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    //To deal with potential
    virtual void get_potentials( int nbParts, FReal * potentialsToFill, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_potentials( int nbParts, FReal * potentialsToRead, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void get_potentials_npart( int nbParts, int* idxOfParticles, FReal * potentialsToFill, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_potentials_npart( int nbParts, int* idxOfParticles, FReal * potentialsToRead, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void apply_on_each_leaf(Callback_apply_on_leaf function){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }


    /** Test ... */
    struct RunContainer{
        template< int nbAttributeToInsert,class ContainerClass,class LeafClass, class CellClass>
        static void Run(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                        int NbPartToInsert,int * strideForEachAtt,
                        FReal* rawDatas){
            generic_tree_abstract_insert<ContainerClass,LeafClass,CellClass,nbAttributeToInsert>(octree,
                                                                                                 NbPartToInsert,strideForEachAtt,rawDatas);
        }
    };

    template<class ContainerClass,class LeafClass, class CellClass, int nbAttributeToInsert>
    void generic_tree_abstract_insert(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                      int NbPartToInsert,int * strideForEachAtt,
                                      FReal* rawDatas){
        for(FSize idxPart = 0; idxPart<NbPartToInsert ; ++idxPart){
            FPoint<FReal> pos = FPoint<FReal>(rawDatas[0],rawDatas[1],rawDatas[2]);
            MortonIndex index = octree->getMortonFromPosition(pos);
            //Insert with how many attributes ???
            octree->insert(pos,idxPart);
            //Get again the container
            ContainerClass * containerToFill = octree->getLeafSrc(index);//cannot be nullptr
            std::array<FReal,nbAttributeToInsert> arrayOfAttribute;
            for(int idxAtt = 0; idxAtt<nbAttributeToInsert ; ++idxAtt){
                arrayOfAttribute[idxAtt] = rawDatas[3+ strideForEachAtt[idxAtt]];
            }
            int idxToRemove = containerToFill->getNbParticles();
            containerToFill->remove(&idxToRemove,1);
            containerToFill->push(pos,idxPart,arrayOfAttribute);
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_forces_xyz(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                int nbParts, FReal * forcesToFill, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            std::cout << "No meaning to retrieve source forces ... " << std::endl;
        }
        else{ //Targets OR Both
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        forcesToFill[indexes[idxPart]*3+0] = targets->getForcesX()[idxPart];
                        forcesToFill[indexes[idxPart]*3+1] = targets->getForcesY()[idxPart];
                        forcesToFill[indexes[idxPart]*3+2] = targets->getForcesZ()[idxPart];
                        checkCount++;
                    }
                });
        }
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" parts has been read (only "<<checkCount<<")"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" has been read"<< std::endl;}
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_forces_xyz_npart(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                      int nbParts, int* idxOfParticles , FReal * forcesToFill, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            std::cout << "No meaning to retrieve source forces ... " << std::endl;
        }
        else{ //Targets OR Both
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbParts && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                forcesToFill[iterPart*3+0] = targets->getForcesX()[idxPart];
                                forcesToFill[iterPart*3+1] = targets->getForcesY()[idxPart];
                                forcesToFill[iterPart*3+2] = targets->getForcesZ()[idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" parts has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" has been read"<< std::endl;}
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_forces(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                            int nbParts, FReal * fX, FReal* fY, FReal* fZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            std::cout << "No meaning to retrieve source forces ... " << std::endl;
        }
        else{ //Targets OR Both
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        fX[indexes[idxPart]] = targets->getForcesX()[idxPart];
                        fY[indexes[idxPart]] = targets->getForcesY()[idxPart];
                        fZ[indexes[idxPart]] = targets->getForcesZ()[idxPart];
                        checkCount++;
                    }
                });
        }
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" parts has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" has been read"<< std::endl;}
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_forces_nbpart(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                   int nbParts, int* idxOfParticles ,FReal * fX, FReal* fY, FReal* fZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            std::cout << "No meaning to retrieve source forces ... " << std::endl;
        }
        else{ //Targets OR Both
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbParts && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                fX[iterPart] = targets->getForcesX()[idxPart];
                                fY[iterPart] = targets->getForcesY()[idxPart];
                                fZ[iterPart] = targets->getForcesZ()[idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" parts has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" has been read"<< std::endl;}
        }
    }

    //Arranger parts : following function provide a way to move parts
    //inside the tree
    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_add_to_positions_xyz(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                      int NbPositions,FReal * updatedXYZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getWPositions()[0][idxPart] += updatedXYZ[indexes[idxPart]*3+0];
                        sources->getWPositions()[1][idxPart] += updatedXYZ[indexes[idxPart]*3+1];
                        sources->getWPositions()[2][idxPart] += updatedXYZ[indexes[idxPart]*3+2];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getWPositions()[0][idxPart] += updatedXYZ[indexes[idxPart]*3+0];
                        targets->getWPositions()[1][idxPart] += updatedXYZ[indexes[idxPart]*3+1];
                        targets->getWPositions()[2][idxPart] += updatedXYZ[indexes[idxPart]*3+2];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
        update_tree();
    }


    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_add_to_positions(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                  int NbPositions,FReal * X, FReal * Y , FReal * Z, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getWPositions()[0][idxPart] += X[indexes[idxPart]];
                        sources->getWPositions()[1][idxPart] += Y[indexes[idxPart]];
                        sources->getWPositions()[2][idxPart] += Z[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getWPositions()[0][idxPart] += X[indexes[idxPart]];
                        targets->getWPositions()[1][idxPart] += Y[indexes[idxPart]];
                        targets->getWPositions()[2][idxPart] += Z[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
        update_tree();
    }

    //Not yet done
    void add_to_positions_xyz_npart( int NbPositions, int* idxOfParticles, FReal * updatedXYZ, PartType type){
        FAssertLF(0,"Not Yet done ...\n");
    }
    void add_to_positions_npart( int NbPositions, int* idxOfParticles,
                                 FReal * X, FReal * Y , FReal * Z, PartType type){
        FAssertLF(0,"Not Yet done ...\n");
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_set_positions_xyz(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                   int NbPositions, FReal * updatedXYZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getWPositions()[0][idxPart] = updatedXYZ[indexes[idxPart]*3+0];
                        sources->getWPositions()[1][idxPart] = updatedXYZ[indexes[idxPart]*3+1];
                        sources->getWPositions()[2][idxPart] = updatedXYZ[indexes[idxPart]*3+2];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getWPositions()[0][idxPart] = updatedXYZ[indexes[idxPart]*3+0];
                        targets->getWPositions()[1][idxPart] = updatedXYZ[indexes[idxPart]*3+1];
                        targets->getWPositions()[2][idxPart] = updatedXYZ[indexes[idxPart]*3+2];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
        update_tree();
    }


    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_set_positions(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                       int NbPositions, FReal * X, FReal * Y, FReal * Z, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getWPositions()[0][idxPart] = X[indexes[idxPart]];
                        sources->getWPositions()[1][idxPart] = Y[indexes[idxPart]];
                        sources->getWPositions()[2][idxPart] = Z[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getWPositions()[0][idxPart] = X[indexes[idxPart]];
                        targets->getWPositions()[1][idxPart] = Y[indexes[idxPart]];
                        targets->getWPositions()[2][idxPart] = Z[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
        update_tree();
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_set_positions_npart(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                     int NbPositions,int* idxOfParticles,FReal * X, FReal * Y , FReal * Z, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getWPositions()[0][idxPart] = X[indexes[idxPart]];
                        sources->getWPositions()[1][idxPart] = Y[indexes[idxPart]];
                        sources->getWPositions()[2][idxPart] = Z[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getWPositions()[0][idxPart] = X[indexes[idxPart]];
                        targets->getWPositions()[1][idxPart] = Y[indexes[idxPart]];
                        targets->getWPositions()[2][idxPart] = Z[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
        update_tree();
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_set_positions_xyz_npart(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                         int NbPositions,int * idxOfParticles,FReal * updatedXYZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getWPositions()[0][idxPart] = updatedXYZ[indexes[idxPart]*3+0];
                        sources->getWPositions()[1][idxPart] = updatedXYZ[indexes[idxPart]*3+1];
                        sources->getWPositions()[2][idxPart] = updatedXYZ[indexes[idxPart]*3+2];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getWPositions()[0][idxPart] = updatedXYZ[indexes[idxPart]*3+0];
                        targets->getWPositions()[1][idxPart] = updatedXYZ[indexes[idxPart]*3+1];
                        targets->getWPositions()[2][idxPart] = updatedXYZ[indexes[idxPart]*3+2];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
        update_tree();
    }

    virtual void set_positions_xyz_npart( int NbPositions, int* idxOfParticles, FReal * updatedXYZ, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions_npart( int NbPositions, int* idxOfParticles,
                                      FReal * X, FReal * Y , FReal * Z, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions_xyz( int NbPositions, FReal * updatedXYZ, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions( int NbPositions, FReal * X, FReal * Y , FReal * Z, PartType type){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }
    virtual void get_positions_xyz(int NbPositions,FReal * updatedXYZ, PartType type){
    }
    virtual void get_positions(int NbPositions,FReal * X,FReal * Y,FReal * Z, PartType type){
    }
    virtual void get_positions_xyz_npart(int NbPositions,int * idxOfPart, FReal * updatedXYZ, PartType type){
    }
    virtual void get_positions_npart(int NbPositions,int * idxOfPart, FReal * X,FReal * Y,FReal * Z, PartType type){
    }


    //Function to update the tree
    virtual void update_tree(){
        FAssertLF(0,"No tree instancied, exiting ...\n");
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_positions_xyz(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                           int NbPositions, FReal * positionsToFill, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        positionsToFill[indexes[idxPart]*3+0] = sources->getPositions()[0][idxPart];
                        positionsToFill[indexes[idxPart]*3+1] = sources->getPositions()[1][idxPart];
                        positionsToFill[indexes[idxPart]*3+2] = sources->getPositions()[2][idxPart];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        positionsToFill[indexes[idxPart]*3+0] = targets->getPositions()[0][idxPart];
                        positionsToFill[indexes[idxPart]*3+1] = targets->getPositions()[1][idxPart];
                        positionsToFill[indexes[idxPart]*3+2] = targets->getPositions()[2][idxPart];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_positions_xyz_npart(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                                 int NbPositions, int * idxOfParticles, FReal * positionsToFill, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < NbPositions && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                positionsToFill[iterPart] =  sources->getPositions()[0][idxPart];
                                positionsToFill[iterPart] =  sources->getPositions()[1][idxPart];
                                positionsToFill[iterPart] =  sources->getPositions()[2][idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }else {//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < NbPositions && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                positionsToFill[iterPart] =  targets->getPositions()[0][idxPart];
                                positionsToFill[iterPart] =  targets->getPositions()[1][idxPart];
                                positionsToFill[iterPart] =  targets->getPositions()[2][idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_positions(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                       int NbPositions, FReal * X, FReal * Y , FReal * Z, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        X[indexes[idxPart]] = sources->getPositions()[0][idxPart];
                        Y[indexes[idxPart]] = sources->getPositions()[1][idxPart];
                        Z[indexes[idxPart]] = sources->getPositions()[2][idxPart];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        X[indexes[idxPart]] = targets->getPositions()[0][idxPart];
                        Y[indexes[idxPart]] = targets->getPositions()[1][idxPart];
                        Z[indexes[idxPart]] = targets->getPositions()[2][idxPart];
                        checkCount++;
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
    }

    template<class ContainerClass,class LeafClass,class CellClass>
    void generic_get_positions_npart(FOctree<FReal,CellClass,ContainerClass,LeafClass> * octree,
                             int NbPositions, int * idxOfParticles,FReal * X, FReal * Y , FReal * Z,PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < NbPositions && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                X[iterPart] =  sources->getPositions()[0][idxPart];
                                Y[iterPart] =  sources->getPositions()[1][idxPart];
                                Z[iterPart] =  sources->getPositions()[2][idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < NbPositions && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                X[iterPart] =  targets->getPositions()[0][idxPart];
                                Y[iterPart] =  targets->getPositions()[1][idxPart];
                                Z[iterPart] =  targets->getPositions()[2][idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }
        if(checkCount < NbPositions){std::cout << "Not all "<<NbPositions <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > NbPositions){std::cout << "More parts than  "<<NbPositions <<" potentials has been read"<< std::endl;}
        }
    }

    virtual void apply_on_cell(Callback_apply_on_cell function){
    }

    template<class ContainerClass, class CellClass, class LeafClass>
    void generic_apply_on_cell(FOctree<FReal,CellClass,ContainerClass,LeafClass> * tree){
        //Reset forces and potentials
        tree->forEachLeaf([&](LeafClass * leaf){
                ContainerClass * targets = leaf->getTargets();
                FSize nbPartTarget = targets->getNbParticles();
                //Set potential to 0
                FReal * potentialsTarget = targets->getPotentialsArray();
                memset(potentialsTarget,0,sizeof(FReal)*nbPartTarget);
                //Set forces to 0
                FReal * forcesX = targets->getForcesXArray();
                FReal * forcesY = targets->getForcesYArray();
                FReal * forcesZ = targets->getForcesZArray();
                memset(forcesX,0,sizeof(FReal)*nbPartTarget);
                memset(forcesY,0,sizeof(FReal)*nbPartTarget);
                memset(forcesZ,0,sizeof(FReal)*nbPartTarget);
            });

        //Reset multipole and local development
        tree->forEachCell([&](CellClass * cell){
                cell->resetToInitialState();
            });
    }

    //User define Kernel Part
    virtual void user_kernel_config( Scalfmm_Kernel_Descriptor userKernel, void * userDatas){
        FAssertLF(0,"No user kernel defined, exiting ...\n");
    }

    virtual void execute_fmm(){
        FAssertLF(0,"No kernel set, cannot execute anything, exiting ...\n");
    }

    virtual void intern_dealloc_handle(Callback_free_cell userDeallocator){
        FAssertLF(0,"No kernel set, cannot execute anything, exiting ...\n");
    }


    virtual void print_everything(){
    }

    /**
     * Monitoring Function, once the FMM has ended, it's possible to
     * get the time spent in each operator.
     */
    virtual void get_timers(FReal * Timers){
        const FTic * timers = algoTimer->getAllTimers();
        int nbTimers = algoTimer->getNbOfTimerRecorded();
        for(int idTimer = 0; idTimer<nbTimers ; ++idTimer){
            Timers[idTimer] = timers[idTimer].elapsed();
        }
    }

    virtual int get_nb_timers(){
        return algoTimer->getNbOfTimerRecorded();
    }

    virtual void set_upper_limit(int upperLimit){
        FAssertLF(0,"This feature is not available with Chebyshev Kernel, please use your own kernel or do not use it.\n Exiting anyways...\n");
    }

    virtual void create_local_partition(int nbPoints, double * particleXYZ, double ** localArrayFilled, FSize ** indexesFilled, FSize * outputNbPoint){
        FAssertLF(0,"Either no MPI used or wrong initiation function called\n");
    }
    virtual void create_global_partition(int nbPoints, double * particleXYZ, double ** localArrayFilled, FSize ** indexesFilled, FSize * outputNbPoint){
        FAssertLF(0,"Either no MPI used or wrong initiation function called\n");
    }
    virtual void generic_partition(FSize nbThings, size_t sizeofthing, void * arrayOfThing, void ** newArray){
        FAssertLF(0,"Either no MPI used or wrong initiation function called\n");
    }
};

template<class FReal>
struct ScalFmmCoreHandle {

    struct ScalFmmCoreConfig {
        // Read/Write parameter
        int treeHeight;     //  Number of level in the octree
        FReal boxWidth;    // Simulation box size (root level)
        FPoint<FReal> boxCenter; // Center position of the box simulation(FReal[3])
    };

    ScalFmmCoreConfig config;
    FScalFMMEngine<FReal>* engine;
};



extern "C" void scalfmm_build_tree(scalfmm_handle Handle,int TreeHeight,double BoxWidth,double* BoxCenter,Scalfmm_Cell_Descriptor user_cell_descriptor){
    ((ScalFmmCoreHandle<double> *) Handle)->engine->build_tree(TreeHeight,BoxWidth, BoxCenter, user_cell_descriptor);
}

extern "C" void scalfmm_tree_insert_particles(scalfmm_handle Handle, int NbPositions, double * arrayX, double * arrayY, double * arrayZ,
                                              PartType type){
    ((ScalFmmCoreHandle<double> *) Handle)->engine->tree_insert_particles(NbPositions, arrayX, arrayY, arrayZ, type);
}

extern "C" void scalfmm_tree_insert_particles_xyz(scalfmm_handle Handle, int NbPositions, double * XYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->tree_insert_particles_xyz(NbPositions, XYZ,type);
}

extern "C" void scalfmm_set_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_physical_values(nbPhysicalValues, physicalValues, type);
}

extern "C" void scalfmm_get_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_physical_values(nbPhysicalValues, physicalValues, type);
}

extern "C" void scalfmm_set_physical_values_npart(scalfmm_handle Handle, int nbPhysicalValues,
                                                  int* idxOfParticles, double * physicalValues, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_physical_values_npart(nbPhysicalValues,
                                                                       idxOfParticles, physicalValues, type);
}
extern "C" void scalfmm_get_physical_values_npart(scalfmm_handle Handle, int nbPhysicalValues,
                                                  int* idxOfParticles, double * physicalValues, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_physical_values_npart(nbPhysicalValues,
                                                                       idxOfParticles, physicalValues, type);
}

//To get the result
extern "C" void scalfmm_get_forces_xyz(scalfmm_handle Handle, int nbParts, double * forcesToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_forces_xyz(nbParts, forcesToFill, type);
}

extern "C" void scalfmm_get_forces_xyz_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * forcesToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_forces_xyz_npart(nbParts, idxOfParticles, forcesToFill, type);
}
extern "C" void scalfmm_get_forces(scalfmm_handle Handle, int nbParts, double * fX, double* fY, double* fZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_forces(nbParts,fX, fY, fZ, type);
}

extern "C" void scalfmm_get_forces_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_forces_npart(nbParts, idxOfParticles, fX, fY, fZ, type);
}

//To set iniital condition
extern "C" void scalfmm_set_forces_xyz(scalfmm_handle Handle, int nbParts, double * forcesToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_forces_xyz(nbParts, forcesToFill, type);
}

extern "C" void scalfmm_set_forces_xyz_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * forcesToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_forces_xyz_npart(nbParts, idxOfParticles, forcesToFill, type);
}

extern "C" void scalfmm_set_forces(scalfmm_handle Handle, int nbParts, double * fX, double* fY, double* fZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_forces(nbParts, fX, fY, fZ, type);
}

extern "C" void scalfmm_set_forces_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_forces_npart(nbParts, idxOfParticles, fX, fY, fZ, type);
}

//To deal with potential
extern "C" void scalfmm_get_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_potentials(nbParts, potentialsToFill, type);
}

extern "C" void scalfmm_set_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_potentials(nbParts, potentialsToFill, type);
}

extern "C" void scalfmm_get_potentials_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * potentialsToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_potentials_npart(nbParts, idxOfParticles, potentialsToFill, type);
}

extern "C" void scalfmm_set_potentials_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * potentialsToFill, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_potentials_npart(nbParts, idxOfParticles, potentialsToFill, type);
}


// //To deal with positions
// //Out of the box behavior
// extern "C" void scalfmm_out_of_the_box_config(scalfmm_handle Handle,scalfmm_out_of_box_behavior config){
//     ((ScalFmmCoreHandle<double> * ) Handle)->engine->out_of_the_box_config(config);
// }

//Update
extern "C" void scalfmm_add_to_positions_xyz(scalfmm_handle Handle, int NbPositions, double * updatedXYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->add_to_positions_xyz(NbPositions, updatedXYZ, type);
}

extern "C" void scalfmm_add_to_positions(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->add_to_positions(NbPositions, X, Y, Z, type);
}

extern "C" void scalfmm_add_to_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->add_to_positions_xyz_npart(NbPositions, idxOfParticles, updatedXYZ, type);
}

extern "C" void scalfmm_add_to_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles,
                                               double * X, double * Y , double * Z, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->add_to_positions_npart(NbPositions, idxOfParticles, X, Y, Z, type);
}
//Set new positions
extern "C" void scalfmm_set_positions_xyz(scalfmm_handle Handle, int NbPositions, double * updatedXYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_positions_xyz(NbPositions, updatedXYZ, type);
}

extern "C" void scalfmm_set_positions(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_positions(NbPositions, X, Y , Z, type);
}

extern "C" void scalfmm_set_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_positions_xyz_npart(NbPositions, idxOfParticles, updatedXYZ, type);
}

extern "C" void scalfmm_set_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles,
                                            double * X, double * Y , double * Z, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_positions_npart(NbPositions, idxOfParticles, X, Y, Z, type);
}

extern "C" void scalfmm_update_tree(scalfmm_handle Handle){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->update_tree();
}

//Get back positions
extern "C" void scalfmm_get_positions_xyz(scalfmm_handle Handle, int NbPositions, double * updatedXYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_positions_xyz(NbPositions, updatedXYZ, type);
}

extern "C" void scalfmm_get_positions(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_positions(NbPositions, X, Y , Z, type);
}

extern "C" void scalfmm_get_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_positions_xyz_npart(NbPositions, idxOfParticles, updatedXYZ, type);
}

extern "C" void scalfmm_get_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles,
                                            double * X, double * Y , double * Z, PartType type){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_positions_npart(NbPositions, idxOfParticles, X, Y, Z, type);
}

//To choose algorithm
extern "C" void scalfmm_algorithm_config(scalfmm_handle Handle, scalfmm_algorithm config){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->algorithm_config(config);
}

//Executing FMM
extern "C" void scalfmm_execute_fmm(scalfmm_handle Handle){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->execute_fmm();
}

extern "C" void scalfmm_user_kernel_config(scalfmm_handle Handle, Scalfmm_Kernel_Descriptor userKernel, void * userDatas){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->user_kernel_config(userKernel,userDatas);
}

//Monitoring functions
extern "C" void scalfmm_get_timers(scalfmm_handle Handle, double * Timers){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_timers(Timers);
}

extern "C" int scalfmm_get_nb_timers(scalfmm_handle Handle){
    return ((ScalFmmCoreHandle<double> * ) Handle)->engine->get_nb_timers();
}

// extern "C" void scalfmm_tree_abstract_insert(scalfmm_handle Handle, int NbPartToInsert, int nbAttributeToInsert, int * strideForEachAtt,
//                                              double* rawDatas){
//     ((ScalFmmCoreHandle<double> * ) Handle)->engine->tree_abstract_insert(NbPartToInsert,nbAttributeToInsert,strideForEachAtt,rawDatas);
// }

/**
 * These functions are just translating functions.
 */

//< This function fill the childFullPosition[3] with [-1;1] to know the position of a child relatively to
//< its position from its parent
extern "C" void scalfmm_utils_parentChildPosition(int childPosition, int* childFullPosition){
    childFullPosition[2] = (childPosition%2 ? 1 : -1);
    childFullPosition[1] = ((childPosition/2)%2 ? 1 : -1);
    childFullPosition[0] = ((childPosition/4)%2 ? 1 : -1);
}

//< This function fill the childFullPosition[3] with [-3;3] to know the position of a interaction
//< cell relatively to its position from the target
extern "C" void scalfmm_utils_interactionPosition(int interactionPosition, int* srcPosition){
    srcPosition[2] = interactionPosition%7 - 3;
    srcPosition[1] = (interactionPosition/7)%7 - 3;
    srcPosition[0] = (interactionPosition/49)%7 - 3;
}


extern "C" void scalfmm_apply_on_cell(scalfmm_handle Handle, Callback_apply_on_cell function){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->apply_on_cell(function);
}

extern "C" void scalfmm_print_everything(scalfmm_handle Handle){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->print_everything();
}

extern "C" void scalfmm_set_upper_limit(scalfmm_handle Handle, int upperLimit){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->set_upper_limit(upperLimit);
}

extern "C" void scalfmm_apply_on_leaf(scalfmm_handle Handle, Callback_apply_on_leaf function){
    ((ScalFmmCoreHandle<double> * ) Handle)->engine->apply_on_each_leaf(function);
}

#ifdef SCALFMM_USE_MPI
extern "C" void scalfmm_create_local_partition(scalfmm_handle handle, int nbPoints, double * particleXYZ, double ** localArrayFilled,
                                               FSize ** indexesFilled, FSize * outputNbPoint){
    ((ScalFmmCoreHandle<double> * ) handle)->engine->create_local_partition(nbPoints,particleXYZ,localArrayFilled,indexesFilled,outputNbPoint);
}

extern "C" void scalfmm_create_global_partition(scalfmm_handle handle, int nbPoints, double * particleXYZ, double ** localArrayFilled,
                                                FSize ** indexesFilled, FSize * outputNbPoint){
    ((ScalFmmCoreHandle<double> * ) handle)->engine->create_global_partition(nbPoints,particleXYZ,localArrayFilled,indexesFilled,outputNbPoint);
}

extern "C" void scalfmm_generic_partition(scalfmm_handle handle, FSize nbThings, size_t sizeofthing, void * arrayOfThing, void ** newArray){
    ((ScalFmmCoreHandle<double> * ) handle)->engine->generic_partition(nbThings,sizeofthing,arrayOfThing,newArray);
}

#endif

#endif
