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
 * @file This file contains a class that inherits from FScalFMMEngine,
 * and will implement the API functions for Interpolations kernels.
 */

#ifndef FINTERENGINE_HPP
#define FINTERENGINE_HPP

#include "FScalFMMEngine.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Components/FTypedLeaf.hpp"

#include "Arranger/FOctreeArranger.hpp"
#include "Arranger/FArrangerPeriodic.hpp"
#include "Arranger/FBasicParticleContainerIndexedMover.hpp"
#include "Arranger/FParticleTypedIndexedMover.hpp"
#include "Extensions/FExtendCellType.hpp"

#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmSectionTask.hpp"
#include "Core/FFmmAlgorithmTask.hpp"
#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmPeriodic.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"



/**
 * @class FInterEngine implements API for Interpolations kernels, its
 * templates can be ChebCell/ChebKernel or UnifCell/UnifKernel
 */
template<class FReal, class InterCell,class InterKernel,class LeafClass,
         class MatrixKernelClass = FInterpMatrixKernelR<FReal> >
class FInterEngine : public FScalFMMEngine<FReal>{
private:
    //Typedefs
    typedef FP2PParticleContainerIndexed<FReal>           ContainerClass;
    typedef FTypedLeaf<FReal,ContainerClass>                    LeafClassTyped;

    //Typedef on the Typed OctreeClass
    typedef FOctree<FReal,InterCell,ContainerClass,LeafClass>           OctreeClass;


    //Pointer to the kernel to be executed
    InterKernel * kernel;
    MatrixKernelClass * matrix;
    //Link to the tree
    OctreeClass * octree;

    // ArrangerClass * arranger;



public:
    /**
     * @brief Constructor : build the tree and the interpolation
     * kernel
     * @param TreeHeight Height of the tree
     * @param BoxWidth box Width
     * @param BoxCenter double[3] coordinate of the center of the
     * simulation box
     */
    FInterEngine(scalfmm_kernel_type KernelType, scalfmm_algorithm algo) :
        kernel(nullptr), matrix(nullptr), octree(nullptr)/*,arranger(nullptr)*/{
        FScalFMMEngine<FReal>::kernelType = KernelType;
        FScalFMMEngine<FReal>::Algorithm = algo;
    }

    void build_tree(int TreeHeight, FReal BoxWidth , FReal * BoxCenter,User_Scalfmm_Cell_Descriptor notUsedHere){
        octree = new OctreeClass(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint<FReal>(BoxCenter));
        this->matrix = new MatrixKernelClass();
        this->kernel = new InterKernel(TreeHeight,BoxWidth,FPoint<FReal>(BoxCenter),matrix);
    }


    //TODO free kernel too
    ~FInterEngine(){
        delete matrix;
        if(octree){
            delete octree;
        }
        if(kernel){
            delete kernel;
        }
        // if(arranger){
        //     delete arranger;
        // }
    }

    //Inserting array of position
    //Need to be disabled if Source/Target is used
    void tree_insert_particles_xyz(int NbPositions, FReal * XYZ, PartType type){
        if(type == BOTH){
            for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                octree->insert(FPoint<FReal>(&XYZ[3*idPart]),idPart);
            }
            FScalFMMEngine<FReal>::nbPart += NbPositions;
        }else{
            if(type==SOURCE){
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(&XYZ[3*idPart]),FParticleTypeSource,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }else{
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(&XYZ[3*idPart]),FParticleTypeTarget,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }
        }
    }

    //Inserting arrayS of position
    //Need to be disabled if Source/Target is used
    void tree_insert_particles(int NbPositions, FReal * X, FReal * Y, FReal * Z, PartType type){
        if(type == BOTH){
            for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),idPart);
            }
            FScalFMMEngine<FReal>::nbPart += NbPositions;
        }else{
            if(type==SOURCE){
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),FParticleTypeSource,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }else{
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),FParticleTypeTarget,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }
        }
    }

    // void tree_abstract_insert(int NbPartToInsert, int nbAttributeToInsert, int * strideForEachAtt,
    //                           FReal* rawDatas){
    //     FAssertLF(nbAttributeToInsert > 2,"Need space to store positions, thus nbAttributeToInsert must be >= 3\nExiting ... \n");
    //     FAssertLF(nbAttributeToInsert < 15,"Cannot instanciate more than 15 Attribute per Particules\n");
    //     FRunIf::Run<int,3,15,1,RunContainer>(nbAttributeToInsert,);
    //     generic_tree_abstract_insert<ContainerClass,LeafClass,InterCell,nbAttributeToInsert>(octree,
    //                                                                                          NbPartToInsert,strideForEachAtt,rawDatas);
    // }

    //Set the physical values
    void set_physical_values(int nbPhysicalValues,FReal * physicalValues, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            octree->forEachLeaf([&] (LeafClass * leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(int idx=0 ; idx<nbPartThere ; ++idx){
                        sources->getPhysicalValues()[idx] = physicalValues[indexes[idx]];
                        ++checkCount;
                    }
                });
        }
        else{ // type must be equal to TARGETS or BOTH
            octree->forEachLeaf([&] (LeafClass * leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(int idx=0 ; idx<nbPartThere ; ++idx){
                        targets->getPhysicalValues()[idx] = physicalValues[indexes[idx]];
                        ++checkCount;
                    }
                });
        }
        if(checkCount < nbPhysicalValues){std::cout << "Not all "<<nbPhysicalValues <<" parts has been set (only "<<checkCount<<")" << std::endl;}
        else{
            if(checkCount > nbPhysicalValues){std::cout << "More parts than  "<<nbPhysicalValues <<" has been set"<< std::endl;}
        }
    }

    //Set only a subpart of physical values
    //Algorithm : loop over each leaf, and then search in user array
    //if any index matches
    void set_physical_values_npart( int nbPhysicalValues, int* idxOfParticles, FReal * physicalValues, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPhysicalValues && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                sources->getPhysicalValues()[idxPart] = physicalValues[iterPart];
                                checkCount++;
                                notFoundYet = false;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }else{//Parts are target
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartTarget = targets->getNbParticles();
                    //Targets part
                    for(FSize idxPart = 0 ; idxPart<nbPartTarget ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPhysicalValues && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                targets->getPhysicalValues()[idxPart] = physicalValues[iterPart];
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
        if(checkCount < nbPhysicalValues){std::cout << "Not all "<<nbPhysicalValues <<" parts has been set"<< std::endl;}
        else{
            if(checkCount > nbPhysicalValues){std::cout << "More parts than  "<<nbPhysicalValues <<" has been set"<< std::endl;}
        }
    }


    //get back the physical values
    void get_physical_values( int nbPhysicalValues, FReal * physicalValues, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        physicalValues[indexes[idxPart]] = sources->getPhysicalValues()[idxPart];
                        checkCount++;
                    }
                });
        }
        else{//Get the targets forces
            octree->forEachLeaf([&](LeafClass* leaf){
                        ContainerClass * targets = leaf->getTargets();
                        const FVector<FSize>& indexes = targets->getIndexes();
                        FSize nbPartThere = targets->getNbParticles();
                        for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                            physicalValues[indexes[idxPart]] = targets->getPhysicalValues()[idxPart];
                            checkCount++;
                        }
                    });
        }
        if(checkCount < nbPhysicalValues){std::cout << "Not all "<<nbPhysicalValues <<" parts has been read"<< std::endl;}
        else{
            if(checkCount > nbPhysicalValues){std::cout << "More parts than  "<<nbPhysicalValues <<" has been read"<< std::endl;}
        }
    }


    //Same algorithm as in set_physical_values_npart
    void get_physical_values_npart( int nbPhysicalValues, int* idxOfParticles, FReal * physicalValues, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPhysicalValues && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                physicalValues[iterPart] = sources->getPhysicalValues()[idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }else{ //Target
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPhysicalValues && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                physicalValues[iterPart] = targets->getPhysicalValues()[idxPart];
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
        if(checkCount < nbPhysicalValues){std::cout << "Not all "<<nbPhysicalValues <<" parts has been read"<< std::endl;}
        else{
            if(checkCount > nbPhysicalValues){std::cout << "More parts than  "<<nbPhysicalValues <<" has been read"<< std::endl;}
        }
    }

    void get_forces_xyz( int nbParts, FReal * forcesToFill, PartType type){
        FScalFMMEngine<FReal>::template generic_get_forces_xyz<ContainerClass,LeafClass,InterCell>(octree,nbParts,forcesToFill,type);
    }

    void get_forces(int nbParts, FReal * fX, FReal* fY, FReal* fZ, PartType type){
        FScalFMMEngine<FReal>::template generic_get_forces<ContainerClass,LeafClass,InterCell>(octree,nbParts,fX,fY,fZ,type);
    }

    void get_forces_nbpart(int nbParts, int* idxOfParticles ,FReal * fX, FReal* fY, FReal* fZ, PartType type){
        FScalFMMEngine<FReal>::template generic_get_forces_xyz_npart<ContainerClass,LeafClass,InterCell>(octree,nbParts,idxOfParticles,fX,fY,fZ,type);
    }

    void get_forces_xyz_nbpart(int nbParts, int* idxOfParticles, FReal * forcesToFill, PartType type){
        FScalFMMEngine<FReal>::template generic_get_forces_xyz_npart<ContainerClass,LeafClass,InterCell>(octree,nbParts,idxOfParticles,forcesToFill,type);
    }


    //To set initial condition
    void set_forces_xyz( int nbParts, FReal * forcesToRead, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<FSize>& indexes = sources->getIndexes();
                FSize nbPartThere = sources->getNbParticles();
                for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getForcesX()[idxPart] = forcesToRead[indexes[idxPart]*3+0];
                    sources->getForcesY()[idxPart] = forcesToRead[indexes[idxPart]*3+1];
                    sources->getForcesZ()[idxPart] = forcesToRead[indexes[idxPart]*3+2];
                    checkCount++;
                }
            });
        }
        else{//Set force on target
            octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * targets = leaf->getTargets();
                const FVector<FSize>& indexes = targets->getIndexes();
                FSize nbPartThere = targets->getNbParticles();
                for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    targets->getForcesX()[idxPart] = forcesToRead[indexes[idxPart]*3+0];
                    targets->getForcesY()[idxPart] = forcesToRead[indexes[idxPart]*3+1];
                    targets->getForcesZ()[idxPart] = forcesToRead[indexes[idxPart]*3+2];
                    checkCount++;
                }
            });
        }
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" forces has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" forces has been read"<< std::endl;}
        }
    }

    void set_forces_xyz_npart( int nbParts, int* idxOfParticles, FReal * forcesToRead, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<FSize>& indexes = sources->getIndexes();
                FSize nbPartThere = sources->getNbParticles();
                for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbParts && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            sources->getForcesX()[idxPart] = forcesToRead[iterPart];
                            sources->getForcesY()[idxPart] = forcesToRead[iterPart];
                            sources->getForcesZ()[idxPart] = forcesToRead[iterPart];
                            notFoundYet = false;
                            checkCount++;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
        }else{ //Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * targets = leaf->getTargets();
                const FVector<FSize>& indexes = targets->getIndexes();
                FSize nbPartThere = targets->getNbParticles();
                for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbParts && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            targets->getForcesX()[idxPart] = forcesToRead[iterPart];
                            targets->getForcesY()[idxPart] = forcesToRead[iterPart];
                            targets->getForcesZ()[idxPart] = forcesToRead[iterPart];
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
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" forces has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" forces has been read"<< std::endl;}
        }
    }

    void set_forces( int nbParts, FReal * fX, FReal* fY, FReal* fZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        sources->getForcesX()[idxPart] = fX[indexes[idxPart]];
                        sources->getForcesY()[idxPart] = fY[indexes[idxPart]];
                        sources->getForcesZ()[idxPart] = fZ[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        targets->getForcesX()[idxPart] = fX[indexes[idxPart]];
                        targets->getForcesY()[idxPart] = fY[indexes[idxPart]];
                        targets->getForcesZ()[idxPart] = fZ[indexes[idxPart]];
                        checkCount++;
                    }
                });
        }
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" forces has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" forces has been read"<< std::endl;}
        }
    }

    void set_forces_npart( int nbParts, int* idxOfParticles, FReal * fX, FReal* fY, FReal* fZ, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbParts && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                sources->getForcesX()[idxPart] = fX[indexes[idxPart]];
                                sources->getForcesY()[idxPart] = fY[indexes[idxPart]];
                                sources->getForcesZ()[idxPart] = fZ[indexes[idxPart]];
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
                        while(iterPart < nbParts && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                targets->getForcesX()[idxPart] = fX[indexes[idxPart]];
                                targets->getForcesY()[idxPart] = fY[indexes[idxPart]];
                                targets->getForcesZ()[idxPart] = fZ[indexes[idxPart]];
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
        if(checkCount < nbParts){std::cout << "Not all "<<nbParts <<" forces has been read"<< std::endl;}
        else{
            if(checkCount > nbParts){std::cout << "More parts than  "<<nbParts <<" forces has been read"<< std::endl;}
        }
    }



    /**
     *  Position related methods
     */
    void get_positions_xyz(int NbPositions, double * positionsToFill, PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions_xyz<ContainerClass,LeafClass,InterCell>(octree,NbPositions,positionsToFill,type);
    }
    void get_positions_xyz_npart(int NbPositions, int * idxOfParticles, double * positionsToFill,PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions_xyz_npart<ContainerClass,LeafClass,InterCell>(octree,NbPositions,idxOfParticles,positionsToFill,type);
    }
    void get_positions( int NbPositions, double *X, double *Y , double *Z, PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions<ContainerClass,LeafClass,InterCell>(octree,NbPositions,X,Y,Z,type);
    }
    void get_positions_npart(int NbPositions, int * idxOfParticles,double * X, double * Y , double * Z,PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions_npart<ContainerClass,LeafClass,InterCell>(octree,NbPositions,idxOfParticles,X,Y,Z,type);
    }
    void set_positions_xyz(int NbPositions, FReal * updatedXYZ, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions_xyz<ContainerClass,LeafClass,InterCell>(octree,NbPositions,updatedXYZ,type);
    }
    void set_positions(int NbPositions, FReal * X, FReal * Y, FReal * Z, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions<ContainerClass,LeafClass,InterCell>(octree,NbPositions,X,Y,Z,type);
    }
    void set_positions_xyz_npart(int NbPositions, int* idxOfParticles, FReal * updatedXYZ, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions_xyz_npart<ContainerClass,LeafClass,InterCell>(octree,NbPositions,idxOfParticles,updatedXYZ,type);
    }
    void set_positions_npart(int NbPositions, int* idxOfParticles, FReal * X, FReal * Y , FReal * Z, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions_npart<ContainerClass,LeafClass,InterCell>(octree,NbPositions,idxOfParticles,X,Y,Z,type);
    }
    void add_to_positions_xyz(int NbPositions,FReal * updatedXYZ,PartType type){
        FScalFMMEngine<FReal>::template generic_add_to_positions_xyz<ContainerClass,LeafClass,InterCell>(octree,NbPositions,updatedXYZ,type);
    }
    void add_to_positions(int NbPositions,FReal * X, FReal * Y , FReal * Z, PartType type){
        FScalFMMEngine<FReal>::template generic_add_to_positions<ContainerClass,LeafClass,InterCell>(octree,NbPositions,X,Y,Z,type);
    }

    //Set the potentials
    void set_potentials(int nbPotentials,FReal * potentialsToRead, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<FSize>& indexes = sources->getIndexes();
                FSize nbPartThere = sources->getNbParticles();
                for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getPotentials()[idxPart] = potentialsToRead[indexes[idxPart]];
                    checkCount++;
                }
            });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * targets = leaf->getTargets();
                const FVector<FSize>& indexes = targets->getIndexes();
                FSize nbPartThere = targets->getNbParticles();
                for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    targets->getPotentials()[idxPart] = potentialsToRead[indexes[idxPart]];
                    checkCount++;
                }
            });
        }
        if(checkCount < nbPotentials){std::cout << "Not all "<<nbPotentials <<" forces has been read"<< std::endl;}
        else{
            if(checkCount > nbPotentials){std::cout << "More parts than  "<<nbPotentials <<" forces has been read"<< std::endl;}
        }
    }

    //Set only a subpart of potentials
    //Algorithm : loop over each leaf, and then search in user array
    //if any index matches
    void set_potentials_npart( int nbPotentials, int* idxOfParticles, FReal * potentialsToRead, PartType type){
        int checkCount = 0;
        if(type == SOURCE || type==BOTH){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPotentials && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                sources->getPotentials()[idxPart] = potentialsToRead[iterPart];
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
                        while(iterPart < nbPotentials && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                targets->getPotentials()[idxPart] = potentialsToRead[iterPart];
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
        if(checkCount < nbPotentials){std::cout << "Not all "<<nbPotentials <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > nbPotentials){std::cout << "More parts than  "<<nbPotentials <<" potentials has been read"<< std::endl;}
        }
    }

    //get back the potentials
    void get_potentials( int nbPotentials, FReal * potentialsToFill, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    FSize nbPartThere = sources->getNbParticles();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        potentialsToFill[indexes[idxPart]] = sources->getPotentials()[idxPart];
                        checkCount++;
                    }
                });
        }else{//Targets
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        potentialsToFill[indexes[idxPart]] = targets->getPotentials()[idxPart];
                        checkCount++;
                    }
                });
        }
        if(checkCount < nbPotentials){std::cout << "Not all "<<nbPotentials <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > nbPotentials){std::cout << "More parts than  "<<nbPotentials <<" potentials has been read"<< std::endl;}
        }
    }

    //Same algorithm as in set_potentials_npart
    void get_potentials_npart( int nbPotentials, int* idxOfParticles, FReal * potentialsToFill, PartType type){
        int checkCount = 0;
        if(type == SOURCE){
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * sources = leaf->getSrc();
                    const FVector<FSize>& indexes = sources->getIndexes();
                    FSize nbPartThere = sources->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPotentials && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                potentialsToFill[indexes[idxPart]] = sources->getPotentials()[idxPart];
                                notFoundYet = false;
                                checkCount++;
                            }
                            else{
                                ++iterPart;
                            }
                        }
                    }
                });
        }else{
            octree->forEachLeaf([&](LeafClass* leaf){
                    ContainerClass * targets = leaf->getTargets();
                    const FVector<FSize>& indexes = targets->getIndexes();
                    FSize nbPartThere = targets->getNbParticles();
                    for(FSize idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                        int iterPart = 0;
                        bool notFoundYet = true;
                        while(iterPart < nbPotentials && notFoundYet){
                            if(indexes[idxPart] == idxOfParticles[iterPart]){
                                potentialsToFill[indexes[idxPart]] = targets->getPotentials()[idxPart];
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
        if(checkCount < nbPotentials){std::cout << "Not all "<<nbPotentials <<" potentials has been read"<< std::endl;}
        else{
            if(checkCount > nbPotentials){std::cout << "More parts than  "<<nbPotentials <<" potentials has been read"<< std::endl;}
        }
    }

    //Simple call to FScalFMMEngine method with good template
    void apply_on_cell(Callback_apply_on_cell /*not used*/){
        //We used this one to clean the Cehb cell in the user defined cheb kernel situation.
        FScalFMMEngine<FReal>::template generic_apply_on_cell<ContainerClass,InterCell,LeafClass>(octree);
    }


    // void update_tree(){
    //     if(arranger){
    //         arranger->rearrange();
    //     }
    //     else{
    //         if(FScalFMMEngine<FReal>::Algorithm == 2){ //case in wich the periodic algorithm is used
    //             arranger = new ArrangerClassPeriodic(octree);
    //             arranger->rearrange();
    //         }
    //         else{
    //             arranger = new ArrangerClass(octree);
    //             arranger->rearrange();
    //         }
    //     }
    // }


    void execute_fmm(){
        switch(FScalFMMEngine<FReal>::Algorithm){
        case 0:
            {
                typedef FFmmAlgorithm<OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassSeq;
                AlgoClassSeq* algoSeq = new AlgoClassSeq(octree,kernel);
                algoSeq->execute();
                FScalFMMEngine<FReal>::algoTimer = algoSeq;
                FScalFMMEngine<FReal>::abstrct = algoSeq;
                break;
            }
        case 1:
            {
                typedef FFmmAlgorithmThread<OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassThread;
                AlgoClassThread* algoThread = new AlgoClassThread(octree,kernel);
                algoThread->execute();
                FScalFMMEngine<FReal>::algoTimer = algoThread;
                FScalFMMEngine<FReal>::abstrct = algoThread;
                break;
            }
        case 2:
            {
                typedef FFmmAlgorithmPeriodic<FReal,OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassPeriodic;
                AlgoClassPeriodic algoPeriod(octree,2);
                algoPeriod.setKernel(kernel);
                algoPeriod.execute();
                break;
            }
        case 3:
            {
                typedef FFmmAlgorithmThreadTsm<OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassTargetSource;
                AlgoClassTargetSource* algoTS = new AlgoClassTargetSource(octree,kernel);
                algoTS->execute();
                FScalFMMEngine<FReal>::algoTimer = algoTS;
                FScalFMMEngine<FReal>::abstrct = algoTS;
                break;
            }
        default :
            std::cout<< "No algorithm found (probably for strange reasons) : "<< FScalFMMEngine<FReal>::Algorithm <<" exiting" << std::endl;
        }
    }

    void intern_dealloc_handle(Callback_free_cell unUsed){
        //this->~FInterEngine();
    }

    void print_everything(){
        octree->forEachLeaf([&](LeafClass * leaf){
                ContainerClass * sources = leaf->getSrc();
                ContainerClass * targets = leaf->getTargets();
                const FVector<FSize>& indexesSources = sources->getIndexes();
                const FVector<FSize>& indexesTargets = targets->getIndexes();
                FSize nbPartSource = sources->getNbParticles();
                FSize nbPartTarget = targets->getNbParticles();
                for(int i=0 ; i<nbPartSource ; ++i){
                    printf("Src : Leaf : %p Part : %lld/%lld, pos: %e,%e,%e phy: %e, forces: %e,%e,%e pot %e\n",
                           leaf,indexesSources[i],nbPartSource,
                           sources->getPositions()[0][i],sources->getPositions()[1][i],sources->getPositions()[2][i],
                           sources->getPhysicalValues()[i],
                           sources->getForcesX()[i],sources->getForcesY()[i],sources->getForcesZ()[i],
                           sources->getPotentials()[i]);
                }
                for(int i=0 ; i<nbPartTarget ; ++i){
                    printf("Tgt : Leaf : %p Part : %lld/%lld, pos %e,%e,%e phy: %e, forces: %e,%e,%e pot %e\n",
                           leaf,indexesTargets[i],nbPartTarget,
                           targets->getPositions()[0][i],targets->getPositions()[1][i],targets->getPositions()[2][i],
                           targets->getPhysicalValues()[i],
                           targets->getForcesX()[i],targets->getForcesY()[i],targets->getForcesZ()[i],
                           targets->getPotentials()[i]);
                }

            });
    }
};


#endif
