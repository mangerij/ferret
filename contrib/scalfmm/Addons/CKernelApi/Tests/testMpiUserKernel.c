#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//For timing monitoring
#include "Timers.h"
#include <mpi.h>
#include "../Src/CScalfmmApi.h"
#include "../../Src/Kernels/Chebyshev/FChebInterface.h"


// ==== CMAKE =====
// @FUSE_MPI
// ================

/**
 * This file is an example of the user defined kernel API, with link
 * to our Chebyshev Kernel
 **/

/**
 * @brief Wrapper to init internal ChebCell
 */
void* cheb_init_cell(int level, long long morton_index, int* tree_position,
                     double* spatial_position, void * KernelDatas){
    return ChebCellStruct_create(morton_index,tree_position);
}

/**
 * @brief Wrapper to free internal ChebCell
 */
void cheb_free_cell(void * inCell){
    ChebCellStruct_free(inCell);
}

/**
 * @brief Wrapper to FMM operators (refer to CScalfmmApi.h to get the
 * detailed descriptions)
 */
void cheb_p2m(void* cellData,void * leafData, FSize nbParticlesInLeaf, const FSize* particleIndexes,
              void* userData){
    ChebKernel_P2M(cellData,nbParticlesInLeaf,particleIndexes,userData);
}
void cheb_m2m(int level, void* parentCell, int childPosition, void* childCell,
              void* userData){
    ChebKernel_M2M(level,parentCell,childPosition,childCell,userData);
}
void cheb_m2l_full(int level, void* targetCell,const int* neighborPosition, const int size, void** sourceCell,
                   void* userData){
    ChebKernel_M2L(level, targetCell, neighborPosition, size, sourceCell, userData);
}
void cheb_l2l(int level, void* parentCell, int childPosition, void* childCell,
              void* userData){
    ChebKernel_L2L( level, parentCell, childPosition, childCell,  userData);
}
void cheb_l2p(void* cellData, void* leafData, FSize nbParticles, const FSize* particleIndexes,
              void* userData){
    ChebKernel_L2P( cellData, nbParticles, particleIndexes, userData);
}
void cheb_p2pFull(void * targetLeaf, FSize nbParticles, const FSize* particleIndexes,
                  void ** sourceLeaves,
                  const FSize ** sourceParticleIndexes, FSize* sourceNbPart,const int * sourcePosition,
                  const int size, void* userData) {
    ChebKernel_P2P(nbParticles, particleIndexes, sourceParticleIndexes, sourceNbPart,sourcePosition,size,
                   userData);
}

void cheb_resetCell(int level, long long morton_index, int* tree_position,
                    double* spatial_position, void * userCell, void * userData){
    ChebCell_reset(level,morton_index,tree_position,
                   spatial_position,userCell,userData);
}

FSize cheb_get_size(int level,void * userData, long long morton_index){
    return ChebCell_getSize(level,userData,morton_index);
}

void cheb_copy_cell(void * userDatas, FSize size, void * memoryAllocated){
    ChebCell_copy(userDatas,size,memoryAllocated);
}

void * cheb_restore_cell(int level, void * arrayTobeRead){
    return ChebCell_restore(level,arrayTobeRead);
}

void on_leaf(int level, FSize nbParts, const FSize * idxParts, long long morton_index, double center[3],
             void * cellDatas,void * leafdata, void * userDatas){

    UserData * ptrToUserData = (UserData*) userDatas;
    //Compute totalEnergy
    int nbThread = omp_get_max_threads();
    int i,j;
    for( i=0 ; i<nbParts ; ++i){
        double pot = 0.0;
        //Small loop over the different arrays
        for( j=0 ; j<nbThread ; ++j){
            pot += ptrToUserData->potentials[j][idxParts[i]];
        }
        ptrToUserData->totalEnergy += pot*(ptrToUserData->myPhyValues[idxParts[i]]);
    }

    /* printf("I'm leaf at %lld pos, of center [%e %e %e], containing %lld parts\n", */
    /*        morton_index,center[0],center[1],center[2],nbParts); */

    /* printf("\n"); */
    /* for(i=0 ; i<nbParts ; ++i){ */
    /*     double fx=0,fy=0,fz=0; */
    /*     for( j=0 ; j<nbThread ; ++j){ */
    /*         fx = ptrToUserData->forcesComputed[j][idxParts[i]*3+0]; */
    /*         fy = ptrToUserData->forcesComputed[j][idxParts[i]*3+1]; */
    /*         fz = ptrToUserData->forcesComputed[j][idxParts[i]*3+2]; */
    /*     } */
    /*     printf("Parts : id: %lld \t %e - %e - %e\n",idxParts[i],fx,fy,fz); */
    /* } */
}

int main(int argc, char ** argv){
    omp_set_num_threads(1);
    if(argc<2){
        printf("Use : %s nb_part (optionnal : TreeHeight) \nexiting\n",argv[0]);
        exit(0);
    }

    int nbPart= atoi(argv[1]);
    int treeHeight = 5 ;
    if(argc>2){
        int treeHeight = atoi(argv[2]);
    }
    double boxWidth = 1.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;
    //Allocation of the locals positions and physical values
    double* particleXYZ = malloc(sizeof(double)*3*nbPart);
    double* physicalValues = malloc(sizeof(double)*nbPart);

    //Initialisation
    int provided;
    //Init MPI
    MPI_Init_thread(&argc,&argv,MPI_THREAD_SERIALIZED,&provided);

    //Classic MPI values
    int my_rank;
    int nbProc;
    int nameLen;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nbProc );
    MPI_Get_processor_name(processorName,&nameLen);

    //For initialising seeds
    int seed = (unsigned)time(NULL) + my_rank*nbProc + nameLen;
    srand(seed);

    {//Create bunch of positions
        printf("Creating Particles:\n");
        FSize idxPart = 0;
        for(idxPart = 0 ; idxPart < nbPart ; ++idxPart){
            particleXYZ[idxPart*3]   = (random()/(double)(RAND_MAX))*boxWidth
                - boxWidth/2 + boxCenter[0];
            particleXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth
                - boxWidth/2 + boxCenter[1];
            particleXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth
                - boxWidth/2 + boxCenter[2];
            physicalValues[idxPart] = 1.0;
        }
    }

    //Init scalfmm
    scalfmm_handle Handle = scalfmm_init_distributed(user_defined_kernel,mpi,
                                                     MPI_COMM_WORLD); //Only difference
    printf("Handle :: %p\n",Handle);
    //This is the cell descriptor
    struct User_Scalfmm_Cell_Descriptor cellDescriptor;
    cellDescriptor.user_init_cell = cheb_init_cell;
    cellDescriptor.user_free_cell = cheb_free_cell;
    cellDescriptor.user_get_size = cheb_get_size;
    cellDescriptor.user_copy_cell = cheb_copy_cell;
    cellDescriptor.user_restore_cell = cheb_restore_cell;

    //Set our callbacks
    struct User_Scalfmm_Kernel_Descriptor kernel;
    kernel.p2m = cheb_p2m;
    kernel.m2m = cheb_m2m;
    //init the other to NULL
    kernel.m2l = NULL;
    kernel.m2l_full = cheb_m2l_full;
    kernel.l2l = cheb_l2l;
    kernel.l2p = cheb_l2p;
    kernel.p2p_sym = NULL;
    kernel.p2p_full = cheb_p2pFull;
    kernel.p2pinner = NULL;
    kernel.p2p = NULL;

    //Set my datas
    UserData userDatas;

    //Give ScalFMM the datas before calling fmm (this will set as well the kernel)
    scalfmm_user_kernel_config(Handle,kernel,&userDatas);

    scalfmm_build_tree(Handle,treeHeight, boxWidth, boxCenter, cellDescriptor);

    printf("Tree built.\n");
    double *  outputArray;     //Will be allocated inside create_local_partition
    double ** outputArrayPtr = &outputArray;
    //In order to store the indexes
    FSize *  outputIndexes;
    FSize ** outputIndexesPtr = &outputIndexes;
    //In order to store the other attributes
    double * outputPhyVal;
    double ** outputPhyValPtr = &outputPhyVal;
    //Get the number of pts
    FSize outputNbPoint;
    FSize * outputNbPointPtr = &outputNbPoint;

    scalfmm_create_local_partition(Handle,nbPart,particleXYZ,outputArrayPtr,outputIndexesPtr,outputNbPointPtr);

    printf("Partinionning done, I'm process %d/%d, I had %d, now i have %lld \n",
           my_rank,nbProc,nbPart,outputNbPoint);


    //Here, we store what we need inside the userData, but we do that
    //accordingly to the partition given by scalfmm
    int nb_threads = omp_get_max_threads();
    int idThreads= 0;

    userDatas.kernelStruct = ChebKernelStruct_create(treeHeight,boxWidth,boxCenter);
    //Only read, so no split needed
    userDatas.insertedPositions = outputArray;                     // Set the position

    //Need to get the good ones...
    double * newPhyValues = malloc(sizeof(double) * outputNbPoint);
    //In this part, we store the physicalvalues

    //Create as many array of forces as there are threads in order to
    //void concurrent write
    double ** forcesToStore = malloc(sizeof(double*)*nb_threads);
    double ** potentialToStore = malloc(sizeof(double*)*nb_threads);
    //For each array, initialisation
    for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        forcesToStore[idThreads] =  malloc(sizeof(double)*outputNbPoint*3); //allocate memory
        memset(forcesToStore[idThreads],0,sizeof(double)*outputNbPoint*3);  //set to zero (IMPORTANT, as operators usually "+=" on forces)
        potentialToStore[idThreads] = malloc(sizeof(double)* outputNbPoint);
        memset(potentialToStore[idThreads],0,sizeof(double)* outputNbPoint);
    }
    userDatas.forcesComputed = forcesToStore;
    userDatas.potentials = potentialToStore;
    userDatas.totalEnergy = 0;

    scalfmm_tree_insert_particles_xyz(Handle,outputNbPoint,outputArray,BOTH);

    printf("Insertion done, I'm process %d, and have inserted %lld parts.\n",my_rank,outputNbPoint);

    scalfmm_generic_partition(Handle,nbPart,sizeof(physicalValues[0]),physicalValues,outputPhyValPtr);

    printf("[%d] Generic Partition Done ! \n",my_rank);

    userDatas.myPhyValues = outputPhyVal;

    //Execution !!
    scalfmm_execute_fmm(Handle);

    scalfmm_apply_on_leaf(Handle,on_leaf);

    printf("Total energy for proc %d is \t%e\n",my_rank,userDatas.totalEnergy);

    printf("Execution done, I'm process %d.\n",my_rank);

    //Dealloc scalfmm handle
    scalfmm_dealloc_handle(Handle,cheb_free_cell);

    {//This part will write generated particles to a file at ScalFMM
     //format in order to verify numercal results
        /* if(my_rank==0){ */
        /*     FILE * fd = fopen("input.fma","a+"); */
        /*     fprintf(fd,"8\t 4\n %d\n %f\t %f\t %f\t %f\n",nbPart,boxWidth/2.0, boxCenter[0],boxCenter[1],boxCenter[2]); */
        /*     FSize idxPart; */
        /*     for(idxPart=0 ; idxPart<nbPart ; ++idxPart){ */
        /*         fprintf(fd,"%e\t %e\t %e\t %e \n", */
        /*                 particleXYZ[idxPart*3], */
        /*                 particleXYZ[idxPart*3+1], */
        /*                 particleXYZ[idxPart*3+2], */
        /*                 physicalValues[idxPart]); */
        /*     } */
        /*     fclose(fd); */
        /*     MPI_Barrier(MPI_COMM_WORLD); */
        /* }else{ */
        /*     MPI_Barrier(MPI_COMM_WORLD); */
        /*     FILE * fd = fopen("input.fma","a+"); */
        /*     fprintf(fd,"8\t 4\n %d\n %f\t %f\t %f\t %f\n",nbPart,boxWidth/2.0, boxCenter[0],boxCenter[1],boxCenter[2]); */
        /*     FSize idxPart; */
        /*     for(idxPart=0 ; idxPart<nbPart ; ++idxPart){ */
        /*         fprintf(fd,"%e\t %e\t %e\t %e \n", */
        /*                 particleXYZ[idxPart*3], */
        /*                 particleXYZ[idxPart*3+1], */
        /*                 particleXYZ[idxPart*3+2], */
        /*                 physicalValues[idxPart]); */
        /*     } */
        /*     fclose(fd); */
        /* } */
    }

    free(outputIndexes);
    free(outputArray);


    MPI_Finalize();
    return 0;
}
