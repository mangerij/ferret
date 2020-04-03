#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//For timing monitoring
#include "Timers.h"

#include "../Src/CScalfmmApi.h"

int globalId = 0;

typedef struct cell_display{
    int Id;      //Unique Id
    int morton;  //Unique Id by level
    int coord[3];//Number of boxes from the corner of the big domain
    int lvl;     //current lvl
} Cell;

typedef struct userData{
    int nbP2M;
    int nbM2M;
    int nbM2L;
    int nbL2L;
    int nbL2P;
    int nbP2PFull;
} UserData;

/**
 * This file is an example of the user defined kernel API, with link
 * to our Chebyshev Kernel
 **/

/**
 * @brief Wrapper to init internal ChebCell
 */
void* displ_init_cell(int level, long long morton_index, int* tree_position, double* spatial_position, void * KernelDatas){
    printf("Level : %d \t morton : %lld \ntree_positions : %d %d %d,\tspatial_position : %e %e %e\n" ,
           level,morton_index,
           tree_position[0],tree_position[1],tree_position[2],
           spatial_position[0],spatial_position[1],spatial_position[2]);

    Cell * current_cell = malloc(sizeof(Cell));
    current_cell->lvl = level;
    current_cell->morton = morton_index;
    current_cell->Id = globalId++;
    //Store tree coordinate
    current_cell->coord[0] = tree_position[0];
    current_cell->coord[1] = tree_position[1];
    current_cell->coord[2] = tree_position[2];
    return current_cell;
}

/**
 * @brief Wrapper to free internal ChebCell
 */
void displ_free_cell(void * inCell){
    free((Cell * )inCell);
}

/**
 * @brief Wrapper to FMM operators (refer to CScalfmmApi.h to get the
 * detailed descriptions)
 */
void displ_p2m(void* cellData, void * leafData, FSize nbParticlesInLeaf, const FSize* particleIndexes, void* userData){
    Cell * current_cell = cellData;
    printf("P2M with %lld parts filling cell %d ::\n",nbParticlesInLeaf,current_cell->Id);

    printf("[p2m]Cell is %d/%d, level is %d, coord are %d %d %d\n\n",
           current_cell->morton,current_cell->Id,
           current_cell->lvl,
           current_cell->coord[0],current_cell->coord[1],current_cell->coord[2]);
    UserData * user_data = userData;
    user_data->nbP2M += 1;
}

void displ_m2m(int level, void* parentCell, int childPosition, void* childCell, void* userData){
    Cell * parent_cell = parentCell;
    Cell * child_cell = childCell;
    printf("M2M with %d cell and %d cell, childPosition is %d ::\n",parent_cell->Id,child_cell->Id,childPosition);
    printf("\t[m2m]Parent Cell is %d, level is %d, coord are %d %d %d, and child cell is %d, level is %d, coord are %d %d %d\n\n",
           parent_cell->morton,
           parent_cell->lvl,
           parent_cell->coord[0],parent_cell->coord[1],parent_cell->coord[2],
           child_cell->morton,
           child_cell->lvl,
           child_cell->coord[0],child_cell->coord[1],child_cell->coord[2]);
    UserData * user_data = userData;
    user_data->nbM2M += 1;
}

void displ_m2l_full(int level, void* targetCell, const int * neighborPosition, const int size, void** sourceCell, void* userData){
    Cell * target_cell = targetCell;
    printf("M2L at lvl %d filling target cell of Id %d and morton : %d with %d neighbors ::\n",target_cell->lvl,target_cell->Id,target_cell->morton,size);
    int idSource = 0;
    int countSource = 0;
    for(idSource = 0 ; idSource<size ; ++idSource){
        if(sourceCell[idSource]){
            Cell * source_cell = sourceCell[idSource];
            printf("\t[m2l]Target : %d, sourceCell[%d] : morton %d, with %d,%d,%d coordinates \n",target_cell->morton,
                   idSource,source_cell->morton,
                   source_cell->coord[0],source_cell->coord[1],source_cell->coord[2]);
            countSource++;
        }
    }
    printf("\tFinished M2L at lvl %d with filling %d cell reading %d sources\n\n",level,target_cell->morton,countSource);
    UserData * user_data = userData;
    user_data->nbM2L += 1;
}
void displ_l2l(int level, void* parentCell, int childPosition, void* childCell, void* userData){
    Cell * parent_cell = parentCell;
    Cell * child_cell = childCell;
    printf("L2L with %d cell and %d cell, childPosition is %d ::\n",parent_cell->Id,child_cell->Id,childPosition);
    printf("\t[l2l]Parent Cell is %d, level is %d, coord are %d %d %d, Child Cell is %d, level is %d, coord are %d %d %d\n\n",
           parent_cell->morton,
           parent_cell->lvl,
           parent_cell->coord[0],parent_cell->coord[1],parent_cell->coord[2],
           child_cell->morton,
           child_cell->lvl,
           child_cell->coord[0],child_cell->coord[1],child_cell->coord[2]);
    UserData * user_data = userData;
    user_data->nbL2L += 1;
}
void displ_l2p(void* leafCell,void * leafData, FSize nbParticles, const FSize* particleIndexes, void* userData){
    Cell * current_cell = leafCell;
    printf("L2P with %lld parts reading cell %d ::\n",nbParticles,current_cell->Id);

    printf("\t [l2p]Cell is %d/%d, level is %d, coord are %d %d %d\n\n",
           current_cell->morton,current_cell->Id,
           current_cell->lvl,
           current_cell->coord[0],current_cell->coord[1],current_cell->coord[2]);
    UserData * user_data = userData;
    user_data->nbL2P += 1;
}
void displ_p2pFull(void * targetLeaf, FSize nbParticles, const FSize* particleIndexes,
                   void ** sourceLeaves,
                   const FSize ** sourceParticleIndexes, FSize * sourceNbPart, const int* sourcePosition, const int size, void* userData) {
    printf("P2P, no cells involved, only leaves\n");
    UserData * user_data = userData;
    user_data->nbP2PFull+= 1;
}

void displ_resetCell(int level, long long morton_index, int* tree_position, double* spatial_position, void * userCell, void * userData){
    printf("Removing %d cell at lvl %d... \n",((Cell*)userCell)->Id,((Cell*)userCell)->lvl);
}



/**
 * @brief Do everything
 * @param number of particle (no default value)
 */
int main(int argc, char ** av){
    omp_set_num_threads(1);
    printf("Start\n");
    if(argc<2){
        printf("Use : %s nb_part (optionnal : TreeHeight) \nexiting\n",av[0]);
        exit(0);
    }
    int nbPart= atoi(av[1]);
    int treeHeight = 5 ;
    if(argc>2){
        int treeHeight = atoi(av[2]);
    }


    double boxWidth = 1.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;

    //Allocation of the positions and physical values
    double* particleXYZ = malloc(sizeof(double)*3*nbPart);
    double* physicalValues = malloc(sizeof(double)*nbPart);


    {
        printf("Creating Particles:\n");
        FSize idxPart;
        for(idxPart = 0 ; idxPart < nbPart ; ++idxPart){
            particleXYZ[idxPart*3]   = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
            particleXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
            particleXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
            physicalValues[idxPart] = 1.0;

        }

    }

    {//This part will write generated particles to a file at ScalFMM
     //format in order to verify numercal results
        /* FILE * fd = fopen("input.fma","w"); */
        /* fprintf(fd,"8\t 4\n %d\n %f\t %f\t %f\t %f\n",nbPart,boxWidth/2.0, boxCenter[0],boxCenter[1],boxCenter[2]); */
        /* FSize idxPart; */
        /* for(idxPart=0 ; idxPart<nbPart ; ++idxPart){ */
        /*     fprintf(fd,"%e\t %e\t %e\t %e \n", */
        /*             particleXYZ[idxPart*3], */
        /*             particleXYZ[idxPart*3+1], */
        /*             particleXYZ[idxPart*3+2], */
        /*             physicalValues[idxPart]); */
        /* } */
        /* fclose(fd); */
    }

    scalfmm_handle handle = scalfmm_init(user_defined_kernel,multi_thread);

    //Struct for user defined kernel
    struct User_Scalfmm_Cell_Descriptor cellDescriptor;
    cellDescriptor.user_init_cell = displ_init_cell;
    cellDescriptor.user_free_cell = displ_free_cell;


    // Init tree and cell
    printf("Building the tree and Initizalizing the cells:\n");

    scalfmm_build_tree(handle,treeHeight, boxWidth, boxCenter, cellDescriptor);

    //Once is the tree built, one must set the kernel before inserting particles

    //Set our callbacks
    struct User_Scalfmm_Kernel_Descriptor kernel;
    kernel.p2m = displ_p2m;
    kernel.m2m = displ_m2m;
    //init the other to NULL
    kernel.m2l = NULL;
    kernel.m2l_full = displ_m2l_full;
    kernel.l2l = displ_l2l;
    kernel.l2p = displ_l2p;
    kernel.p2p_full = displ_p2pFull;
    kernel.p2pinner = NULL;
    kernel.p2p = NULL;

    //Set my datas
    UserData userDatas;
    userDatas.nbP2M = 0;
    userDatas.nbM2M = 0;
    userDatas.nbM2L = 0;
    userDatas.nbL2L = 0;
    userDatas.nbL2P = 0;
    userDatas.nbP2PFull = 0;

    //Give ScalFMM the datas before calling fmm (this will set as well the kernel)
    scalfmm_user_kernel_config(handle,kernel,&userDatas);


    // Insert particles
    printf("Inserting particles...\n");
    scalfmm_tree_insert_particles_xyz(handle, nbPart, particleXYZ,BOTH);

    //Set timers
    Timer interface_timer,ref_timer;
    int ite=0, max_ite=1;
    while(ite<max_ite){
        //Execute FMM
        printf("\n\tStart FMM \n \n");
        tic(&interface_timer);
        scalfmm_execute_fmm(handle/*, kernel, &my_data*/);
        tac(&interface_timer);

        printf("User defined Chebyshev done\n");
        print_elapsed(&interface_timer);

        printf("FMM finished: \n\tNumber of P2M : %d\n\tNumber of M2M : %d\n\tNumber of M2L : %d\n\tNumber of L2L : %d\n\tNumber of L2P : %d\n\tNumber of P2P : %d\n",
               userDatas.nbP2M,
               userDatas.nbM2M,
               userDatas.nbM2L,
               userDatas.nbL2L,
               userDatas.nbL2P,
               userDatas.nbP2PFull);

        scalfmm_apply_on_cell(handle,displ_resetCell);

        ite++;
    }
    printf("Free the kernels\n");

    printf("Free the Handles ...\n");
    scalfmm_dealloc_handle(handle,displ_free_cell);

    free(particleXYZ);

    return EXIT_SUCCESS;
}
