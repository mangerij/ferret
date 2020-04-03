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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
  * In this file we implement a example kernel which is simply printing information
  * about the cells and particles.
  * It is recommanded to compile and execute this code in order to understand the behavior
  * of the application.
  * We mark some part with "JUST-PUT-HERE:" to give instruction to user to create its own kernel.
  **/

// Include the FMM API (should be in a different folder for you)
#include "../Src/CScalfmmApi.h"

// Uncomment the next line to avoid the print in the kernel and verbosity
#define NOT_TOO_MUCH_VERBOSE
#ifdef NOT_TOO_MUCH_VERBOSE
#define VerbosePrint(X...)
#else
#define VerbosePrint(X...) printf(X)
#endif

///////////////////////////////////////////////////////////////////////////
/// Cell struct and functions
///////////////////////////////////////////////////////////////////////////

// This represent a cell
struct MyCellDescriptor{
    long long int dataUp,dataDown;
    int level;
    long long mortonIndex;
    int coord[3];
    double position[3];
    // JUST-PUT-HERE:
    // You local and multipole arrays
};

// This is our function that init a cell (struct MyCellDescriptor)
void* my_Callback_init_cell(int level, long long mortonIndex, int* coord, double* position, void * kernel){
    VerbosePrint("\tAllocating cell for level %d, morton index %lld, coord %d/%d/%d\n", level, mortonIndex, coord[0], coord[1], coord[2]);
    struct MyCellDescriptor* cellData = (struct MyCellDescriptor*)malloc(sizeof(struct MyCellDescriptor));
    memset(cellData, 0, sizeof(struct MyCellDescriptor));
    cellData->level       = level;
    cellData->mortonIndex = mortonIndex;

    //Count the interactions
    cellData->dataUp = 0;
    cellData->dataDown = 0;

    memcpy(cellData->coord, coord, sizeof(int)*3);
    memcpy(cellData->position, position, sizeof(double)*3);
    // JUST-PUT-HERE:
    // Fill your structure
    return cellData;
}

// To dealloc a cell
void my_Callback_free_cell(void* cellData){
    free(cellData);
}


///////////////////////////////////////////////////////////////////////////
/// Kernel struct and functions
///////////////////////////////////////////////////////////////////////////

// This is the data passed to our kernel
struct MyData {
    // We simply strore the number of call the and particle indexes
    double* insertedPositions;
    int countP2M;
    int countM2M;
    int countM2L;
    int countL2L;
    int countL2P;
    int countP2PInner;
    int countP2P;
    int countM2L_ext;
    // JUST-PUT-HERE:
    // everything your kernel will need for its computation
    // pre-computed matrices, etc....
    long long int * arrayOfInfluence;
};


// Our P2M
void my_Callback_P2M(void* cellData, void * leafData, FSize nbParticlesInLeaf, const FSize* particleIndexes, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countP2M += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    my_cell->dataUp += nbParticlesInLeaf;

    VerbosePrint("Cell morton %lld is doing P2M with %lld particles :\n", my_cell->mortonIndex, nbParticlesInLeaf);
    FSize idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        VerbosePrint("\t particle idx %lld position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
        // JUST-PUT-HERE:
        // Your real P2M computation
    }
}

void my_Callback_M2M(int level, void* cellData, int childPosition, void* childData, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countM2M += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    struct MyCellDescriptor* my_child = (struct MyCellDescriptor*) childData;

    my_cell->dataUp += my_child->dataUp;

    int childFullPosition[3];
    scalfmm_utils_parentChildPosition(childPosition, childFullPosition);

    VerbosePrint("Doing a M2M at level %d, between cells %lld and %lld (position %d/%d/%d)\n",
           level, my_cell->mortonIndex, my_child->mortonIndex,
           childFullPosition[0], childFullPosition[1], childFullPosition[2]);
    // JUST-PUT-HERE: your real M2M computation
}

void my_Callback_M2L(int level, void* cellData, int srcPosition, void* srcData, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countM2L += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    struct MyCellDescriptor* my_src_cell = (struct MyCellDescriptor*) srcData;

    my_cell->dataDown += my_src_cell->dataUp;

    int interactionFullPosition[3];
    scalfmm_utils_interactionPosition(srcPosition, interactionFullPosition);

    VerbosePrint("Doing a M2L at level %d, between cells %lld and %lld (position %d/%d/%d)\n",
           level, my_cell->mortonIndex, my_src_cell->mortonIndex,
           interactionFullPosition[0], interactionFullPosition[1], interactionFullPosition[2]);
    // JUST-PUT-HERE: Your M2L
}

void my_Callback_M2L_ext(int level, void * cellTgt, void * cellSrc, int transfer[3], void * userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countM2L_ext += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellTgt;
    struct MyCellDescriptor* my_src_cell = (struct MyCellDescriptor*) cellSrc;

    my_cell->dataDown += my_src_cell->dataUp;

    VerbosePrint("Doing a M2L_ext at level %d, between cells %lld and %lld (transfer %d/%d/%d)\n",
                 level, my_cell->mortonIndex, my_src_cell->mortonIndex,
                 transfer[0], transfer[1], transfer[2]);
}


void my_Callback_L2L(int level, void* cellData, int childPosition, void* childData, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countL2L += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    struct MyCellDescriptor* my_child = (struct MyCellDescriptor*) childData;

    my_child->dataDown += my_cell->dataDown;

    int childFullPosition[3];
    scalfmm_utils_parentChildPosition(childPosition, childFullPosition);

    VerbosePrint("Doing a L2L at level %d, between cells %lld and %lld (position %d/%d/%d)\n",
           level, my_cell->mortonIndex, my_child->mortonIndex,
           childFullPosition[0], childFullPosition[1], childFullPosition[2]);
    // JUST-PUT-HERE: Your L2L
}

void my_Callback_L2P(void* cellData, void * leafData, FSize nbParticlesInLeaf, const FSize* particleIndexes, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countL2P += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;

    //printf("[%lld] Far Field Influence : \t%lld parts have been summed on %lld locals parts \n",my_cell->mortonIndex,my_cell->dataDown,nbParticlesInLeaf);

    VerbosePrint("Cell morton %lld is doing L2P with %lld particles :\n", my_cell->mortonIndex, nbParticlesInLeaf);
    FSize idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        //Store the number of cells that contribute to the far field of current cell.
        my_data->arrayOfInfluence[particleIndexes[idxPart]] += my_cell->dataDown;
        VerbosePrint("\t particle idx %lld position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
        // JUST-PUT-HERE: Your L2P
    }
}

void my_Callback_P2P(void * targetLeaf, FSize nbParticlesInLeaf, const FSize* particleIndexes,void * sourceLeaf,
                     FSize nbParticlesInSrc, const FSize* particleIndexesSrc, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countP2P += 1;

    VerbosePrint("Doing P2P between two leaves of %lld particles and %lld particles :\n", nbParticlesInLeaf, nbParticlesInSrc);
    FSize idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];

        VerbosePrint("\t Target >> particle idx %lld position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
    }
    for(idxPart = 0 ; idxPart < nbParticlesInSrc ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexesSrc[idxPart]*3];
        VerbosePrint("\t Target >> Src idx %lld position %e/%e%e\n", particleIndexesSrc[idxPart],
               position[0], position[1], position[2]);
    }

    for(idxPart = 0; idxPart < nbParticlesInLeaf ; ++idxPart){
        my_data->arrayOfInfluence[particleIndexes[idxPart]] += nbParticlesInSrc;
    }
    // JUST-PUT-HERE:
    // Put one loop into the other to make all particles from source
    // interacting with the target particles
}

void my_Callback_P2PInner(void * targetLeaf, FSize nbParticlesInLeaf, const FSize* particleIndexes, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countP2PInner += 1;

    VerbosePrint("Doing P2P inside a leaf of %lld particles :\n", nbParticlesInLeaf);
    FSize idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        my_data->arrayOfInfluence[particleIndexes[idxPart]] += nbParticlesInLeaf-1;
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        VerbosePrint("\t particle idx %lld position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
        // JUST-PUT-HERE: another loop to have all particles interacting with
        // the other from the same leaf
    }
}

///////////////////////////////////////////////////////////////////////////
/// Main
///////////////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    omp_set_num_threads(1);
    // The properties of our tree
    int treeHeight = 5;
    double boxWidth = 1.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;

    // Create random particles
    FSize nbParticles = 100;
    int particleIndexes[nbParticles];
    double * particleXYZ = malloc(sizeof(double)*nbParticles*3);
    {
        printf("Creating Particles:\n");
        FSize idxPart;
        for(idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            particleIndexes[idxPart] = idxPart;
            particleXYZ[idxPart*3]   = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
            particleXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
            particleXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
            VerbosePrint("\t %lld] %e/%e/%e\n", idxPart, particleXYZ[idxPart*3], particleXYZ[idxPart*3+1], particleXYZ[idxPart*3+2]);
        }
    }

    // Init the handle
    scalfmm_handle handle = scalfmm_init(user_defined_kernel,multi_thread);

    //Build our own call backs
    struct User_Scalfmm_Cell_Descriptor cellDescriptor;
    cellDescriptor.user_init_cell = my_Callback_init_cell;
    cellDescriptor.user_free_cell = my_Callback_free_cell;
    // Init tree and cell
    printf("Building the tree and Initizalizing the cells:\n");

    scalfmm_build_tree(handle,treeHeight, boxWidth, boxCenter, cellDescriptor);
    // Insert particles
    printf("Inserting particles...\n");
    scalfmm_tree_insert_particles_xyz(handle, nbParticles, particleXYZ,BOTH);
    printf("Particles Inserted ...\n");

    // Init our callback struct
    struct User_Scalfmm_Kernel_Descriptor kernel;
    kernel.p2m = my_Callback_P2M;
    kernel.m2m = my_Callback_M2M;
    kernel.m2l = my_Callback_M2L;
    kernel.m2l_ext = my_Callback_M2L_ext;
    kernel.m2l_full = NULL;
    kernel.l2l = my_Callback_L2L;
    kernel.l2p = my_Callback_L2P;
    kernel.p2pinner = my_Callback_P2PInner;
    kernel.p2p = my_Callback_P2P;
    kernel.p2p_full = NULL;

    // Init the data to pass to all our callbacks
    struct MyData my_data;
    memset(&my_data, 0, sizeof(struct MyData));
    my_data.insertedPositions = particleXYZ;
    //Set my datas before calling fmm (this will set as well the kernel)
    my_data.arrayOfInfluence = malloc(sizeof(long long int)* nbParticles);
    memset(my_data.arrayOfInfluence,0,sizeof(long long int)* nbParticles);

    scalfmm_user_kernel_config(handle,kernel,&my_data);

    printf("Kernel set ... \n");


    //loop to multiples runs of the fmm
    int nb_ite = 1;
    int curr_ite = 0;
    //array to store positions
    double new_positions[nbParticles*3];
    memset(new_positions,0,3*nbParticles*sizeof(double));
    scalfmm_set_upper_limit(handle,4);

    while(curr_ite < nb_ite){
        printf("Start FMM number %d/%d ... \n", curr_ite,nb_ite);
        // Execute the FMM
        scalfmm_execute_fmm(handle/*, kernel, &my_data*/);
        printf("FMM finished ... \n");
        scalfmm_get_positions_xyz(handle,nbParticles,new_positions,BOTH);
        //Computation on those positions
        //Here it's a random
        int id;
        for(id = 0 ; id < nbParticles ; ++id){
            new_positions[id*3  ] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
            new_positions[id*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
            new_positions[id*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
        }
        printf("Positions changed \n");

        //scalfmm_set_positions_xyz(handle,nbParticles,new_positions,BOTH);
        //scalfmm_update_tree(handle);
        curr_ite++;
    }

    //Check the result:
    int idPart;

    for(idPart = 0; idPart < nbParticles ; ++idPart){
        if(my_data.arrayOfInfluence[idPart] != nbParticles-1){
            printf("Probleme with part %d \t:: %lld\n",idPart,my_data.arrayOfInfluence[idPart]);
        }
    }

    // Print the results store in our callback
    printf("There was %d \tP2M\n", my_data.countP2M);
    printf("There was %d \tM2M\n", my_data.countM2M);
    printf("There was %d \tM2L\n", my_data.countM2L);
    printf("There was %d \tM2L_ext\n", my_data.countM2L_ext);
    printf("There was %d \tL2L\n", my_data.countL2L);
    printf("There was %d \tL2P\n", my_data.countL2P);
    printf("There was %d \tP2PInner\n", my_data.countP2PInner);
    printf("There was %d \tP2P\n", my_data.countP2P);

    free(particleXYZ);
    free(my_data.arrayOfInfluence);
    // Dealloc the handle
    scalfmm_dealloc_handle(handle,my_Callback_free_cell);

    return 0;
}
