#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//For timing monitoring
#include "Timers.h"

#include "../Src/CScalfmmApi.h"

#include "../../Src/Kernels/Chebyshev/FChebInterface.h"


/**
 * This file is an example of the user defined kernel API, with link
 * to our Chebyshev Kernel
 **/

/**
 * @brief Wrapper to init internal ChebCell
 */
void* cheb_init_cell(int level, long long morton_index, int* tree_position, double* spatial_position, void * KernelDatas){
    return ChebCellStruct_create(morton_index,tree_position);
}

/**
 * @brief Wrapper to free internal ChebCell
 */
void cheb_free_cell(void * inCell){
    ChebCellStruct_free(inCell);
}

/**
 * No need for leaf function
 */
void * cheb_init_leaf(int level, FSize nbParts, const FSize * idxParts, long long morton_index, double center[3],
                      void * cellDatas, void * userDatas){
    //Do nothing
    int * A = malloc(sizeof(double) * nbParts);
    return A;
}

/**
 * No need for leaf function
 */
void cheb_free_leaf(void * cellDatas, FSize nbParts, const FSize * idxParts, void * leafData, void * userDatas){
    free(leafData);
    //Do nothing
}

/**
 * @brief Wrapper to FMM operators (refer to CScalfmmApi.h to get the
 * detailed descriptions)
 */
void cheb_p2m(void* cellData, void * leafData, FSize nbParticlesInLeaf, const FSize* particleIndexes,
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
void cheb_l2p(void* leafCell, void * leafData, FSize nbParticles, const FSize* particleIndexes,
              void* userData){
    ChebKernel_L2P( leafCell, nbParticles, particleIndexes, userData);
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
    ChebCell_reset(level,morton_index,tree_position,spatial_position,userCell,userData);
}

/**
 * This function is mainly a display of its args...
 */
void on_leaf(int level, FSize nbParts, const FSize * idxParts, long long morton_index, double center[3],
             void * cellDatas, void * leafData, void * userDatas){
    /* printf("I'm leaf at %lld pos, of center [%e %e %e], containing %lld parts\n", */
    /*        morton_index,center[0],center[1],center[2],nbParts); */
}



/**
 * @brief Do everything
 * @param number of particle (no default value)
 */
int main(int argc, char ** av){
    //omp_set_num_threads(1);
    printf("Start\n");
    if(argc<2){
        printf("Use : %s nb_part (optionnal : TreeHeight) \nexiting\n",av[0]);
        exit(0);
    }
    int nbPart= atoi(av[1]);;
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
    scalfmm_set_upper_limit(handle,2);
    //For Reference
    scalfmm_handle handle_ref = scalfmm_init(chebyshev,multi_thread);


    //Struct for user defined kernel
    struct User_Scalfmm_Cell_Descriptor cellDescriptor;
    cellDescriptor.user_init_cell = cheb_init_cell;
    cellDescriptor.user_free_cell = cheb_free_cell;
    cellDescriptor.user_init_leaf = cheb_init_leaf;
    cellDescriptor.user_free_leaf = cheb_free_leaf;

    //Struct for ref cheb kernel
    struct User_Scalfmm_Cell_Descriptor user_descr;
    user_descr.user_init_cell = NULL;
    user_descr.user_free_cell = NULL;


    // Init tree and cell
    printf("Building the tree and Initizalizing the cells:\n");

    scalfmm_build_tree(handle,treeHeight, boxWidth, boxCenter, cellDescriptor);
    scalfmm_build_tree(handle_ref,treeHeight, boxWidth, boxCenter, user_descr);

    //Once is the tree built, one must set the kernel before inserting particles

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

    int nb_threads = omp_get_max_threads();
    int idThreads= 0;

    //Create as many kernels as there are threads in order to void concurrent write
    userDatas.kernelStruct = ChebKernelStruct_create(treeHeight,boxWidth,boxCenter);
    /* malloc(sizeof(void *)*nb_threads); */
    /* for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){ */
    /*     userDatas.kernelStruct[idThreads] = ChebKernelStruct_create(treeHeight,boxWidth,boxCenter); // Set my kernel inside userDatas */
    /* } */

    //Only read, so no split needed
    userDatas.insertedPositions = particleXYZ;                                       // Set the position
    userDatas.myPhyValues = physicalValues;                                          // Set the physical values

    //Create as many array of forces as there are threads in order to
    //void concurrent write
    double ** forcesToStore = malloc(sizeof(double*)*nb_threads);
    //For each array, initialisation
    for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        forcesToStore[idThreads] =  malloc(sizeof(double)*nbPart*3); //allocate memory
        memset(forcesToStore[idThreads],0,sizeof(double)*nbPart*3);  //set to zero (IMPORTANT, as operators usually "+=" on forces)
    }
    userDatas.forcesComputed = forcesToStore;

    //Same idea with potential -> could be merge with loop just above,
    //but it's kept apart in order to simply
    double ** potentialToStore = malloc(sizeof(double*)*nb_threads);
    //For each array, initialisation
    for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        potentialToStore[idThreads] = malloc(sizeof(double)* nbPart);
        memset(potentialToStore[idThreads],0,sizeof(double)* nbPart);
    }
    userDatas.potentials = potentialToStore;

    //Give ScalFMM the datas before calling fmm (this will set as well the kernel)
    scalfmm_user_kernel_config(handle,kernel,&userDatas);


    // Insert particles
    printf("Inserting particles...\n");
    scalfmm_tree_insert_particles_xyz(handle, nbPart, particleXYZ,BOTH);
    scalfmm_tree_insert_particles_xyz(handle_ref, nbPart, particleXYZ,BOTH);


    //Set physical values for Cheb_ref
    scalfmm_set_physical_values(handle_ref,nbPart,physicalValues,BOTH);



    //Set timers
    Timer interface_timer,ref_timer;
    int ite=0, max_ite=1;
    while(ite<max_ite){
        //Execute FMM

        scalfmm_apply_on_leaf(handle,on_leaf);

        tic(&interface_timer);
        scalfmm_execute_fmm(handle/*, kernel, &my_data*/);
        tac(&interface_timer);

        { //Temporary
            int nbTimers = scalfmm_get_nb_timers(handle);
            double * timersArray = malloc(sizeof(double)*nbTimers);
            scalfmm_get_timers(handle,timersArray);
            int i;
            for(i=0; i<nbTimers ; ++i){
                printf("ScalFMM Operands : %d : \t %e\n",i,timersArray[i]);
            }
            free(timersArray);
        }


        //Reduction on forces & potential arrays
        {
            FSize idxPart;
            for(idThreads=1 ; idThreads<nb_threads ; ++idThreads){
                for(idxPart=0 ; idxPart<nbPart ; ++idxPart){
                    //Everything is stored in first array

                    forcesToStore[0][3*idxPart+0] += forcesToStore[idThreads][3*idxPart+0];
                    forcesToStore[0][3*idxPart+1] += forcesToStore[idThreads][3*idxPart+1];
                    forcesToStore[0][3*idxPart+2] += forcesToStore[idThreads][3*idxPart+2];
                    potentialToStore[0][idxPart] += potentialToStore[idThreads][idxPart];
                }
            }
        }

        printf("User defined Chebyshev done\n");
        print_elapsed(&interface_timer);

        tic(&ref_timer);
        scalfmm_execute_fmm(handle_ref/*, kernel, &my_data*/);
        tac(&ref_timer);

        /* { //Temporary */
        /*     int nbTimers = scalfmm_get_nb_timers(handle_ref); */
        /*     double * timersArray = malloc(sizeof(double)*nbTimers); */
        /*     scalfmm_get_timers(handle_ref,timersArray); */
        /*     int i; */
        /*     for(i=0; i<nbTimers ; ++i){ */
        /*         printf("ScalFMM Operands : %d : \t %e\n",i,timersArray[i]); */
        /*     } */
        /*     free(timersArray); */
        /* } */


        printf("Intern Chebyshev done\n");
        print_elapsed(&ref_timer);

        //Print time results
        print_difference_elapsed(&interface_timer,&ref_timer);

        //get back the forces & potentials for ref_cheb execution
        double * forcesRef     = malloc(sizeof(double)*3*nbPart);
        double * potentialsRef = malloc(sizeof(double)*nbPart);

        memset(forcesRef,0,sizeof(double)*3*nbPart);
        memset(potentialsRef,0,sizeof(double)*nbPart);

        scalfmm_get_forces_xyz(handle_ref,nbPart,forcesRef,BOTH);
        scalfmm_get_potentials(handle_ref,nbPart,potentialsRef,BOTH);
        //scalfmm_print_everything(handle_ref);

        {//Comparison part
            FSize idxPart;
            int nbPartOkay = 0;
            for(idxPart=0 ; idxPart<nbPart ; ++idxPart ){
                double diffX,diffY,diffZ,diffPot;
                diffX = fabs( forcesToStore[0][idxPart*3+0]-forcesRef[idxPart*3+0] );
                diffY = fabs( forcesToStore[0][idxPart*3+1]-forcesRef[idxPart*3+1] );
                diffZ = fabs( forcesToStore[0][idxPart*3+2]-forcesRef[idxPart*3+2] );
                diffPot = fabs( potentialToStore[0][idxPart]-potentialsRef[idxPart] );

                //THERE

                if(diffX < 0.00000001 && diffY < 0.00000001 && diffZ < 0.00000001 && diffPot < 0.00000001){
                    nbPartOkay++;
                }
                else{
                    printf("id : %lld : %e, %e, %e, %e, ChebInterf Pot : %e  Cheb Pot : %e \n",
                           idxPart,diffX,diffY,diffZ,diffPot,
                           potentialToStore[0][idxPart],
                           potentialsRef[idxPart]);
                }
                //That part is to verify with our usual exec' if everything is alright
                if(idxPart == 0 || idxPart == nbPart/2 || idxPart == nbPart-1){
                    printf("User one's id : %lld : %e, %e, %e, %e\n",idxPart,
                           forcesToStore[0][idxPart*3+0],
                           forcesToStore[0][idxPart*3+1],
                           forcesToStore[0][idxPart*3+2],
                           potentialToStore[0][idxPart]);
                    printf("Chebyshev one's id : %lld : %e, %e, %e, %e\n",idxPart,
                           forcesRef[idxPart*3+0],
                           forcesRef[idxPart*3+1],
                           forcesRef[idxPart*3+2],
                           potentialsRef[idxPart]);
                }
            }
            printf("End of simulation -- \t %d\n \t Percentage of good parts : %d/%d (%f %%) \n",ite,
                   nbPartOkay,nbPart,(((double) nbPartOkay)/(double)nbPart)*100);
        }

        free(forcesRef);
        free(potentialsRef);

        //Reset
        scalfmm_apply_on_cell(handle,cheb_resetCell);
        scalfmm_apply_on_cell(handle_ref,NULL);

        printf("Internal resets done \n");

        {//Reset User's datas
            FSize idThreads;
            for(idThreads=0;idThreads<nb_threads;++idThreads){
                memset(potentialToStore[idThreads],0,sizeof(double)*nbPart);
                memset(forcesToStore[idThreads],0,sizeof(double)*nbPart*3);
            }
        }
        printf("External resets done ...\n");

        ite++;
    }
    printf("Free the kernels\n");

    printf("Free the Handles ...\n");
    scalfmm_dealloc_handle(handle,cheb_free_cell);
    scalfmm_dealloc_handle(handle_ref,NULL);

    free(particleXYZ);
    free(physicalValues);
    //free the thread' specific datas
    for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        free(forcesToStore[idThreads]);
        free(potentialToStore[idThreads]);
    }
    free(forcesToStore);
    free(potentialToStore);

    ChebKernelStruct_free(userDatas.kernelStruct);
    return EXIT_SUCCESS;
}
