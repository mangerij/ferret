#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../Src/CScalfmmApi.h"



void compute_displacement_from_forces(double * fXYZ, double * charge,double* dXYZ, int nb_xyz){
    int idPart;
    for(idPart = 0 ; idPart < nb_xyz ; ++idPart){
        dXYZ[idPart*3+0] = fXYZ[idPart*3+0]*0.001/charge[idPart];
        dXYZ[idPart*3+1] = fXYZ[idPart*3+1]*0.001/charge[idPart];
        dXYZ[idPart*3+2] = fXYZ[idPart*3+2]*0.001/charge[idPart];
    }
}



int main(int argc, char ** av){

    scalfmm_kernel_type myChoice = chebyshev;

    //Octree configuration
    int TreeHeight = 7;
    double boxWidth = 2.0;
    double boxCenter[3] = {0.0,0.0,0.0};

    //Init our lib
    scalfmm_handle handle = scalfmm_init(/* TreeHeight,boxWidth,boxCenter, */myChoice,multi_thread); //The tree is built
    struct User_Scalfmm_Cell_Descriptor user_descr;
    user_descr.user_init_cell = NULL;
    user_descr.user_free_cell = NULL;

    scalfmm_build_tree(handle,TreeHeight,boxWidth,boxCenter,user_descr);
    scalfmm_algorithm_config(handle,periodic);
    //Creation of an array of particles
    int nb_of_parts = 2;
    FSize idxPart;
    double * positionsXYZ = malloc(sizeof(double)*3*nb_of_parts);
    memset(positionsXYZ,0,sizeof(double)*3*nb_of_parts);

    for(idxPart = 0; idxPart<nb_of_parts ; ++idxPart){
        positionsXYZ[idxPart*3+0] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
        positionsXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
        positionsXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
    }

    //Creation of charge for each part
    double * array_of_charge = malloc(sizeof(double)*nb_of_parts);
    for(idxPart = 0; idxPart<nb_of_parts ; ++idxPart){
        array_of_charge[idxPart] = 1.0;
    }

    //Inserting the array in the tree
    scalfmm_tree_insert_particles_xyz(handle,nb_of_parts,positionsXYZ,BOTH);

    //Set the charge
    scalfmm_set_physical_values(handle,nb_of_parts,array_of_charge,BOTH);


    //Computation Part
    int nb_iteration = 3;//atoi(av[1]);
    int curr_iteration = 0;

    //Array to store the output
    double * array_of_forces = malloc(sizeof(double)*3*nb_of_parts);
    memset(array_of_forces,0,sizeof(double)*3*nb_of_parts);
    double * array_of_pot = malloc(sizeof(double)*nb_of_parts);
    memset(array_of_pot,0,sizeof(double)*nb_of_parts);

    //Array to store the displacement
    double * array_of_displacement = malloc(sizeof(double)*3*nb_of_parts);
    memset(array_of_displacement,0,sizeof(double)*3*nb_of_parts);

    //Start of an iteration loop
    while(curr_iteration < nb_iteration){
        //Execute
        scalfmm_execute_fmm(handle);

        //Get the resulting forces
        scalfmm_get_forces_xyz(handle,nb_of_parts,array_of_forces,BOTH);
        //Get the resulting potential
        scalfmm_get_potentials(handle,nb_of_parts,array_of_pot,BOTH);



        //User function to compute the movement of each part
        compute_displacement_from_forces(array_of_forces,array_of_charge,array_of_displacement,nb_of_parts);

        //get position in order to display
        scalfmm_get_positions_xyz(handle,nb_of_parts,positionsXYZ,BOTH);

        //Display forces :
        {
            printf("Iteration n° %d\n",curr_iteration);
            for(idxPart = 0 ; idxPart< nb_of_parts ; ++idxPart){
                printf("Pos : %e, %e, %e, Fxyz %e %e %e, \n \t displ : %e, %e, %e \n",
                       positionsXYZ[idxPart*3+0],
                       positionsXYZ[idxPart*3+1],
                       positionsXYZ[idxPart*3+2],
                       array_of_forces[idxPart*3],
                       array_of_forces[idxPart*3+1],
                       array_of_forces[idxPart*3+2],
                       array_of_displacement[idxPart*3+0],
                       array_of_displacement[idxPart*3+1],
                       array_of_displacement[idxPart*3+2]);
            }
        }
        //Apply displacement computed
        scalfmm_add_to_positions_xyz(handle,nb_of_parts,array_of_displacement,BOTH);

        //Update Consequently the tree
        scalfmm_update_tree(handle);
        //Continue the loop
        curr_iteration++;
    }

    //Free memory
    free(positionsXYZ);
    free(array_of_charge);
    free(array_of_forces);
    free(array_of_pot);
    free(array_of_displacement);

    //End of Computation, useer get the position after some iterations
    scalfmm_dealloc_handle(handle,NULL);


    return 0;
}
