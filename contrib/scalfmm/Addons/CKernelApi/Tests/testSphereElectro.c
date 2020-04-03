#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>




//For timing monitoring
#include "Timers.h"

#include "../Src/CScalfmmApi.h"

#include "../../Src/Kernels/Chebyshev/FChebInterface.h"

double getRandom(){
    return (random()/(double)(RAND_MAX));
}

void generateSurfacePointOnUnitSphere(int N, double * points){
    double u, v, theta, phi, sinPhi ;
    //
    int j = 0,i=0 ;
    for ( i = 0 ; i< N ; ++i, j+=3)  {
        //
        u = getRandom() ;  v = getRandom() ;
        theta  = 2*M_PI*u ;
        phi     = acos(2*v-1);
        sinPhi = sin(phi);
        //
        points[j]   =     cos(theta)*sinPhi ;
        points[j+1] =   sin(theta)*sinPhi ;
        points[j+2] =   2*v-1 ;
        //
    }
}

void generateSurfacePoints(double rayon, double* centre, int nbDePoints, double* points){

    generateSurfacePointOnUnitSphere(nbDePoints , points) ;
    int j =0,i=0 ;
    for ( i = 0 ; i< nbDePoints ; ++i, j+=3)  {
        points[j]    *= rayon + centre[0];
        points[j+1]  *= rayon + centre[1];
        points[j+2]  *= rayon + centre[2];
    }
}

void generateInsidePoints(double rayon, double*centre, int nbDePoints,  double* points){
    generateSurfacePointOnUnitSphere(nbDePoints, points);
    int j=0;
    double u;
    for(j=0 ; j<nbDePoints ; ++j){
        u = getRandom();
        points[j]    *= (rayon + centre[0])*u;
        points[j+1]  *= (rayon + centre[1])*u;
        points[j+2]  *= (rayon + centre[2])*u;
    }
}

void displayPoints(int nbPoints, double * points){
    int i = 0;
    for(i=0 ; i<nbPoints ; ++i){
        printf("%e %e %e \n",points[i*3],points[i*3+1],points[i*3+2]);
    }
}

void displayArray(int nbValue, double * array){
    int i = 0;
    for(i=0 ; i<nbValue ; ++i){
        printf("%e \n",array[i]);
    }
}

void getNormal(double * positions, double * normeToFill){
    int i;
    double norme = sqrt(positions[0]*positions[0] + positions[1]*positions[1] + positions[2]*positions[2]);
    for(i=0 ; i<3 ; ++i){
        normeToFill[i] = positions[i]/norme;
    }
    /* printf("Tgt Norme %e - %e - %e\n", */
    /*        normeToFill[0], */
    /*        normeToFill[1], */
    /*        normeToFill[2]); */
}

void computeNormalXForces(int nbPoints, double * forcesToRead, double * positionsToRead, double * arrayToFill){
    double * currentNormal = malloc(sizeof(double)*3);
    int idxPart,i;
    for(idxPart = 0 ; idxPart<nbPoints ; ++idxPart){
        getNormal(&positionsToRead[idxPart],currentNormal); //get the norme
        for(i=0 ; i<3 ; ++i){
            arrayToFill[idxPart] += currentNormal[i]*forcesToRead[idxPart+i];
        }
    }
    free(currentNormal);
}

int main(int argc, char ** av){
    omp_set_num_threads(4);
    printf("Start %s nb_targets nb_sources \n",av[0]);
    if(argc<3){
        printf("Use : %s nb_part(cible) nb_part(source) (optionnal : TreeHeight) \nexiting\n",av[0]);
        exit(0);
    }
    int nbPartTarget= atoi(av[1]);
    int treeHeight = 5 ;
    if(argc>3){
        int treeHeight = atoi(av[3]);
        printf("Tree Heigth Input %d\n",treeHeight );
    }

    double boxWidth = 2.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;

    int i;
    //Allocation of the target points
    double * targetsXYZ = malloc(sizeof(double)* 3*nbPartTarget);
    double * targetsPhiValues = malloc(sizeof(double)* nbPartTarget);
    //Memset (au cas ou)
    memset(targetsXYZ,0,sizeof(double)*3*nbPartTarget);
    memset(targetsPhiValues,0,sizeof(double)*nbPartTarget);
    //Fill
    for(i=0 ; i<nbPartTarget ; ++i){
        targetsPhiValues[i] = 1.0;
    }

    generateSurfacePoints(1.0,boxCenter,nbPartTarget,targetsXYZ);
    printf("%d Surface points generated : Targets\n",nbPartTarget);

    //Allocation of the sources points
    int nbPartSource = atoi(av[2]);
    double * sourceXYZ = malloc(sizeof(double)* 3*nbPartSource);
    double * sourcePhiValues = malloc(sizeof(double)* nbPartSource);
    //Set to Zero
    memset(sourceXYZ,0,3*sizeof(double)*nbPartSource);
    memset(sourcePhiValues,0,sizeof(double)*nbPartSource);
    //Fill
    for(i=0 ; i<nbPartSource ; ++i){
        sourcePhiValues[i] = 1.0;
    }

    generateInsidePoints(1.0,boxCenter,nbPartSource,sourceXYZ);
    //displayPoints(nbPartTarget,targetsXYZ);

    printf("%d Inside points generated :  Sources\n",nbPartSource);

    //displayPoints(nbPartSource,sourceXYZ);
    //Creation of arrays to store forces
    double * arrayOfForces = malloc(sizeof(double )* 3 * nbPartTarget);
    double * arrayOfPot    = malloc(sizeof(double ) * (nbPartTarget));
    memset(arrayOfForces,0,sizeof(double)* 3 * (nbPartTarget));
    memset(arrayOfPot,0,sizeof(double)* (nbPartTarget));

    {//Start of computation

        //For handling the library
        scalfmm_handle handle = scalfmm_init(chebyshev,source_target);

        //Struct for ref cheb kernel
        struct User_Scalfmm_Cell_Descriptor user_descr;
        user_descr.user_init_cell = NULL;
        user_descr.user_free_cell = NULL;
        //Set algorithm to source target
        scalfmm_algorithm_config(handle,source_target);
        //Build the tree
        scalfmm_build_tree(handle,treeHeight, boxWidth, boxCenter, user_descr);

        //Insert Sources and targets
        scalfmm_tree_insert_particles_xyz(handle,nbPartSource,sourceXYZ,SOURCE);
        printf("Sources inserted \n");
        scalfmm_tree_insert_particles_xyz(handle,nbPartTarget,targetsXYZ,TARGET);
        printf("Targets inserted \n");

        //        scalfmm_print_everything(handle);

        int * arrayofIndicesSource = malloc(sizeof(int)*nbPartSource);
        int * arrayofIndicesTarget = malloc(sizeof(int)*nbPartTarget);
        {//Set physical values
            printf("Setting physical values ... \n");
            Timer phy_timer;
            tic(&phy_timer);
            //SRC
            int idPart;
            for(idPart = 0 ; idPart<nbPartSource ; ++idPart){
                arrayofIndicesSource[idPart] = idPart;
            }
            scalfmm_set_physical_values_npart(handle,nbPartSource,arrayofIndicesSource,sourcePhiValues,SOURCE);
            //TGT
            for(idPart = 0 ; idPart<nbPartTarget ; ++idPart){
                arrayofIndicesTarget[idPart] = idPart; // here, we add the number of sources previously inserted
            }
            scalfmm_set_physical_values_npart(handle,nbPartTarget,arrayofIndicesTarget,targetsPhiValues,TARGET);
            tac(&phy_timer);
            print_elapsed(&phy_timer);
            printf("Physical values setted ... \n");

        }
        printf("Executing FMM ...\n");
        Timer fmm_timer;
        tic(&fmm_timer);
        //Computation
        scalfmm_execute_fmm(handle/*, kernel, &my_data*/);
        tac(&fmm_timer);
        printf("FMM finished ...\n");
        print_elapsed(&fmm_timer);

        //Get back the forces
        scalfmm_get_forces_xyz(handle,nbPartTarget,arrayOfForces,TARGET);
        //Get back the potentials, too
        scalfmm_get_potentials(handle,nbPartTarget,arrayOfPot,TARGET);

        //No need to get Source forces, since it will be 0 anyway
        //scalfmm_get_forces_xyz(handle,nbPartSource,&arrayOfForces[nbPartTarget],SOURCE);

        printf("Forces computed : \n");
        //Display array of forces ::
        //displayPoints(nbPartTarget+nbPartSource,arrayOfForces);

        //Release memory used :
        free(arrayofIndicesSource);
        free(arrayofIndicesTarget);

        scalfmm_dealloc_handle(handle,NULL);

    }

    {//Let's check the result, we computed fr each target part its forces
        //Storage of reference forces + potential
        double * arrayRefForces = malloc(sizeof(double)*nbPartTarget*3);
        memset(arrayRefForces,0,sizeof(double)*nbPartTarget*3);
        double * arrayOfRefPot = malloc(sizeof(double)*nbPartTarget);
        memset(arrayOfRefPot,0,sizeof(double)*nbPartTarget);

        int idTgt;
        Timer check_timer;
        tic(&check_timer);
        for(idTgt = 0 ; idTgt<nbPartTarget ; ++idTgt){
            int idSrc;
            double dx,dy,dz;

            for(idSrc = 0 ; idSrc<nbPartSource ; ++idSrc){
                //First compute dist.
                dx = sourceXYZ[idSrc*3+0] - targetsXYZ[idTgt*3+0];
                dy = sourceXYZ[idSrc*3+1] - targetsXYZ[idTgt*3+1];
                dz = sourceXYZ[idSrc*3+2] - targetsXYZ[idTgt*3+2];

                //Secondly, compute coeff
                double coeffs = targetsPhiValues[idTgt] * sourcePhiValues[idSrc];
                double one_over_r = 1.0/(sqrt(dx*dx+dy*dy+dz*dz));
                double one_over_r3 = one_over_r * one_over_r * one_over_r;

                arrayRefForces[idTgt*3+0] += dx*coeffs*one_over_r3;
                arrayRefForces[idTgt*3+1] += dy*coeffs*one_over_r3;
                arrayRefForces[idTgt*3+2] += dz*coeffs*one_over_r3;
                arrayOfRefPot[idTgt] += one_over_r * sourcePhiValues[idSrc];
            }
        }
        tac(&check_timer);
        print_elapsed(&check_timer);

        {//Then, we compare

            //For L2 norm
            double errorCumulXSquared = 0,
                errorCumulYSquared = 0,
                errorCumulZSquared = 0,
                errorCumulPotSquared = 0;
            //For Inf Norm
            double maxErrorX = 0,
                maxErrorY = 0,
                maxErrorZ = 0,
                maxErrorPot = 0;
            int idArr;
            for(idArr = 0 ; idArr<nbPartTarget ; ++idArr){
                double deltaX = fabs(arrayRefForces[idArr*3+0]-arrayOfForces[idArr*3+0]);
                double deltaY = fabs(arrayRefForces[idArr*3+1]-arrayOfForces[idArr*3+1]);
                double deltaZ = fabs(arrayRefForces[idArr*3+2]-arrayOfForces[idArr*3+2]);
                double deltaPot = fabs(arrayOfRefPot[idArr]-arrayOfPot[idArr]);

                errorCumulXSquared += deltaX*deltaX;
                errorCumulYSquared += deltaY*deltaY;
                errorCumulZSquared += deltaZ*deltaZ;
                errorCumulPotSquared += deltaPot*deltaPot;

                if(maxErrorX < deltaX){ maxErrorX = deltaX ;}
                if(maxErrorY < deltaY){ maxErrorY = deltaY ;}
                if(maxErrorZ < deltaZ){ maxErrorZ = deltaZ ;}
                if(maxErrorPot < deltaPot) {maxErrorPot = deltaPot;}
            }
            //Last check aim to verify the difference between energy
            //total directly computed and energy total FMM's computed
            double energy = 0,energyFmm =0;
            for(idArr = 0; idArr<nbPartTarget ; ++idArr){
                energy += arrayOfRefPot[idArr] * targetsPhiValues[idArr];
                energyFmm += arrayOfPot[idArr] * targetsPhiValues[idArr];
            }



            printf("\n \t\t X error \t Y error \t Z error \t Pot error\n");
            printf("Norme Sup :\t %e \t %e \t %e \t %e\n",maxErrorX,maxErrorY,maxErrorZ,maxErrorPot);
            printf("Norme L2  :\t %e \t %e \t %e \t %e\n",
                   sqrt(errorCumulXSquared),
                   sqrt(errorCumulYSquared),
                   sqrt(errorCumulZSquared),
                   sqrt(errorCumulPotSquared));
            printf("Norme Rms :\t %e \t %e \t %e \t %e\n",
                   sqrt(errorCumulXSquared/((double) nbPartTarget)),
                   sqrt(errorCumulYSquared/((double) nbPartTarget)),
                   sqrt(errorCumulZSquared/((double) nbPartTarget)),
                   sqrt(errorCumulPotSquared/((double) nbPartTarget)));
            printf(" \n Energy Error  : \t %e \n Energy FMM   : \t %e \n Energy DIRECT : \t %e\n",
                   fabs(energy-energyFmm),
                   energyFmm,
                   energy);
        }
        free(arrayRefForces);
        free(arrayOfRefPot);
    }


    //Part where we apply normal on target's forces vector
    //Copying each target's parts forces,
    double * targetsForces = malloc(sizeof(double) * 3 * nbPartTarget);
    memcpy(targetsForces,arrayOfForces,sizeof(double)*3*nbPartTarget);

    double * normeXForces =  malloc(sizeof(double) * nbPartTarget);
    memset(normeXForces,0,sizeof(double) * nbPartTarget);

    computeNormalXForces(nbPartTarget,targetsForces,targetsXYZ,normeXForces);
    printf("For each target, we display [Normal To Sphere] . [Force product] \n");
    //displayArray(nbPartTarget,normeXForces);


    //Free memory
    free(normeXForces);
    free(targetsForces);
    free(sourceXYZ);
    free(sourcePhiValues);
    free(targetsXYZ);
    free(targetsPhiValues);
    free(arrayOfForces);
    free(arrayOfPot);

    return EXIT_SUCCESS;
}
