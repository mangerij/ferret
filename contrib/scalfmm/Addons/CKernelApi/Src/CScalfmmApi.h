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
//
#ifndef CKERNELAPI_H
#define CKERNELAPI_H

typedef long long int FSize;

/**
 * @file
 * This file defines the API for the C USER.  The objective is to
 * provide an unified way to compute the Fast Multipole Method with
 * ScalFMM's algorithm. The mathematical kernel used can be defined by
 * the user (we refer to user_kernel) or one of the ScalFMM kernels
 * (Chebyshev Interpolation or Lagrange Interpolation).
 *
 * @section Scalfmm Handle
 *
 * @section Particle indexes
 * A index is assign to each particle inserted in the tree.
 * It correspond to its position (the first particle has index 0 and so on).
 * It is then possible to interact with particles using these indexes.
 * Moreover, these indexes are used when get/set function are called,
 * (the data for a particle are read/copy using its index).
 *
 * @section Function names
 * Several functions are made to insert, set/get particles data.
 * Most these functions have the following syntaxe rule:
 * @code {void,int} name_{get,set,add}_properties_[xyz]_[nparts](Handle [, parameters]);
 * {void,int} : The function return "int" in case of error code returned.
 * {get,set,add} : In case of "get" the data are copied from the particles,
 * in case of "set" the data are copied to the particles, and
 * in case of "add" the values are add to the particles ones (+=)
 * [xyz] : When the values to work on are in shape of an unique array composed
 * of 3 values per particles x1y1z1.x2y2z2, ...
 * [nparts] : Should be used to work on a part of particles and not all of them.
 * In this case an array of int is given in parameter to give the indexes of the particles
 * to work on.
 */



/////////////////////////////////////////////////////////////////////
//////////////////     Init  Part                      //////////////
/////////////////////////////////////////////////////////////////////

/**
 * Enum over the different kernel usable
 */
typedef enum kernel_type {
    user_defined_kernel = 0,  /** Case if user provides a complete
                                * kernel, ie all the needed function to
                                * compute FMM */
    chebyshev = 1,            /** Case if user uses ScalFMM's implementation of
                                *  Chebyshev Interpolation */
    lagrange = 2,             /** Case if user uses ScalFMM's implementation of
                                *  Lagrange Interpolation*/
} scalfmm_kernel_type;


/* /\** */
/*  * Enum over the different way scalfmm can handle a particule moving */
/*  * out of the simulation box */
/*  *\/ */
/* typedef enum out_of_the_box_config { */
/*     exiting = 0,  /\*If a particule move outside the simulation box, */
/*                     the simulation stop*\/ */
/*     periodic = 1, /\*If a particule move outside the simulation box, */
/*                     the particules is inserted back at the other side */
/*                     of the simulation box*\/ */
/*     erasing = 2   /\*If a particule move outside the simulation box, it */
/*                     simply disappears from the simulation *\/ */
/* } scalfmm_out_of_box_behavior; */

/**
 * Enum over the different algorithm that scalfmm can use
 */
typedef enum scalfmm_algorithm_config {
    sequential = 0,  /* Use the sequential version of Scalfmm*/
    multi_thread = 1, /* Use the Multi thread version of Scalfmm*/
    periodic = 2,    /* Use the periodic version of Scalfmm*/
    source_target = 3, /* USe the source/target algorithm */
    adaptiv = 4, /*Temporary*/
    mpi = 5
} scalfmm_algorithm;


/**
 * Handle for the user
 */
typedef void* scalfmm_handle;



/**
 * @brief This function initialize scalfmm datas.
 * @param KernelType kernel to be used
 * @return scalfmm_handle (ie void *). This handle will be given to
 * every other scalfmm functions.
 *
 * Every data will be stored in order to crete later through a builder
 * what is needed for the simulation
 */
scalfmm_handle scalfmm_init( scalfmm_kernel_type KernelType,scalfmm_algorithm algo);


/////////////////////////////////////////////////////////////////////
//////////////////       Tree  Part                    //////////////
/////////////////////////////////////////////////////////////////////

//In order to initate the tree for user defined cell and kernel, we
//define here the call_backs to deal with cells

/**
 * @brief Function to init the cells (should be given by the user when
 * calling Scalfmm_init_cell)
 * @param level  level of the cell.
 * @param morton_index morton index of the cell to be allocated.
 * @param tree_position int[3] position inside the tree (number of boxes in
 * each direction)
 * @param spatial_position double[3] lower left corner of the cell
 * @param inDatas user generic pointer to kernel.
 */
typedef void* (*Callback_init_cell)(int level, long long morton_index, int* tree_position, double* spatial_position, void * inDatas);

/**
 * @brief Function to destroy what have bee initialized by the user
 * (should be give in Scalfmm_dealloc_handle)
 */
typedef void (*Callback_free_cell)(void*);


/**
 * @brief Callback used to know the size of userData.
 * @param level current level of current cell
 * @param userData Datas that will be serialize
 * @param morton_index of the current cell
 */
typedef FSize (*Callback_get_cell_size)(int level, void * userDatas, long long morton_index);

/**
 * @brief Callback used to serialize userdata inside an array of size
 * given above.
 * @param level current level of current cell
 * @param userData Datas that will be serialize
 * @param morton_index of the current cell
 */
typedef void (*Callback_copy_cell)(void * userDatas, FSize size, void * memoryAllocated);


/**
 * @brief Callback called if scalfmm_finalize_cell is called.
 * @param level current level of leaves (ie height of the tree)
 * @param nbParts Number of particles inside that leaf
 * @param idxParts array of size nbParts, containing the indices of each parts
 * @param morton_index of the current cell
 * @param lower left corner of the current leaf (3 double)
 * @param userData cell user data
 * @param userData Kernel user data
 */
typedef void (*Callback_apply_on_leaf)(int level, FSize nbParts, const FSize * idxParts, long long morton_index, double llc[3],
                                       void * cellDatas,void * leafDatas, void * userDatas);

/**
 * @brief Callback to initialise data inside the Leaves
 * @param level current level of leaves (ie height of the tree)
 * @param nbParts Number of particles inside that leaf"
 * @param idxParts array of size nbParts, containing the indices of each parts
 * @param morton_index of the current cell
 * @param lower left corner of the current leaf (3 double)
 * @param userData leaf user data
 * @param userData Kernel user data
 */
typedef void* (*Callback_init_leaf)(int level, FSize nbParts, const FSize * idxParts, long long morton_index, double llc[3],
                                    void * cellDatas, void * userDatas);

/**
 * @brief Callback to free data inside the Leaves
 * @param nbParts Number of particles inside that leaf
 * @param idxParts array of size nbParts, containing the indices of each parts
 * @param cellDatas ptr to the cell.
 * @param leafData leaf to destroy
 * @param ptr to userData.
 */
typedef void (*Callback_free_leaf)(void * cellDatas, FSize nbParts, const FSize * idxParts, void * leafData, void * userDatas);


/**
 * @brief Callback used to initialize again userDat from what's have
 * been stored inside an array with Callback_copy_cell.
 */
typedef void * (*Callback_restore_cell)(int level, void * arrayTobeRead);


/**
 * @brief Structure containing user's call_backs in order to
 * initialize/free the cell's user data.
 */
typedef struct User_Scalfmm_Cell_Descriptor{
    Callback_free_cell user_free_cell;
    Callback_init_cell user_init_cell;
    Callback_get_cell_size user_get_size;
    Callback_copy_cell user_copy_cell;
    Callback_restore_cell user_restore_cell;
    Callback_init_leaf user_init_leaf;
    Callback_free_leaf user_free_leaf;
}Scalfmm_Cell_Descriptor;


/**
 * @brief This function build the tree. If scalfmm_init has been
 * called with a user defined kernel, then user_cell_descriptor need
 * to be provided
 *
 * @param TreeHeight Height of the octree.
 * @param BoxWidth Width of the entire simulation box.
 * @param BoxCenter Coordinate of the center of the box (ie array)
 */
void scalfmm_build_tree(scalfmm_handle handle,int TreeHeight,double BoxWidth,double* BoxCenter,Scalfmm_Cell_Descriptor user_cell_descriptor);


/**
 * @brief This enum flag is to know if function calling will deal with
 * source, target or both
 */
typedef enum particule_type{
    SOURCE=0,
    TARGET=1,
    BOTH=2
}PartType;


/**
 * @brief This function insert alongside to position an arbitrary
 * number of attributes.
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param NbPartToInsert number of particles to be inserted
 * @param nbAttributeToInsert number of attribute to insert (this
 * number will be > 3 (because we need at least 3 doubles for
 * position))
 * @param strideForEachAtt How to get each attribute for each particles
 * @param rawDatas datas to be read
 *
 * Example :
 struct part{
   double[3] position;
   double charge;
   double test; //not used
   double coeff;
 };
 Then nbAttributeToInsert will be 3+1+1
 and strideForEachAtt will be : [0,1,2,3,5]
 */
void scalfmm_tree_abstract_insert(scalfmm_handle Handle, int NbPartToInsert, int nbAttributeToInsert, int * strideForEachAtt,
                                  double* rawDatas);


/**
 * @brief This function insert an array of position into the
 * octree. THis fonction will insert particules with no SOURCE/TARGET
 * type.
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param NbPositions Number of position to be inserted
 * @param arrayX Array containing the X coordinate for all the parts, size : NbPositions
 * @param arrayY Array containing the Y coordinate for all the parts, size : NbPositions
 * @param arrayZ Array containing the Z coordinate for all the parts, size : NbPositions
 * @param type : type to insert
 * The parts will be inserted with their indices inside the
 * array. Each index will be unique.
 *
 *
 * Default physical values, potential and forces are set to 0.
 */
void scalfmm_tree_insert_particles(scalfmm_handle Handle, int NbPositions, double * arrayX, double * arrayY, double * arrayZ, PartType type);


/**
 * This function is equivalent to scalfmm_tree_insert_particles
 * but the given array XYZ should contains a triple value per paticles.
 */
void scalfmm_tree_insert_particles_xyz(scalfmm_handle Handle, int NbPositions, double * XYZ, PartType type);


/**
 * @brief This function set the physical values of all the particles
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbPhysicalValues Number of physical values to be
 * inserted. Must be equal to the number of positions previously
 * inserted.
 * @param physicalValues Array containing the physical values to be
 * associated to each parts.
 * @param type : type of the particules to be setted.
 * The physical values will be stored according to their indices in
 * the array. First particle inserted will take value physicalValues[0].
 */
void scalfmm_set_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues, PartType type);
/**
 * @brief get the physical values.
 *
 * WARNING : the user must allocate (and initialize) the array given
 */
void scalfmm_get_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues, PartType type);

/**
 *
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbPhysicalValues the number of particles to set the physical values to
 * @param idxOfParticles an array of indexes of size nbPhysicalValues to know which particles
 * to set the values to.
 * @param physicalValues the physical values.
 * @param type : type of the particules to be setted.
 *
 * For example to set the physical values to particles 0 and 1 to values 1.1 and 1.4:
 * @code nbPhysicalValues = 2;
 * @code idxOfParticles = {0 , 1};
 * @code physicalValues = {1.1 , 1.4};
 *
 * Be aware that such approach requiere to find particles in the tree which can have high cost.
 */
void scalfmm_set_physical_values_npart(scalfmm_handle Handle, int nbPhysicalValues,
                                       int* idxOfParticles, double * physicalValues, PartType type);
void scalfmm_get_physical_values_npart(scalfmm_handle Handle, int nbPhysicalValues,
                                       int* idxOfParticles, double * physicalValues, PartType type);


/**
 * @brief This function give back the resulting forces
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbParts total number of particles to retrieve (must be equal
 * to the number of parts inserted)
 * @param forcesToFill array of size nbParts*3, that will contains the
 * forces. WARNING : User must allocate the array before call.
 * @param type : type of the particules to be setted.
 * Forces will be stored sequentially, according to the indices in the
 * array. (ie fx1,fy1,fz1,fx2,fy2,fz2,fx3 ....)
 * @param idxOfParticles : array of indices of the particles
 * wanted. If used, then the parts will be given in the order of
 * idxOfParticles
 */
void scalfmm_get_forces_xyz(scalfmm_handle Handle, int nbParts, double * forcesToFill, PartType type);
void scalfmm_get_forces_xyz_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * forcesToFill, PartType type);
void scalfmm_get_forces(scalfmm_handle Handle, int nbParts, double * fX, double* fY, double* fZ, PartType type);
void scalfmm_get_forces_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ, PartType type);


/**
 * @brief This function give back the resulting forces
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbParts total number of particles to retrieve (must be equal
 * to the number of parts inserted)
 * @param forcesX array of size nbParts, that will contains the
 * forces . WARNING : User must allocate the array before call.
 * @param forcesY array of size nbParts, that will contains the
 * forces . WARNING : User must allocate the array before call.
 * @param forcesZ array of size nbParts, that will contains the
 * forces . WARNING : User must allocate the array before call.
 * @param type : type of the particules to be setted.
 */
void scalfmm_set_forces_xyz(scalfmm_handle Handle, int nbParts, double * forcesToFill, PartType type);
void scalfmm_set_forces_xyz_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * forcesToFill, PartType type);
void scalfmm_set_forces(scalfmm_handle Handle, int nbParts, double * fX, double* fY, double* fZ, PartType type);
void scalfmm_set_forces_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ, PartType type);


/**
 * @brief This function give back the resulting potentials
 * @param Handle scalfmm_handle provided by scalfmm_init
 * @param nbParts total number of particles to retrieve (must be equal
 * to the number of parts inserted)
 * @param potentialsToFill array of potentials to be filled. WARNING :
 * User must allocate the array before call.
 * @param type : type of the particules to be setted.
 *
 * Potentials will be stored sequentially, according to the indices in the
 * array.
 */
void scalfmm_get_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToFill, PartType type);
void scalfmm_set_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToRead, PartType type);
void scalfmm_get_potentials_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * potentialsToFill, PartType type);
void scalfmm_set_potentials_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * potentialsToRead, PartType type);


/**
 * @brief This function update the positions inside the tree, in case
 * of multiple runs of the FMM.
 * @param Handle scalfmm_handle provided by scalfmm_init
 * @param NbPositions total number of particles (must be equal to the
 * number of parts inserted)
 * @param updatedXYZ array of displacement (ie
 * dx1,dy1,dz1,dx2,dy2,dz2,dx3 ...)
 * @param type : type of the particules to be setted.
 */
void scalfmm_add_to_positions_xyz(scalfmm_handle Handle, int NbPositions, double * updatedXYZ, PartType type);
void scalfmm_add_to_positions(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z, PartType type);
void scalfmm_add_to_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ, PartType type);
void scalfmm_add_to_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * X, double * Y , double * Z, PartType type);


/**
 * @brief This function set again the positions inside the tree, in case
 * of multiple runs of the FMM.
 * @param Handle scalfmm_handle provided by scalfmm_init
 * @param NbPositions total number of particles (must be equal to the
 * number of parts inserted)
 * @param newXYZ array of new positions (ie
 * dx1,dy1,dz1,dx2,dy2,dz2,dx3 ...)
 * @param type : type of the particules to be setted.
 * @return Error code, a parts may move out of the simulation
 * box. ScalFMM cannot deals with that specific case. Error code :
 * 0. Success code 1. Could be an arg in order to be Fortran
 * compliant.
 *
 */
void scalfmm_set_positions_xyz(scalfmm_handle Handle, int NbPositions, double * updatedXYZ, PartType type);
void scalfmm_set_positions(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z, PartType type);
void scalfmm_set_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ, PartType type);
void scalfmm_set_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * X, double * Y , double * Z, PartType type);

/**
 *@brief This is function is to be called after a call modifying some
 *of the particles positions (add_to_position or set_position)
 *
 */
void scalfmm_update_tree(scalfmm_handle handle);

void scalfmm_get_positions_xyz(scalfmm_handle Handle, int NbPositions, double * positionsToFill, PartType type);
void scalfmm_get_positions(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z, PartType type);
void scalfmm_get_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * positionsToFill, PartType type);
void scalfmm_get_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * X, double * Y , double * Z, PartType type);

/* /\** */
/*  * @brief This function provides a way for the user to define scalfmm */
/*  * behavior in case a particule get out of the box after a */
/*  * displacement */
/*  * @param  Handle scalfmm_handle provided by scalfmm_init */
/*  * @param  Member of enum scalfmm_out_of_box_behavior */
/*  *\/ */

/* void scalfmm_out_of_the_box_config(scalfmm_handle Handle,scalfmm_out_of_box_behavior config); */

/**
 * @brief This function provides a way for choosing the algorithm to
 * be used
 * @param  Handle scalfmm_handle provided by scalfmm_init
 * @param  Member of enum scalfmm_algorithm
 */

void scalfmm_algorithm_config(scalfmm_handle Handle,scalfmm_algorithm config);


/////////////////////////////////////////////////////////////////////
//////////////////       Kernel  Part                  //////////////
/////////////////////////////////////////////////////////////////////


///////////////// User kernel part :


/**
 * @brief Function to be filled by user's P2M
 * @param nbParticles number of particle in current leaf
 * @param cellData current cell
 * @param leafData currentLeaf
 * @param particleIndexes indexes of particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2M)(void* cellData, void * leafData, FSize nbParticles, const FSize* particleIndexes, void* userData);

/**
 * @brief Function to be filled by user's M2M
 * @param level current level in the tree
 * @prama parentCell cell to be filled
 * @param childPosition number of child (child position in the tree can inferred from its number (refer to doc))
 * @param userData datas specific to the user's kernel
 * @param childCell array of cells to be read
 */
typedef void (*Callback_M2M)(int level, void* parentCell, int childPosition, void* childCell, void* userData);

/**
 * @brief Function to be filled by user's M2L
 * @param level current level in the tree
 * @param targetCell pointer to cell to be filled
 * @param sourceCellPosition number of source cell (source cell
 * position in the tree can inferred from its number (refer to doc))
 * @param sourceCell cell to be read
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_M2L)(int level, void* targetCell, int sourceCellPosition, void* sourceCell, void* userData);

/**
 * @brief Function to be filled by user's M2L
 * @param level current level in the tree
 * @param targetCell pointer to cell to be filled
 * @param sourceCell cell to be read
 * @param transfer array of 3 int, displaying the number of boxes at
 * the current level in along X,Y and Z axis.
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_M2L_Ext)(int level, void* targetCell, void* sourceCell, int tranfer[3], void* userData);

/**
 * @brief Function to be filled by user's M2L
 * @param level current level in the tree
 * @param targetCell pointer to cell to be filled
 * @param sourceCell array of cell to be read
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_M2LFull)(int level, void* targetCell, const int * neighborPosition, const int size, void** sourceCell, void* userData);

/**
 * @brief Function to be filled by user's L2L
 * @param level current level in the tree
 * @param parentCell cell to be read
 * @param childPosition number of child (child position in the tree can inferred from its number (refer to doc))
 * @param childCell cell to be filled
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_L2L)(int level, void* parentCell, int childPosition, void* childCell, void* userData);

/**
 * @brief Function to be filled by user's L2P
 * @param cellData cell to be read
 * @param leafData leaf to be written
 * @param nbParticles number of particles in the current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_L2P)(void* cellData, void * leafData, FSize nbParticles,const FSize* particleIndexes, void* userData);

/**
 * @brief Function to be filled by user's P2P
 * @param targetLeaf ptr to user target leaf
 * @param nbParticles number of particle in current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param sourceLeaf ptr to user source leaf
 * @param nbSourceParticles number of particles in source leaf
 * @param sourceParticleIndexes indexes of cource particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2P)(void * targetLeaf, FSize nbParticles, const FSize* particleIndexes,
                             void * sourceLeaf, FSize nbSourceParticles, const FSize* sourceParticleIndexes, void* userData);

/**
 * @brief Function to be filled by user's P2P
 * @attention This function is symmetric, thus when a call is done
 * between target and neighbors cell, the user needs to apply the target
 * field onto the neighbors cells too.
 * @param targetLeaf ptr to user target leaf
 * @param nbParticles number of particle in current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param sourceLeaves array of ptr to user source leaves.
 * @param sourceParticleIndexes array of indices of source particles currently computed
 * @param sourceNbPart array containing the number of part in each neighbors
 * @param sourcePosition array containing relative position of the neighbor
 * @param size : size of the arrays (thus, number of existing neighbor cell)
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2PFull)(void * targetLeaf, FSize nbParticles, const FSize* particleIndexes,
                                 void ** sourceLeaves,
                                 const FSize ** sourceParticleIndexes,FSize * sourceNbPart,
                                 const int * sourcePosition, const int size, void* userData);


/**
 * @brief Function to be filled by user's P2P inside the leaf
 * @param targetLeaf ptr to user target leaf
 * @param nbParticles number of particle in current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2PInner)(void * targetLeaf,FSize nbParticles, const FSize* particleIndexes, void* userData);

/**
 * @brief Function to be filled by user's P2P
 * @param targetLeaf ptr to user target leaf
 * @param nbParticles number of particle in current leaf
 * @param particleIndexes indexes of particles currently computed
 * @praram sourceLeaf ptr to user source leaf
 * @param nbSourceParticles number of particles in source leaf
 * @param sourceParticleIndexes indexes of cource particles currently computed
 * @param userData datas specific to the user's kernel
 * @attention The fact that this callback is defined (and thus
 * different from NULL or nullptr) means that it will be used, and
 * thus if a call is done with cell1 as source and cell2 as target, no
 * call will be done with cell2 as source and cell1 as target.
 */
typedef void (*Callback_P2PSym)(void * targetLeaf, FSize nbParticles, const FSize* particleIndexes,
                                void * sourceLeaf, FSize nbSourceParticles, const FSize* sourceParticleIndexes, void* userData);

/**
 * @brief Function to be filled by user's method to reset a user's cell
 * @param level  level of the cell.
 * @param morton_index morton index of the cell to be allocated.
 * @param tree_position int[3] position inside the tree (number of boxes in
 * each direction)
 * @param spatial_position double[3] lower left corner of the cell
 * @param usercell ptr to user's cell
 */
typedef void (*Callback_apply_on_cell)(int level, long long morton_index, int* tree_position, double* spatial_position, void * userCell, void * userData);

/**
 * @brief Structure containing callbacks to fill in order to define
 * user kernel.
 *
 */
typedef struct User_Scalfmm_Kernel_Descriptor {
    Callback_P2M p2m;
    Callback_M2M m2m;
    Callback_M2L m2l;
    Callback_M2L_Ext m2l_ext;
    Callback_M2LFull m2l_full;
    Callback_L2L l2l;
    Callback_L2P l2p;
    Callback_P2P p2p;
    Callback_P2PFull p2p_full;
    Callback_P2PInner p2pinner;
    Callback_P2PSym p2p_sym;
}Scalfmm_Kernel_Descriptor;


/**
 * @param user_kernel Scalfmm_Kernel_Descriptor containing callbacks
 * to user Fmm function. Meaningless if using one of our kernel
 * (Chebyshev or Interpolation).

 * @param userDatas Data that will be passed to each FMM
 * function. Can be anything, but allocation/deallocation is user
 * side.
 *
 */
void scalfmm_user_kernel_config(scalfmm_handle Handle, Scalfmm_Kernel_Descriptor userKernel, void * userDatas);




/* void scalfmm_init_cell(scalfmm_handle Handle, Callback_init_cell user_cell_initializer); */

/* void scalfmm_free_cell(scalfmm_handle Handle, Callback_free_cell user_cell_deallocator); */
/* /\** */
/*  *@param Struct defining how scalfmm will handle user's cell data. */
/*  *\/ */
/* void scalfmm_user_cell_config(scalfmm_handle Handle, Scalfmm_Cell_Descriptor user_cell_descriptor); */

///////////////// Common kernel part : /////////////////

/**
 *
 * @brief This function launch the fmm on the parameters given
 * @param Handle scalfmm_handle provided by scalfmm_init
 */
void scalfmm_execute_fmm(scalfmm_handle Handle);

/**
 * @brief This function apply the call back on each leaf. Should be
 * called after insert_parts.
 * @param
 */
void scalfmm_apply_on_leaf(scalfmm_handle Handle, Callback_apply_on_leaf function);


/////////////////////////////////////////////////////////////////////
//////////////////       Dealloc  Part              /////////////////
/////////////////////////////////////////////////////////////////////




/**
 * @brief This function dealloc the scalfmm handle.
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param cellDestroyer Function to be called on the user cell to
 * dealloc
 *
 * The last param cellDestroyer is meaningless in case the user uses
 * one of the provided kernel. (i.e. Chebyshev, Lagrange)
 */
void scalfmm_dealloc_handle(scalfmm_handle handle, Callback_free_cell cellDestroyer);

/**
 * @brief This function apply the function param on each cell of the octree
 * @param Handle scalfmm_handle provided by scalfmm_init.
 */
void scalfmm_apply_on_cell(scalfmm_handle handle, Callback_apply_on_cell function);

/////////////////////////////////////////////////////////////////////
///////////////       Monitoring functions          /////////////////
/////////////////////////////////////////////////////////////////////

/**
 * @brief Scalfmm has a built in feature to get the time elapsed in
 * each operator. This function returns the number of different timers.
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @return Number of timers
 */
int scalfmm_get_nb_timers(scalfmm_handle handle);

/**
 * @brief Scalfmm has a built in feature to get the time elapsed in
 * each operator. This function fill the array with elapsed time for
 * each operator.
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param Array of Timers, to be allocated by the user (using
 * scalfmm_get_nb_timers)
 * Order inside the array : P2M, M2M, M2L, L2L, L2P, P2P, NearField
 * (P2P+L2P).
 */
void scalfmm_get_timers(scalfmm_handle handle,double * Timers);

/////////////////////////////////////////////////////////////////////
///////////////        Algorithms functions         /////////////////
/////////////////////////////////////////////////////////////////////


/**
 * @brief Set the upper limit int the tree for applying FMM : standard
 * value is 2. If used, then user MUST provide a M2L source to target
 * to be applyied on a wider range than the usual 189 cells.
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param upperLimit : int  than 2.
 */
void scalfmm_set_upper_limit(scalfmm_handle handle, int upperLimit);

/////////////////////////////////////////////////////////////////////
///////////////        Distributed Version         //////////////////
/////////////////////////////////////////////////////////////////////

#ifdef SCALFMM_USE_MPI
#warning "IS_THAT_REALLY_WORKING"
//YES

/**
 * @brief Init scalfmm library with MPI
 * @param Same as in Init
 * @param comm Mpi Communicator
 */
scalfmm_handle scalfmm_init_distributed( scalfmm_kernel_type KernelType,scalfmm_algorithm algo, const MPI_Comm comm);

/**
 * @brief Those function are to be called before the insert method
 * @param Handle scalfmm_handle provided by scalfmm_init_distributed.
 * @param nbPoints Number of particles (if local, then it's the number
 * of particles given to that proc, if global, then it's the total
 * number of particles)
 * @param particleXYZ Array of Position. The size is in both case nbPoints.
 * @param localArrayFilled Array that will be filled with particles
 * once the partitionning done. Can be inserted with no changes.
 * @param indexesFilled Array that store the global index of each part
 * in the localArrayFilled.
 * @param stride stride between two attributes inside attr.
 * @param attr array of attribute to be distributed alongside the positions.
 */
void scalfmm_create_local_partition(scalfmm_handle handle, int nbPoints, double * particleXYZ, double ** localArrayFilled,
                                    FSize ** indexesFilled, FSize * outputNbPoint);

void scalfmm_create_global_partition(scalfmm_handle handle, int nbPoints, double * particleXYZ, double ** localArrayFilled,
                                     FSize ** indexesFilled, FSize * outputNbPoint);

/**
 * @brief Once the partition done, one can call this fct in order to
 * partition "things" following the same scheme. Note that arrays must
 * be in the same order as the original parts.
 * @param handle scalfmm_handle provided by scalfmm_init_distributed.
 * @param nbThings number of items
 * @param sizeofthing size of ONE item
 * @param arrayOfThing array of items to be sorted/partitionned
 * @param newArray output array
 */
void scalfmm_generic_partition(scalfmm_handle handle, FSize nbThings, size_t sizeofthing, void * arrayOfThing, void ** newArray);

/**
 * @brief This fct will call delete on its arg, in order to free the
 * memory allocated inside scalfmm, but given back to the user.
 */
void scalfmm_call_delete(void * array);

#endif



#endif
