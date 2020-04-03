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
#ifndef FMMAPI_H
#define FMMAPI_H

enum FmmApiErrors {
    FMMAPI_NO_ERROR,
    FMMAPI_SUPPORTED_PARAMETER,
    FMMAPI_UNSUPPORTED_PARAMETER,
    FMMAPI_UNKNOWN_PARAMETER
};

////////////////////// Opérateurs FMM Core : //////////////////////////

enum FmmApiCoreParameters {
    FMMCORE_TREE_HEIGHT, // hombre de niveaux de l'arbre (int)
    FMMCORE_ROOT_BOX_WIDTH, // taille de la boîte racine (FReal)
    FMMCORE_ROOT_BOX_CENTER, // position du centre de la boîte racine (FReal[3])
    FMMCORE_LEAF_BOX_WIDTH, // taille des boîtes feuilles (FReal)
    FMMCORE_POINTS_PER_LEAF, // nombre moyen de points par feuille (FReal)
    FMMCORE_MPI_COMMUNICATOR, // communicateur MPI (MPI_Comm)
    FMMCORE_THREADS_NUMBER, // nombre de threads (int)
    FMMCORE_THREAD_ID, // id du thread (int)
    FMMCORE_RHS_NUMBER, // nombre de seconds membres (int)
    //paramètres en lecture seule :
    FMMCORE_HANDLES_P2P // renvoie 0 ou 1 pour dire si le FmmCore gère ou pas le P2P.
};

int FmmCore_init(void **fmmCore) ; /*alloue et initialise le FmmCore*/
int FmmCore_free(void *fmmCore) ; /*libère le FmmCore*/
int FmmCore_isParameterUsed(void */*fmmCore*/, int *name, int *flag);
int FmmCore_setParameter(void *fmmCore, int *name, void*value);
int FmmCore_setParameter(void *fmmCore, int name, void*value);
int FmmCore_getParameter(void *fmmCore, int *name, void*value);
int FmmCore_getParameter(void *fmmCore, int name, void*value);


int FmmCore_getRadius(void*fmmCore, void *boxId, FReal *radius);/*Renvoie le rayon de la boîte*/
int FmmCore_getCentre(void*/*fmmCore*/, void *boxId, FReal **centre); /*Renvoie dans le FReal[3] centre les coordonnées du centre*/
int FmmCore_getLevel(void*/*fmmCore*/, void *boxId, int *level); /*Renvoie dans level le niveau de la boîte boxId dans son arbre*/
int FmmCore_getMultipoleArray(void* /*fmmCore*/, void *boxId, void **F); /*Renvoie dans F l'adresse où stocker l'expansion multipôle associée à la boîte boxId*/
int FmmCore_getLocalArray(void* /*fmmCore*/, void *boxId, void **F); /*Renvoie dans F l'adresse où stocker l'expansion locale associée à la boîte boxId*/
int FmmCore_getCoord(void*/*fmmCore*/, void *boxId, int *coord); /*Renvoie dans coord la position dans l'arbre*/


/* Données potentiel/champ */
int FmmCore_getSource(void* /*fmmCore*/, void *boxId, FReal** position, void** potential, int *number); /* Appelé par P2P et P2M pour obtenir le nombre, la position et le potentiel des sources.
Les différents tableaux sont (éventuellement) alloués par le FmmCore. */
int FmmCore_releaseSource(void*fmmCore, void *boxId, void* potential, FReal* position); /* si le core veut libérer ces tableaux potentiel et position.*/
int FmmCore_getTargetPoints(void* /*fmmCore*/, void *boxId, FReal** position, int *number) ; /* : Appelé par P2P et L2P pour obtenir le nombre et la position des points cibles.*/
int FmmCore_releaseTargetPoints(void*fmmCore, void *boxId, FReal* position); /* si le core veut libérer ce tableau "position".*/
int FmmCore_getTargetField(void* /*fmmCore*/, void *boxId, FReal* F); /* obtient dans un tableau F alloué/libéré par L2P/P2P les valeurs des champs aux points cible
(pour le cas où P2P et L2P doivent sommer leurs résultats).*/
int FmmCore_setTargetField(void* /*fmmCore*/, void *boxId, FReal* F); /* transmets au FmmCore dans F les valeurs des champs aux points cibles mis à jour.*/

/* Entrée/sortie principale */
int FmmCore_setKernelData(void *fmmCore, void *fmmKernel); /* stocke l'identifiant du FmmKernel dans le FmmCore. Cela permet par la suite aux
opérateurs FMM d'accéder au FmmKernel, et donc d'accéder aux données spécifiques au kernel
(p.ex. fréquence dans le cas Helmholtz, …)*/
int FmmCore_getKernelData(void*fmmCore, void **fmmKernel); /* récupère l'identifiant du FmmKernel. */
int FmmCore_setPositions(void *fmmCore, int *nb, FReal *position) ; /* transmet au FmmCore les potentiels associés aux points sources.
Le tableau potential est alloué et libéré par la routine appelant, le FmmCore doit donc en faire une copie.*/
int FmmCore_setPotentials(void *fmmCore, void *potentials); /* transmet au FmmCore les potentiels associés aux points sources.
Le tableau potential est alloué et libéré par la routine appelant, le FmmCore doit donc en faire une copie.*/
int FmmCore_doComputation(void *fmmCore) ; /* réalise le produit multipôle. */
/* !!! Warning use *filed and not **field */
int FmmCore_getField(void *fmmCore, void *fields) ;/* récupère après le produit multipôle la valeur des champs en chaque point.
Le tableau field doit être alloué et libéré par la routine appelante. */

////////////////////// Opérateurs FMM Kernel : //////////////////////////

enum FmmApiKernelParameters {
    FMMKERNEL_ACCURACY,             // précision demandée à la FMM : 1e-3, 1e-6, ... (FReal)
    FMMKERNEL_POTENTIAL_DATA_SIZE,  // taille en octet de la donnée "potentiel" pour 1 point source (int)
    FMMKERNEL_FIELD_DATA_SIZE,      // taille en octet de la donnée "field" pour 1 point cible (int)
    //paramètres en lecture seule :
    FMMKERNEL_HANDLES_P2P           // renvoie 0 ou 1 pour dire si le FmmKernel gère ou pas le P2P.
};

/******* Allocation : ******/
int FmmKernel_init(void *fmmCore, void **fmmKernel);/* : alloue et initialise le FmmKernel */
int FmmKernel_free(void *fmmKernel); /* libére le FmmKernel */

/******* Configuration : ***/

int FmmKernel_isParameterUsed(void * /*fmm*/, int *name, int *flag);
int FmmKernel_setParameter(void *fmmKernel, int *name, void*value);
int FmmKernel_setParameter(void *fmmKernel, int name, void*value);
int FmmKernel_getParameter(void *fmmKernel, int *name, void*value);
int FmmKernel_getParameter(void *fmmKernel, int name, void*value);

/****** Données FMM : *****/
int FmmKernel_getMultipoleArraySize(void *fmmCore, int *size); /* Renvoie dans size la taille (en octets) de l'expansion multipôle associée à la boîte boxId */
int FmmKernel_getLocalArraySize(void *fmmCore, int *size); /* Renvoie dans size la taille (en octets) de l'expansion locale associée à la boîte boxId*/

/******* Opérateurs FMM : **/
int FmmKernel_P2M(void *fmmCore, void* boxId);
int FmmKernel_L2P(void *fmmCore, void* boxId);
int FmmKernel_M2M(void *fmmCore, void *boxIdFather, void *boxIdSon);
int FmmKernel_L2L(void *fmmCore, void *boxIdFather, void *boxIdSon);
int FmmKernel_M2L(void *fmmCore, void *boxIdSrc, void *boxIdDest);
int FmmKernel_P2P_inner(void *fmmCore, void *boxIdSrcDest);
int FmmKernel_P2P(void *fmmCore, void *boxIdSrc, void *boxIdDest); /* pas mutuel, i.e. on fait seulement dans 1 sens. */


#endif // FMMAPI_H
