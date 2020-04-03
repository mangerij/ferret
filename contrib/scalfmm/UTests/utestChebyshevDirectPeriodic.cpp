// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include "Utils/FGlobal.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FRandomLoader.hpp"
#include "Files/FTreeIO.hpp"

#include "Core/FFmmAlgorithmPeriodic.hpp"
#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

#include "FUTester.hpp"
#include "Utils/FMath.hpp"


#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/*
  In this test we compare the Chebyshev fmm results and the direct results.
 */


/** the test class
 *
 */
class TestChebyshevDirect : public FUTester<TestChebyshevDirect> {
    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
              class LeafClass, class OctreeClass, class FmmClass>
    void RunTest(const int PeriodicDeep){
        // Configs
        const int NbLevels      = 4;
        const int SizeSubLevels = 2;

        const FSize NbParticles     = 16;
        FRandomLoader<FReal> loader(NbParticles);
        FReal BoxWidth = loader.getBoxWidth();
        FPoint<FReal> CenterOfBox = loader.getCenterOfBox();

        // Create octrees
        OctreeClass tree(NbLevels, SizeSubLevels,BoxWidth,CenterOfBox);

        struct TestParticle{
            FPoint<FReal> position;
            FReal physicalValue;
            FReal potential;
            FReal forces[3];
        };

        // Insert the particles
        FReal coeff = -1.0, value = 0.10, sum = 0.0;
        TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> position;
            loader.fillParticle(&position);
            // put in tree
            value *= coeff ;
            sum   += value ;
            // put in tree
            tree.insert(position, idxPart, value);
            // get copy
            particles[idxPart].position         = position;
            particles[idxPart].physicalValue = value;
            particles[idxPart].potential        = 0.0;
            particles[idxPart].forces[0]        = 0.0;
            particles[idxPart].forces[1]        = 0.0;
            particles[idxPart].forces[2]        = 0.0;
        }
        if (FMath::Abs(sum)> 0.00001){
            std::cerr << "Sum of charges is not equal zero!!! (sum=<<"<<sum<<" )"<<std::endl;
            exit(-1);
        }
        coeff = FMath::Abs(coeff)* static_cast<int>(NbParticles)  ;
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run FMM computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Fmm...");
        const MatrixKernelClass MatrixKernel;
        FmmClass algo(&tree, PeriodicDeep );
        {

            printf("Levels %d\n", NbLevels);
            printf("PeriodicDeep %d\n", PeriodicDeep);
            printf("Extended tree height %d\n", algo.extendedTreeHeight());
            printf("Width %e\n", BoxWidth);
            printf("Extended Width %e\n", algo.extendedBoxWidth());

            KernelClass *kernels = new KernelClass(algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter(),&MatrixKernel);
            algo.setKernel(kernels);
            algo.execute();
            delete kernels ;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run direct computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Direct...");
        FReal energy= 0.0 , energyD = 0.0 ;
        {
            FTreeCoordinate min, max;
            algo.repetitionsIntervals(&min, &max);
            const int totalRepetitions = (max.getX()-min.getX()+1) * (max.getY()-min.getY()+1) * (max.getZ()-min.getZ()+1);
            printf("Repetitions [%d , %d]x[%d , %d]x[%d , %d] = %d (%lld)\n",min.getX(), max.getX(),
                   min.getY(), max.getY(),min.getZ(), max.getZ(),
                   totalRepetitions, algo.theoricalRepetition()*algo.theoricalRepetition()*algo.theoricalRepetition());
            uassert(totalRepetitions == algo.theoricalRepetition()*algo.theoricalRepetition()*algo.theoricalRepetition());

            #pragma omp parallel for
            for(FSize idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
                for(FSize idxOther =  0 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                    if(idxTarget != idxOther){
                        FP2P::NonMutualParticles(
                                particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
                                particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
                                &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
                                &particles[idxTarget].forces[2],&particles[idxTarget].potential,
                                &MatrixKernel);
                    }

                }
                for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
                    for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
                        for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
                            if(idxX ==0 && idxY == 0 && idxZ == 0) continue;

                            const FPoint<FReal> offset(loader.getBoxWidth() * FReal(idxX),
                                                       loader.getBoxWidth() * FReal(idxY),
                                                       loader.getBoxWidth() * FReal(idxZ));

                            for(int idxSource = 0 ; idxSource < NbParticles ; ++idxSource){
                                TestParticle source;
                                source.position = particles[idxSource].position;
                                source.physicalValue = particles[idxSource].physicalValue;
                                source.potential = particles[idxSource].potential;
                                source.forces[0] = particles[idxSource].forces[0];
                                source.forces[1] = particles[idxSource].forces[1];
                                source.forces[2] = particles[idxSource].forces[2];

                                source.position += offset;
                                FP2P::NonMutualParticles(
                                            source.position.getX(), source.position.getY(),
                                            source.position.getZ(),source.physicalValue,
                                            particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                            particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
                                            &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
                                        &particles[idxTarget].forces[2],&particles[idxTarget].potential,&MatrixKernel);
                            }
                        }
                    }
                }
            }
            for(FSize idx = 0 ; idx < loader.getNumberOfParticles()  ; ++idx){
                energyD +=  particles[idx].potential*particles[idx].physicalValue ;
            }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Compare
        /////////////////////////////////////////////////////////////////////////////////////////////////

        Print("Compute Diff...");
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;

        { // Check that each particle has been summed with all other

            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                    potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                    fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                    fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                    fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
                    energy+=potentials[idxPart]*physicalValues[idxPart];
                }
            });
        }
        // Print for information

        Print("Potential diff is = ");
        printf("         Pot L2Norm     %e\n",potentialDiff.getL2Norm());
        printf("         Pot RL2Norm   %e\n",potentialDiff.getRelativeL2Norm());
        printf("         Pot RMSError   %e\n",potentialDiff.getRMSError());
        Print("Fx diff is = ");
        printf("         Fx L2Norm     %e\n",fx.getL2Norm());
        printf("         Fx RL2Norm   %e\n",fx.getRelativeL2Norm());
        printf("         Fx RMSError   %e\n",fx.getRMSError());
        Print("Fy diff is = ");
        printf("        Fy L2Norm     %e\n",fy.getL2Norm());
        printf("        Fy RL2Norm   %e\n",fy.getRelativeL2Norm());
        printf("        Fy RMSError   %e\n",fy.getRMSError());
        Print("Fz diff is = ");
        printf("        Fz L2Norm     %e\n",fz.getL2Norm());
        printf("        Fz RL2Norm   %e\n",fz.getRelativeL2Norm());
        printf("        Fz RMSError   %e\n",fz.getRMSError());
        FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
        printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
        printf("  Energy Error  =   %.12e\n",FMath::Abs(energy-energyD));
        printf("  R Energy Error  =   %.12e\n",FMath::Abs(energy-energyD)/FMath::Abs(energyD));
        printf("  Energy FMM    =   %.12e\n",FMath::Abs(energy));
        printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyD));

        // Assert
        const FReal MaximumDiffPotential = FReal(9e-3);
        const FReal MaximumDiffForces     = FReal(9e-2);

        uassert(potentialDiff.getL2Norm() < MaximumDiffPotential);    //1
        //     uassert(potentialDiff.getRMSError() < MaximumDiffPotential);  //2
        uassert(fx.getL2Norm()  < MaximumDiffForces);                       //3
        //      uassert(fx.getRMSError() < MaximumDiffForces);                      //4
        uassert(fy.getL2Norm()  < MaximumDiffForces);                       //5
        //       uassert(fy.getRMSError() < MaximumDiffForces);                      //6
        uassert(fz.getL2Norm()  < MaximumDiffForces);                      //8
        //      uassert(fz.getRMSError() < MaximumDiffForces);                                           //8
        uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
        uassert(FMath::Abs(energy-energyD) / FMath::Abs(energyD) < coeff*MaximumDiffPotential);                     //10  Total Energy

        delete[] particles;
    }


    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
              class LeafClass, class OctreeClass, class FmmClass,
              class FmmClassNonPer>
    void RunTestUpward(const int PeriodicDeep)	{
        // Configs
        const int NbLevels      = 4;
        const int SizeSubLevels = 2;

        const FSize NbParticles     = 150;
        FRandomLoader<FReal> loader(NbParticles);
        FReal BoxWidth = loader.getBoxWidth();
        FPoint<FReal> CenterOfBox = loader.getCenterOfBox();

        // Create octrees
        OctreeClass tree(NbLevels, SizeSubLevels,BoxWidth,CenterOfBox);
        OctreeClass treeNonper(NbLevels, SizeSubLevels,BoxWidth,CenterOfBox);

        struct TestParticle{
            FPoint<FReal> position;
            FReal physicalValue;
            FReal potential;
            FReal forces[3];
        };

        // Insert the particles
        FReal coeff = -1.0, value = 0.10, sum = 0.0;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> position;
            loader.fillParticle(&position);
            // put in tree
            value *= coeff ;
            sum   += value ;
            // put in tree
            tree.insert(position, idxPart, value);
            treeNonper.insert(position, idxPart, value);
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run FMM computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Fmm...");
        const MatrixKernelClass MatrixKernel;
        {
            FmmClass algo(&tree, PeriodicDeep );

            printf("Levels %d\n", NbLevels);
            printf("PeriodicDeep %d\n", PeriodicDeep);
            printf("Extended tree height %d\n", algo.extendedTreeHeight());
            printf("Width %e\n", BoxWidth);
            printf("Extended Width %e\n", algo.extendedBoxWidth());

            FTreeCoordinate min, max;
            algo.repetitionsIntervals(&min, &max);
            const int totalRepetitions = (max.getX()-min.getX()+1) * (max.getY()-min.getY()+1) * (max.getZ()-min.getZ()+1);
            printf("Repetitions [%d , %d]x[%d , %d]x[%d , %d] = %d (%lld)\n",min.getX(), max.getX(),
                   min.getY(), max.getY(),min.getZ(), max.getZ(),
                   totalRepetitions, algo.theoricalRepetition()*algo.theoricalRepetition()*algo.theoricalRepetition());

            KernelClass *kernels = new KernelClass(algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter(),&MatrixKernel);
            algo.setKernel(kernels);
            algo.execute();
            delete kernels ;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run upper non periodic FMM computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Upper non periodic Fmm...");
        {
            KernelClass *kernels = new KernelClass(NbLevels, BoxWidth, CenterOfBox,&MatrixKernel);
            FmmClassNonPer algoNonPer(&treeNonper, kernels);
            algoNonPer.execute(FFmmP2M | FFmmM2M);
            delete kernels ;

            typename OctreeClass::Iterator iter(&tree);
            typename OctreeClass::Iterator iterNonPer(&treeNonper);

            iter.gotoBottomLeft();
            iterNonPer.gotoBottomLeft();


            for(int idxLevel = NbLevels-1 ; idxLevel >= 2 ; --idxLevel){
                while(true){
                    uassert(iter.getCurrentGlobalIndex() == iterNonPer.getCurrentGlobalIndex());

                    CellClass* cl = iter.getCurrentCell();
                    CellClass* clnonper = iterNonPer.getCurrentCell();

                    FMath::FAccurater<FReal> accurater;

                    for(int idx = 0; idx < cl->getVectorSize(); ++idx){
                        accurater.add(cl->getMultipole(0)[idx],
                                      clnonper->getMultipole(0)[idx]);
                    }
                    uassert(accurater.getInfNorm() < 1E-10);
                    uassert(accurater.getL2Norm() < 1E-10);

                    const bool hasNext = iter.moveRight();
                    const bool hasNextNonPer = iterNonPer.moveRight();
                    uassert(hasNext == hasNextNonPer);
                    if(!hasNext || !hasNextNonPer){
                        break;
                    }
                }
                iter.moveUp();
                iter.gotoLeft();
                iterNonPer.moveUp();
                iterNonPer.gotoLeft();
            }
        }
    }


    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
              class LeafClass, class OctreeClass, class FmmClass,
              class FmmClassNonPer>
    void RunTestFake(const int PeriodicDeep)	{
        // Configs
        const int NbLevels      = 4;
        const int SizeSubLevels = 2;

        const FSize NbParticles     = 150;
        FRandomLoader<FReal> loader(NbParticles);
        FReal BoxWidth = loader.getBoxWidth();
        FPoint<FReal> CenterOfBox = loader.getCenterOfBox();

        // Create octrees
        OctreeClass tree(NbLevels, SizeSubLevels,BoxWidth,CenterOfBox);

        struct TestParticle{
            FPoint<FReal> position;
            FReal physicalValue;
            FReal potential;
            FReal forces[3];
        };

        // Insert the particles
        FReal coeff = -1.0, value = 0.10, sum = 0.0;
        TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> position;
            loader.fillParticle(&position);
            // put in tree
            value *= coeff ;
            sum   += value ;
            // put in tree
            tree.insert(position, idxPart, value);
            // get copy
            particles[idxPart].position         = position;
            particles[idxPart].physicalValue    = value;
            particles[idxPart].potential        = 0.0;
            particles[idxPart].forces[0]        = 0.0;
            particles[idxPart].forces[1]        = 0.0;
            particles[idxPart].forces[2]        = 0.0;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run FMM computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Fmm...");
        const MatrixKernelClass MatrixKernel;
        FmmClass algo(&tree, PeriodicDeep );
        {

            printf("Levels %d\n", NbLevels);
            printf("PeriodicDeep %d\n", PeriodicDeep);
            printf("Extended tree height %d\n", algo.extendedTreeHeight());
            printf("Width %e\n", BoxWidth);
            printf("Extended Width %e\n", algo.extendedBoxWidth());

            KernelClass *kernels = new KernelClass(algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter(),&MatrixKernel);
            algo.setKernel(kernels);
            algo.execute();
            delete kernels ;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Run direct computation
        /////////////////////////////////////////////////////////////////////////////////////////////////
        Print("Direct...");
        printf("Extended tree height boundary %d\n", algo.extendedTreeHeightBoundary());
        printf("Extended Width boundary %e\n", algo.extendedBoxWidthBoundary());
        OctreeClass treeNonper(algo.extendedTreeHeightBoundary(), SizeSubLevels, algo.extendedBoxWidthBoundary(), algo.extendedBoxCenterBoundary());
        {
            FTreeCoordinate min, max;
            algo.repetitionsIntervals(&min, &max);
            const int totalRepetitions = (max.getX()-min.getX()+1) * (max.getY()-min.getY()+1) * (max.getZ()-min.getZ()+1);
            printf("Repetitions [%d , %d]x[%d , %d]x[%d , %d] = %d (%lld)\n",min.getX(), max.getX(),
                   min.getY(), max.getY(),min.getZ(), max.getZ(),
                   totalRepetitions, algo.theoricalRepetition()*algo.theoricalRepetition()*algo.theoricalRepetition());
            uassert(totalRepetitions == algo.theoricalRepetition()*algo.theoricalRepetition()*algo.theoricalRepetition());

            for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
                for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
                    for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
                        const FPoint<FReal> offset(loader.getBoxWidth() * FReal(idxX),
                                                   loader.getBoxWidth() * FReal(idxY),
                                                   loader.getBoxWidth() * FReal(idxZ));

                        for(int idxSource = 0 ; idxSource < NbParticles ; ++idxSource){
                            TestParticle source = particles[idxSource];
                            source.position += offset;

                            treeNonper.insert(source.position, idxSource, source.physicalValue);
                        }
                    }
                }
            }
        }
        {
            KernelClass *kernels = new KernelClass(algo.extendedTreeHeightBoundary(), algo.extendedBoxWidthBoundary(), algo.extendedBoxCenterBoundary(),&MatrixKernel);
            FmmClassNonPer algoNonPer(&treeNonper, kernels);
            algoNonPer.execute();
            delete kernels ;
        }
        {
            typename OctreeClass::Iterator iter(&tree);

            iter.gotoBottomLeft();

            const int diffLevel = algo.extendedTreeHeightBoundary()-NbLevels;
            const int idxLevel = NbLevels-1;
            const int offsetAtLevel = (PeriodicDeep == -1? 1<<idxLevel : (3 << idxLevel));
            do {
                const FTreeCoordinate realCoord = iter.getCurrentGlobalCoordinate();
                const FTreeCoordinate boundaryCoord(realCoord.getX() + offsetAtLevel,
                                                    realCoord.getY() + offsetAtLevel,
                                                    realCoord.getZ() + offsetAtLevel);
                const MortonIndex boundaryIndex = boundaryCoord.getMortonIndex(idxLevel + diffLevel);

                uassert((boundaryIndex & ~((~0LL)<<3*idxLevel))== iter.getCurrentGlobalIndex());

                ContainerClass* lfnonper = treeNonper.getLeafSrc(boundaryIndex);
                uassert(lfnonper);
                if(lfnonper){
                    ContainerClass* lf = iter.getCurrentLeaf()->getSrc();

                    uassert(lfnonper->getNbParticles() == lf->getNbParticles());

                    if(lfnonper->getNbParticles() == lf->getNbParticles()){
                        FMath::FAccurater<FReal> potentialDiff;
                        FMath::FAccurater<FReal> fx, fy, fz;

                        const FReal*const potentials        = lf->getPotentials();
                        const FReal*const forcesX = lf->getForcesX();
                        const FReal*const forcesY = lf->getForcesY();
                        const FReal*const forcesZ = lf->getForcesZ();
                        const FSize nbParticlesInLeaf = lf->getNbParticles();
                        const FVector<FSize>& indexes = lf->getIndexes();


                        const FReal*const potentialsnonper        = lfnonper->getPotentials();
                        const FReal*const forcesXnonper = lfnonper->getForcesX();
                        const FReal*const forcesYnonper = lfnonper->getForcesY();
                        const FReal*const forcesZnonper = lfnonper->getForcesZ();
                        const FVector<FSize>& indexesnonper = lfnonper->getIndexes();

                        for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                            uassert(indexes[idxPart] == indexesnonper[idxPart]);
                            potentialDiff.add(potentials[idxPart],potentialsnonper[idxPart]);
                            fx.add(forcesX[idxPart],forcesXnonper[idxPart]);
                            fy.add(forcesY[idxPart],forcesYnonper[idxPart]);
                            fz.add(forcesZ[idxPart],forcesZnonper[idxPart]);
                        }
                        uassert(potentialDiff.getInfNorm() < 1E-10);
                        uassert(potentialDiff.getL2Norm() < 1E-10);
                        uassert(fx.getInfNorm() < 1E-10);
                        uassert(fx.getL2Norm() < 1E-10);
                        uassert(fy.getInfNorm() < 1E-10);
                        uassert(fy.getL2Norm() < 1E-10);
                        uassert(fz.getInfNorm() < 1E-10);
                        uassert(fz.getL2Norm() < 1E-10);
                    }
                }
            } while(iter.moveRight());
        }
        {
            typename OctreeClass::Iterator iter(&tree);

            iter.gotoBottomLeft();

            const int diffLevel = algo.extendedTreeHeightBoundary()-NbLevels;

            for(int idxLevel = NbLevels-1 ; idxLevel >= 1 ; --idxLevel){
                const int offsetAtLevel = (PeriodicDeep == -1? 1<<idxLevel : (3 << idxLevel));
                do {
                    const FTreeCoordinate realCoord = iter.getCurrentGlobalCoordinate();
                    const FTreeCoordinate boundaryCoord(realCoord.getX() + offsetAtLevel,
                                                        realCoord.getY() + offsetAtLevel,
                                                        realCoord.getZ() + offsetAtLevel);
                    const MortonIndex boundaryIndex = boundaryCoord.getMortonIndex(idxLevel + diffLevel);

                    uassert((boundaryIndex & ~((~0LL)<<3*idxLevel))== iter.getCurrentGlobalIndex());

                    CellClass* clnonper = treeNonper.getCell(boundaryIndex, idxLevel + diffLevel);
                    uassert(clnonper);
                    if(clnonper){
                        CellClass* cl = iter.getCurrentCell();

                        FMath::FAccurater<FReal> accurater;

                        for(int idx = 0; idx < cl->getVectorSize(); ++idx){
                            accurater.add(cl->getMultipole(0)[idx],
                                          clnonper->getMultipole(0)[idx]);
                        }
                        uassert(accurater.getInfNorm() < 1E-10);
                        uassert(accurater.getL2Norm() < 1E-10);

                        FMath::FAccurater<FReal> accuraterLocal;

                        for(int idx = 0; idx < cl->getVectorSize(); ++idx){
                            accuraterLocal.add(cl->getLocal(0)[idx],
                                               clnonper->getLocal(0)[idx]);
                        }
                        uassert(accuraterLocal.getInfNorm() < 1E-8);
                        uassert(accuraterLocal.getL2Norm() < 1E-8);
                    }

                } while(iter.moveRight());
                iter.moveUp();
                iter.gotoLeft();
            }
        }
        delete[] particles;
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                      << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated()
                      << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated()
                      << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }


    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////


    /** TestChebSymKernel */
    void TestChebSymKernel(){
        const unsigned int ORDER = 7;
        typedef double FReal;
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FFmmAlgorithmPeriodic<FReal,OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        //typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClassNonPer;
        typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClassNonPer;
        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(-1);
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(0);
        //RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(1);
        RunTestUpward<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass,FmmClassNonPer>(-1);
        RunTestUpward<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass,FmmClassNonPer>(0);
        RunTestUpward<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass,FmmClassNonPer>(1);
        RunTestFake<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass,FmmClassNonPer>(-1);
        RunTestFake<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass,FmmClassNonPer>(0);
    }



    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestChebyshevDirect::TestChebSymKernel,"Test Chebyshev Kernel with 16 small SVDs and symmetries");
    }
};


// You must do this
TestClass(TestChebyshevDirect)




