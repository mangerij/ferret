// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Cyrille Piacibello
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

#include <string>

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FLog.hpp"
#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"
#include "../../Src/Utils/FParameterNames.hpp"

int main(int argc,char** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Compile a Taylor Kernel (but do nothing).");
    const int P = 10;
    const int order = 1;
    typedef double FReal;
    FPoint<FReal> centerBox = FPoint<FReal>(0,0,0);

    typedef FTaylorCell<FReal, P,order> CellClass;
    typedef FP2PParticleContainer<FReal> ContainerClass;
    typedef FTaylorKernel<FReal,CellClass,ContainerClass,P,order> KernelClass;
    //typedef FSimpleLeaf<FReal, ContainerClass > LeafClass;
    //typedef FOctree<FReal, CellClass, ContainerClass , LeafClass > OctreeClass;

    KernelClass kernel(9,1.0,centerBox);

    return 0;
}
