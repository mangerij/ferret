/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//
#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FPoint.hpp"
#include "Files/FGenerateDistribution.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Files/FExportWriter.hpp"

#include "Utils/FParameterNames.hpp"

//
/// \file  generateDistributions.cpp
//!
//! \brief generateDistributions: Driver to generate N points (non)uniformly distributed on a given geometry
//!
//! The goal of this driver is to generate uniform or non uniform points on the following geometries
//!
//!   Uniform : cube, cuboid, sphere, prolate,
//!
//!   Non uniform : ellipsoid, prolate
//!
//!  You can set two kind of physical values depending of your problem. By default all values are between 0 and 1.
//!   If you select the argument -charge (see bellow) the values are between -1 and 1.
//!  The arguments available are
//!
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param  -N     The number of points in the distribution (default 20000)
//!     \param   -fout name: generic name for files (with extension) and save data
//!                  with following format in name.fma or name.bfma in -bin is set"
//!      \param  -fvisuout Filename for the visu file (vtk, vtp, cvs or cosmo). vtp is the default
//!      \param -extraLength   value    extra length to add to the boxWidth (default 0.0)

//!  <b> Geometry arguments:</b>
//!      \param  -unitCube uniform distribution on unit cube
//!      \param  -cube uniform distribution on a cube
//!          \arg         -length  R - default value for R is 2.0
//!      \param  -unitSphere uniform distribution on unit sphere
//!      \param  -sphere  uniform distribution on  sphere of radius given by
//!          \arg         -radius  R - default value for R is 2.0
//!        \param   -ellipsoid non uniform distribution on  an ellipsoid of aspect ratio given by
//!              \arg          -ar a:b:c   with a, b and c > 0
//!         \param  -prolate ellipsoid with aspect ratio a:a:c  given by
//!                \arg             -ar a:a:c   with  c > a > 0
//!          \param   -plummer (Highly non uniform) plummer distribution (astrophysics)
//!                   \arg         -radius  R - default value 10.0"
//!
//!
//!  <b> Physical values argument:</b>
//!         \param -charge generate physical values between -1 and 1 otherwise generate between 0 and 1
//!         \param -zeromean  the average of the physical values is zero
//!
//!
//! \b examples
//!
//!   generateDistributions -prolate -ar 2:2:4   -N 20000 -fout prolate
//!
//! or
//!
//!  generateDistributions -cuboid 2:2:4 -N 100000 -fout cuboid.bfma  -fvisuout cuboid.vtp -charge  -zeromean
//!



int main(int argc, char ** argv){
    const FParameterNames LocalOptionEllipsoid = {
        {"-ellipsoid"} ,
        " non uniform distribution on  an ellipsoid of aspect ratio given by \n  -ar a:b:c   with a, b and c > 0\n"
    };
    FHelpDescribeAndExit(argc, argv,
                         ">> Driver to generate N points (non)uniformly distributed on a given geometry.\n"
                         "Options  \n"
                         "   -help       to see the parameters    \n"
                         "   -N       The number of points in the distribution    \n"
                         "   -extraLength   value    extra length to add to the boxWidth\n"
                         "    Distributions   \n"
                         "        Uniform on    \n"
                         "             -unitCube  uniform distribution on unit cube\n"
                         "             -cuboid  uniform distribution on rectangular cuboid of size  a:b:c\n"
                         "                     -lengths   a:b:c - default values are 1.0:1.0:2.0\n"
                         "             -unitSphere  uniform distribution on unit sphere\n"
                         "             -sphere  uniform distribution on  sphere of radius given by\n"
                         "                     -radius  R - default value for R is 2.0\n"
                         "             -prolate ellipsoid with aspect ratio a:a:c\n"
                         "                     -ar a:a:c   with  c > a > 0\n"
                         "        Non Uniform on    \n"
                         "             -ellipsoid non uniform distribution on  an ellipsoid of aspect ratio given by\n"
                         "                     -ar a:b:c   with a, b and c > 0\n"
                         "             -plummer (Highly non unuiform) plummer distrinution (astrophysics)\n"
                         "                     -radius  R - default value 10.0\n"
                         "    Physical values\n"
                         "             -charge generate physical values between -1 and 1 otherwise generate between 0 and 1		\n"
                         "             -zeromean  the average of the physical values is zero		\n",
//                         "     Output \n"
 //                        "             -filename name: generic name for files (without extension) and save data\n"
  //                       "                     with following format in name.fma or name.bfma in -bin is set\n"
   //                      "             -visufmt  vtk, vtp, cosmo or cvs format.",
                          FParameterDefinitions::OutputFile,
                         FParameterDefinitions::NbParticles,FParameterDefinitions::OutputVisuFile,LocalOptionEllipsoid);


    
    typedef double FReal;
    FReal       extraRadius = 0.000 ;

    const FSize NbPoints  = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options,   FSize(20000));
    const std::string genericFileName(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "unifPointDist"));

    FReal BoxWith = 0.0;
    FPoint<FReal> Centre(0.0, 0.0,0.0);
	//
	// Allocation
	//
	FReal * particles ;
	particles = new FReal[4*NbPoints] ;
	memset(particles,0,4*NbPoints*sizeof(FReal));
    FmaRWParticle<FReal, 4,4> *ppart = (FmaRWParticle<FReal, 4,4>*)(&particles[0]);

	//
	// Generate physical values
	//

	FReal phyVal, sum,a,b ;
	if(FParameters::existParameter(argc, argv, "-charge")){
		a= 2.0 ; b = -1.0 ;
	}
	else {
		a= 1.0 ; b = 0.0 ;
	}
	sum = 0.0 ;
	int j = 3 ;
	for(int i = 0 ; i< NbPoints; ++i, j+=4){
        phyVal            = a*getRandom<FReal>() +b  ;
		sum              += phyVal ;
		particles[j]       = phyVal ;
	}
	if(FParameters::existParameter(argc, argv, "-zeromean")){
        FReal  rm = FReal(sum)/FReal(NbPoints) ; sum = 0.0 ;
		j = 3 ;
		for(int i = 0 ; i< NbPoints; ++i, j+=4){
			particles[j]    -= rm ;
			sum              += particles[j]   ;
		}
	}
    std::cout << "Sum physical value "<< sum << "   Mean Value " << sum/FReal(NbPoints)<<std::endl ;
	//
	// Point  generation
	//
	if(FParameters::existParameter(argc, argv, "-unitCube")){
		unifRandonPointsOnUnitCube(NbPoints, particles) ;
		Centre.setPosition(0.5,0.5,0.5);
		BoxWith = 1.0 ;
	}
	else if(FParameters::existParameter(argc, argv, "-cuboid")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-lengths",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		         aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
		unifRandonPointsOnCube(NbPoints, A,B,C,particles) ;
		BoxWith = FMath::Max(A,FMath::Max(B,C) );
		FReal halfBW = BoxWith*0.5;
		Centre.setPosition(halfBW,halfBW,halfBW);
		std::cout << "Cuboid "<< A << ":"<< B<<":"<<C<<std::endl;
	}
	else if(FParameters::existParameter(argc, argv, "-unitSphere")){
		unifRandonPointsOnUnitSphere(NbPoints, particles) ;
		BoxWith = 2.0 ;
	}
	else if(FParameters::existParameter(argc, argv, "-sphere")){
		const FReal Radius  = FParameters::getValue(argc,argv,"-radius",  2.0);
		unifRandonPointsOnSphere(NbPoints, Radius,particles) ;
		BoxWith = 2.0*Radius ;
	}
	else if(FParameters::existParameter(argc, argv, "-prolate")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
		if(A != B){
			std::cerr << " A /= B in prolate ellipsoide A =B. Your aspect ratio: "<< aspectRatio<<std::endl;
		}
		std::cout << "A: "<<A<<" B: "<< B << " C: " << C<<std::endl;
		unifRandonPointsOnProlate(NbPoints,A,C,particles);
		BoxWith =  2.0*C;
	}    //const FSize NbPoints  = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options,   FSize(20000));
    else if(FParameters::existParameter(argc, argv, "-hyperpara")){
        std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
        FReal A,B,C ;
        size_t pos = aspectRatio.find(":");     aspectRatio.replace(pos,1," ");
        pos = aspectRatio.find(":");        aspectRatio.replace(pos,1," ");
        std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
        std::cout << "A: "<<A<<" B: "<< B << " C: " << C<<std::endl;
        unifRandonPointsOnHyperPara(NbPoints,A,B,C,particles);
        BoxWith =  2.0*FMath::Max( A,FMath::Max( B,C)) ;
        std::cout << "BoxWith: " << BoxWith<<std::endl;

    }
	else if(FParameters::existParameter(argc, argv, "-ellipsoid")){
//		else if(FParameters::existParameter(argc, argv, "-ellipsoid")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
//		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
			std::cout << "A: "<<A<<" B: "<< B << " C: " << C<<std::endl;
		nonunifRandonPointsOnElipsoid(NbPoints,A,B,C,particles);
		BoxWith =  2.0*FMath::Max( A,FMath::Max( B,C)) ;
	}
	else if(FParameters::existParameter(argc, argv, "-plummer")){
		const FReal Radius  = FParameters::getValue(argc,argv,"-radius",  10.0);
		unifRandonPlummer(NbPoints, Radius, sum, particles) ;
		BoxWith = 2.0*Radius ;
	}

	else {
		std::cout << "Bad geometry option"<< std::endl;
		exit(-1) ;
	}
    /////////////////////////////////////////////////////////////////////////
	//                                           Save data
    /////////////////////////////////////////////////////////////////////////
	//
    //  Generate FMA file for FFmaGenericLoader<FReal> Loader
	//
	if(FParameters::existParameter(argc, argv, "-extraLength")){
		extraRadius  = FParameters::getValue(argc,argv,"-extraLength",  0.0);
		BoxWith += 2*extraRadius ;
	}
	std::string name(genericFileName);
	std::cout << "Write "<< NbPoints <<" Particles in file " << name << std::endl;
    FFmaGenericWriter<FReal>  writer(name) ;
	writer.writeHeader(Centre,BoxWith, NbPoints, *ppart) ;
	writer.writeArrayOfParticles(ppart, NbPoints);
	std::cout << "    End of writing "<<std::endl;

	//
	//  Generate  file for visualization
	//
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputVisuFile.options)){
        std::string visufile(FParameters::getStr(argc,argv,FParameterDefinitions::OutputVisuFile.options,   "output.vtp"));
         driverExportData(visufile, particles , NbPoints);
	}
	//
	delete [] particles ;

	//
	return 1;
}
