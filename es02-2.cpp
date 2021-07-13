#include <iostream>
#include <fstream>
#include "random.hpp"
#include <cmath>
#include "TVector3.h"
#include <cassert>
#include <vector>
#include "TVectorT.h"
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

template<typename Number>
int computeRunningAverage (const Number &addend, Number *RunningMeanValue,  unsigned int nValues, Number *RunningSqMeanValue){
	
		//updates mean of square
		if(RunningSqMeanValue){	
				
			if (nValues == 0) {
				*RunningSqMeanValue += addend;
			}
			else {
				*RunningSqMeanValue = *RunningSqMeanValue * nValues/(nValues+1) + addend*addend/(nValues+1);
			}
		}

		//updating mean (always done)
		if (nValues == 0) {
			*RunningMeanValue += addend;
		//	++nValues;
		}
		else {
			*RunningMeanValue = *RunningMeanValue * nValues/(nValues+1) + addend/(nValues+1);
		//	++nValues;
		}

		return 0;
}

double statUncertainty(const double &mean, const double &sqMean, const unsigned int N){
	if (N <= 1) return 0.;

	return sqrt((sqMean-pow(mean,2))/(N-1));
}

/*
template <typename Number>
void computeBlockAverage(const std::vector<Number>& blocks, std::vector <double> & means, std::vector <double> & sqMeans){

   means.clear(); 
   sqMeans.clear();

   double Mean = 0;
   double sqMean = 0;
	   
   for (size_t i = 0; i < blocks.size(); ++i){
     computeRunningAverage (blocks[i], &Mean, i, &sqMean );

     means.push_back(Mean); 
     sqMeans.push_back( statUncertainty(mean, sqMean, i) );
   }  

}
*/
//class for sorting the direction of following step
class RandomAngle3D {
      TVector3 versor;
      Random * pRnd;

   public: 
      RandomAngle3D( Random & random) : versor(1, 0, 0), pRnd (&random) {}
      ~RandomAngle3D( ) {}

      inline void SortAngle () {

         versor.SetTheta (acos (1 - 2 * pRnd->Rannyu() ));
	 versor.SetPhi(pRnd->Rannyu() * 2 * M_PI);
         versor.SetMag(1);
//	 versor.Print();
      }

      const TVector3* GetVersor() const {
         return  & versor;
      }
};

//class for moving around in 3D space/lattice space and tracking position
class walker3D {
      
      TVector3 position;
      int latticeX, latticeY, latticeZ;

      Random * pRnd;
      RandomAngle3D angle;

      

   public:
      walker3D(Random & random) : position (0, 0, 0), latticeX(0), latticeY(0), latticeZ(0), pRnd (&random), angle(random){}
      ~walker3D(){}

      
      void GetContinuousIncrement(){
         angle.SortAngle();
	 position += *(angle.GetVersor());
      }

      void GetLatticeIncrement(){
	     
	 int iCoordinate = std::floor(pRnd->Rannyu() * 3);
	 int increment = std::round(pRnd->Rannyu()) * 2 - 1;
	 
	 switch (iCoordinate){
	    case 0:
		 latticeX += increment;
		 break;
            case 1:
		 latticeY += increment;
		 break;
	    case 2:
		 latticeZ += increment;
		 break;
            default:
		 throw "invalid random index generated";
		 break;
	 }
      
      }

      void Restart(){
         position.SetMag(0);
         latticeX = 0;
	 latticeY = 0;
	 latticeZ = 0;
      }
      
      void PrintPosition (){
         position.Print();
      }

      void PrintLatticePosition () {
         std::cout << latticeX << ' ' << latticeY << ' ' << latticeZ << ' ' << std::endl; 
      }

      double GetSqDistance (){
         return position.Mag2();
      }

      double GetLatticeSqDistance (){
         return latticeX * latticeX + latticeY * latticeY + latticeZ * latticeZ;
      }

};


using namespace std;

int main (int argc, char **argv){

   Random rnd;
   rnd.DoAll();

   size_t nBlocks = 100;
   size_t walksPerBlock = 100;
   size_t stepsPerWalk = 100;

   TVectorT <double> continuousMeans (stepsPerWalk);
   TVectorT <double> continuousSqMeans (stepsPerWalk);
   TVectorT <double> contDisplacements (stepsPerWalk);


   TVectorT <double> latticeMeans (stepsPerWalk);
   TVectorT <double> latticeSqMeans (stepsPerWalk);
   TVectorT <double> lattDisplacements (stepsPerWalk);



   continuousMeans.Zero();
   continuousSqMeans.Zero();

   latticeMeans.Zero();
   latticeSqMeans.Zero();
   
   walker3D DrunkBird (rnd);


   //creating directory for storing 2 output files
   system ("mkdir output/dir02-2");
   ofstream outputCont("output/dir02-2/continuous");
   ofstream outputLatt("output/dir02-2/lattice");


   for(size_t i = 0; i < nBlocks; ++i){

      for (size_t j = 0; j < walksPerBlock; ++j) {        
	 //reset all 'work' variables     
         DrunkBird.Restart();
	 contDisplacements.Zero();
	 lattDisplacements.Zero();


	 for(size_t k = 0; k < stepsPerWalk; ++k){
            DrunkBird.GetContinuousIncrement();
            DrunkBird.GetLatticeIncrement();

	    contDisplacements [k] += DrunkBird.GetSqDistance();
            lattDisplacements [k] += DrunkBird.GetLatticeSqDistance (); 

	 }
      }
      
      for (size_t k = 0; k < stepsPerWalk; ++k ){
         //important: must divide by number of walks in order to obtain mean values
         computeRunningAverage (contDisplacements [k] / walksPerBlock, continuousMeans.GetMatrixArray() + k, i, continuousSqMeans.GetMatrixArray() + k);  
         computeRunningAverage (lattDisplacements [k] / walksPerBlock, latticeMeans.GetMatrixArray() + k, i, latticeSqMeans.GetMatrixArray() + k);  

      }
   }
 	 //write output to file   
   for (size_t k = 0; k < stepsPerWalk; ++k){

	 outputCont << k << ' ' << continuousMeans [k] << ' ' << statUncertainty (continuousMeans [k], continuousSqMeans [k], nBlocks) << ' ';
	 outputLatt << k << ' ' << latticeMeans [k] << ' ' << statUncertainty (latticeMeans [k], latticeSqMeans [k], nBlocks) << ' ';
	  
         outputCont << endl; outputLatt << endl;
   }

   return 0;
}
