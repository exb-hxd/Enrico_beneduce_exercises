/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.hpp"
#include <cmath> 
#include "tools.hpp"



template<typename Number>
int computeRunningAverage (const Number &addend, Number *RunningMeanValue,  unsigned int nValues, Number *RunningSqMeanValue){
	
		//updates mean of square
		if(RunningSqMeanValue){	
				
			if (nValues == 0) {
				*RunningSqMeanValue += addend * addend;
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




using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   rnd.DoAll();


// we create N blocks of M/N throws each


   unsigned int N = 100;
   unsigned int M = 100000;
   unsigned int nsteps = M/N;


 //later on we may try implement a smarter way of addressing elements in \"means\"

   double Mean = 0;
   double sqMean = 0;

   double Sigma = 0;
   double sqSigma = 0;


   ofstream output("output/dir01-1/01-1a");
   if(!output.is_open()) cerr << "unable to create output file"; 


   //loop over blocks 
   for (unsigned int i =0; i < N; ++i){
   	double blockSum =0; 
	double blockSqSum=0;

        // loop over single block and compute sum and square sum (faster than running averaging)	
	for(unsigned int n =0; n < nsteps; ++n){
	   double value = rnd.Rannyu();
	   blockSum += value;
	   blockSqSum += pow(value,2);

	}
	double blockMean = blockSum/nsteps;
	double blockSigma = blockSqSum/nsteps - pow(blockMean, 2);

	//update 'Mean' and 'Sigma' values to current block
	computeRunningAverage( blockMean, &Mean, i, &sqMean);
	computeRunningAverage( blockSigma, &Sigma, i,  &sqSigma);


	output << i << ' ' << Mean-0.5 << ' ' << statUncertainty(Mean, sqMean, i) << ' ' <<Sigma - 1./12 << ' ' << statUncertainty(Sigma, sqSigma, i) << endl;

   }
  
 

   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
