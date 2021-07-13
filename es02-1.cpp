
#include <iostream>
#include "random.hpp"
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "TF1.h"


template<typename Number>
int computeRunningAverage (const Number &addend, Number *RunningMeanValue,  unsigned int nValues, Number *RunningSqMeanValue){
	
		//updates mean of square
		if(RunningSqMeanValue){	
				
			if (nValues == 0) {
				*RunningSqMeanValue += addend*addend;
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



void computePrintBlockAverage(const std::vector<double>& blocks, const std::string & filename  ){

   std::ofstream output(filename);

   double Mean = 0;
   double sqMean = 0;
	   
   for (size_t i = 0; i < blocks.size(); ++i){
     computeRunningAverage (blocks[i], &Mean, i, &sqMean );
     output << i << ' ' << Mean << ' ' << statUncertainty(Mean, sqMean, i) << std::endl; 
      
   
   }  


}


/*
struct ImportanceSampler {

//inverts uniform distribution using the inversion function f 
   double ImportanceDist (double y ){
      return 1 - sqrt (1-y);
   }
// function we want to integrate   
   double integrand (double x){
      return M_PI/2 * std::cos(x) 
   }
   

};

class PDFinverter{
      TF1 PDF;
      double normFactor;

   public

       PDFinverter(double intMin, intMax, TFormula f):  PDF("PDF", f, intMin, intMax) {

	 PDF.SetNpx(4);
	 if (unnormPDF.GetMinimum(intMin, intMax) < 0 ) {
	    throw "error: negative value in PDF" 
	 }
         if(PDF.Integral() != 1) {
	     throw "error: unnormalized PDF"
	 }
	  


      
      }

};
*/

double inversionFunc (TF1* f, double val, double min, double max){

// underliyng line may be edited according to sampling distribution
//   double result = 1 - sqrt(1 - val);
//

   
      	
   double result = f->Eval(val);

   if (result < min or result > max) { 
      std::cerr << "error: out of range";
      throw 1;
   } 
   return result;
}

using namespace std;

int main (int argc, char **argv){
   
   Random rnd;
   rnd.DoAll();

   const size_t M = 2000;
   const size_t N = 100;

   vector <double> Integrals, ImportanceIntegrals;
   for (size_t n = 0; n < N; ++n){
      Integrals.push_back(0);
      ImportanceIntegrals.push_back(0);
   }

   TF1 PDF ("d", "2-2*x", 0, 1);
   TF1 inversePDF ("d-1","1-sqrt(1-x*x)", 0, 1);
   //uniform ditribution
   for (size_t i = 0; i < N; ++i){
      for (size_t j = 0; j < M; ++j){
	 //compute \"random\" number      
	 double value = rnd.Rannyu();



	 //increase integrals
         Integrals[i] += M_PI/2*cos(M_PI/2*value);
	 double ImportanceX = inversionFunc(&inversePDF, value, 0, 1);
	 ImportanceIntegrals[i] += M_PI/2*cos(M_PI/2*ImportanceX)/PDF.Eval(ImportanceX);

   
      }
      //oops... we should have taken the average over M throws of the previous result... let's mend:
        Integrals[i] *= 1./M;
        ImportanceIntegrals[i]  *= 1./M;
   
   }
   
   //print output to file
   string rootname = "output/dir02-1/";
   string filename = rootname + "Uniform";

   computePrintBlockAverage (Integrals, filename);

   filename.clear();
   filename = rootname + "Importance";

   computePrintBlockAverage (ImportanceIntegrals,  filename);

return 0; 
}
