#include <iostream>
#include <cmath>
#include "random.hpp"
#include <fstream>
#include "TVector2.h"
#include <vector>




template<typename Number>
int computeRunningAverage (const Number addend, Number *RunningMeanValue,  unsigned int nValues, Number *RunningSqMeanValue){
	
		//updates mean of square
		if(RunningSqMeanValue){	
//	std::cout << *RunningSqMeanValue << std::endl;
			if (nValues == 0) {
				*RunningSqMeanValue += addend*addend;
//				std::cout << *RunningSqMeanValue;
			}
			else {
				*RunningSqMeanValue = ((*RunningSqMeanValue) * nValues/(nValues+1.)) + addend*addend/(nValues+1.); 
//				std::cout << *RunningSqMeanValue;
			}
		}

		//updating mean (always done)
		if (nValues == 0) {
			*RunningMeanValue += addend;
		//	++nValues;
		}
		else {
			*RunningMeanValue = *RunningMeanValue * nValues/(nValues+1.) + addend/(nValues+1.);
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
     output<< i << ' ' << Mean << ' ' << statUncertainty(Mean, sqMean, i) << std::endl; 
      
   
   }  


}


//if point P =(xP,yP) lies in the unitary circle, updates x and y to xP and yP respectevely. 
//xP and yP are sampled from uniform distribution 
void rejectInCircle (double &x, double &y, Random & rnd ){
   double xP= 1; double yP = 1;

   while (sqrt (xP*xP + yP*yP) > 1) {
      xP = rnd.Rannyu();
      yP = rnd.Rannyu();
   }

   x = xP, y = yP; 
}

TVector2 getVersor (double x, double y){
   TVector2 vec (x, y);
   return vec.Ort();
}

class Needle {
      double length, span, xCenter, yCenter;
      TVector2 vector;
   
     
   public:
      Needle(double L, double d) : length(L), span(d), xCenter(0), yCenter(0), vector (0, length)  {}
      ~Needle(){} 

      void GetVector (double x, double y) {
         vector.SetX(x); vector.SetY(y);
	 vector = vector.Unit() * length/2;
 
      }   
   
      void GetCenter(double x, double y){
         xCenter = x; yCenter = y;

      }


      int DoesIntersect() {
         if( abs(std::floor( (xCenter - vector.X())/span ) - std::floor( (xCenter + vector.X() )/span ) ) >= 1  ) return 1;
	 else return 0;
      }
      
   
};


double getBuffonPiEstimate(double L , double d, double statProb){
   return 2 * L / (d * statProb);

}


using namespace std; 



int main (int argc, char **argv){
   
//creating N blocks, nsteps throws each; total number of throws: M	
   const size_t N = 100;
   const size_t M = 10000000;
   const size_t nsteps = M/N;

// choosing parameters for the needle/ the lines:
   const double Length = 4.5;
   const double Distance = 5;


// create and initialize random number generator
   Random rnd;
   rnd.DoAll();
   rnd.SetParameter(MIN, - Distance/2);
   rnd.SetParameter(MAX, Distance/2);

//creating vector to store block results ... and needle to compute them!
   vector<double> PiValues;
   Needle  ned (Length, Distance);

   for (size_t j = 0; j < N; ++j){

      unsigned int nintersect = 0;

      for (size_t i = 0; i < nsteps; ++i ){
   
	 //sort x coordinate for first extreme point
         ned.GetCenter(rnd.RannyuInRange(), 0);

	 //sort needle direction
         double x, y;
         rejectInCircle (x, y, rnd);
         ned.GetVector (x, y);
   
         //check intersection
	 nintersect += ned.DoesIntersect();


      }

      double prob = double(nintersect)/nsteps; 

      PiValues.push_back(getBuffonPiEstimate(Length, Distance, prob)); 
   }
 //    for (auto it = PiValues.begin(); it != PiValues.end(); ++it ){
   //     cout << *it << endl;
    // }
 
   string filename = "output/01-3";
   computePrintBlockAverage(PiValues, filename);

   //appending information about parameters used
   ofstream moreOutput;
   moreOutput.open(filename.c_str(), ios_base::app);
   moreOutput << "Length: " << Length << " Distance: " << Distance << endl ;  


return 0;
}
