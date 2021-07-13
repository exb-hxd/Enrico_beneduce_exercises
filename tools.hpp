
#ifndef __Tools__
#define __Tools__

#include <random.hpp>
#include <iostream>
#include <cmath>
#include <vector>

template<typename Number>
int computeRunningAverage (Number addend, Number & RunningMeanValue, Number & RunningSqMeanValue, size_t nValues){
	
		
		
   if (nValues == 0) {
      RunningSqMeanValue += addend*addend;
   }
   else {
      RunningSqMeanValue = RunningSqMeanValue * nValues/(nValues+1.) + addend*addend/(nValues+1.);
   }
    
   //updating mean 
   if (nValues == 0) {
      RunningMeanValue += addend;
   
   }
   else {
      RunningMeanValue = RunningMeanValue * nValues/(nValues+1.) + addend/(nValues+1.);
   }
   
   return nValues+1;
}

template <typename NumberH, typename NumberV>
double chiSquared(const std::vector<NumberH> &, const std::vector<NumberV> &  );



double statUncertainty(double, double, size_t );

void computeBlockAverage(const std::vector<double> &, std::vector<double> &, std::vector<double> &);

class BlockMean {
      std::vector <double> val, mean, sqMean;
      size_t n, nPerBlock;

   public: 
      BlockMean( size_t, size_t );
      BlockMean( size_t );
      ~BlockMean();

     inline void AddValue (double v){
        AddValue (&v, 1);
     }
     void AddValue (double *, size_t);

     void AddMean();  

     inline double GetMean(void){
        return GetMean (0);
     }
     inline double GetSqMean(void){
        return GetSqMean (0);
     }
     inline size_t GetN(void){
        return n;
     }

     double GetMean(size_t);
     double GetSqMean(size_t);

     void Reset (size_t, size_t);
};

#endif
