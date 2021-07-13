#include "tools.hpp"

#include <vector>
#include <cmath>
#include <iostream>

/*
template<typename Number>
int computeRunningAverage (Number addend, Number & RunningMeanValue, Number & RunningSqMeanValue, size_t nValues){
	
		
		
   if (nValues == 0) {
      RunningSqMeanValue += addend;
   }
   else {
      RunningMeanValue = RunningMeanValue * nValues/(nValues+1) + addend/(nValues+1);
   }
    
   //updating mean 
   if (nValues == 0) {
      RunningMeanValue += addend;
   
   }
   else {
      RunningMeanValue = RunningMeanValue * nValues/(nValues+1) + addend/(nValues+1);
   }
   
   return nValues+1;
}
*/
double statUncertainty(double mean, double sqMean, size_t N){
	if (N <= 1) return 0.;

	return sqrt((sqMean-mean*mean)/(N-1));
}

//start BlockMean member functions
BlockMean :: BlockMean( size_t m , size_t size ) :val(size, 0.), mean(size, 0.), sqMean(size, 0.), n(0), nPerBlock(m) {}
BlockMean :: BlockMean(size_t m): BlockMean (m, 1) {}
BlockMean :: ~BlockMean(){}

void BlockMean :: Reset (size_t m, size_t size) {
   val = std::vector <double> (size, 0.);
   mean = std::vector <double> (size, 0.);
   sqMean = std::vector <double> (size, 0.);
   n = 0;
   nPerBlock = m;
}

void BlockMean :: AddValue (double * Values, size_t nvals) {
   for(size_t i = 0; i < nvals; ++i){
      val[i] += Values[i];
   }
}

void BlockMean :: AddMean() {
   for(auto & el : val){
      el *= 1./nPerBlock;
   }
   for(size_t i = 0; i < val.size(); ++i){
      mean[i] = mean[i] * (n/double(n+1)) + val[i] * (1./(n+1));
      sqMean[i] = sqMean[i] * (n/double(n+1)) + val[i]*val[i] * (1./(n+1)); 
   }
   ++n;
   for(auto & el : val){
      el = 0.;
   }
}  

double BlockMean :: GetMean( size_t i ){return mean[i]; }
double BlockMean :: GetSqMean( size_t i) {return sqMean[i];}

//end BlockMean memeber functions



template <typename NumberH, typename NumberV>
double chiSquared(const std::vector<NumberH> &histogram, const std::vector<NumberV> &meanValues  ){
	
	if (histogram.size() != meanValues.size()) {
		std::cerr << "error: 'histogram' and 'meanValues' must have same size\n";
		return -1.;
	}

	double Chi=0;
	auto itv = meanValues.begin();
	auto ith = histogram.begin(); 	

	while (ith != histogram.end()){
		
		if (*itv < 0)  {
			std::cerr << "error: 'meanValues is not a probability distribution\n";
			return -2;
		}
		Chi += pow(*ith - *itv, 2) / *itv; 

		++ith;  ++itv;
	}
	return Chi;
}

void computeBlockAverage(const std::vector<double>& blocks, std::vector<double> & means, std::vector<double>& uncertainties){

   double Mean = 0;
   double sqMean = 0;
   means.clear();
   uncertainties.clear();
	   
   for (size_t i = 0; i < blocks.size(); ++i){
     computeRunningAverage (blocks[i], Mean, sqMean, i);
     means.push_back(Mean);
     uncertainties.push_back(statUncertainty(Mean, sqMean, i));
   }  
}


      
