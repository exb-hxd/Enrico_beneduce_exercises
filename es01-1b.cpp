

#include <iostream>
#include <fstream>
#include <string>
#include "random.hpp"
#include <cmath> 
#include <vector>
#include "tools.hpp"
#include <algorithm>

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


using namespace std;

int main (int argc, char *argv[]){

   cout << "es01-1b: ... \n";

   Random rnd;
   rnd.DoAll();

   //we create a vector filled with expected values and another one 
   //to be filled with the result of simulation

   size_t N = 1000; //number of times we are going to compute chi squared
   size_t M = 100;  //number of interval for each histogram
   size_t nthrows = 10000;  //number of thrwos for each histogram
   double expected_value = double(nthrows)/M;


   vector<double> expected;
   vector<unsigned int> observed;
   for(size_t i = 0; i < M; ++i){
   	expected.push_back(expected_value);
	observed.push_back(0);
   } 
 
   //creating output file
   ofstream output("output/dir01-1/01-1b");
   
   if(!output.is_open()) cerr << "unable to create output file"; 
   
   
   for (size_t h = 0; h < N; ++h){

      //put each element of 'observed' to 0 
      for_each(observed.begin(), observed.end(), [] (unsigned int &n){n = 0;});

      //filling the histogram
      for (size_t i = 0; i < nthrows; ++i){
   	++ (observed.at( size_t ( rnd.Rannyu() * M)));     
      }
   
   
      
      //computing Chi squared:
      double Chi = chiSquared(observed, expected);
   
      if (Chi == -1 ) {
         cerr << "chiSquared: error -1 occurred";
         return -1;
      }
      if (Chi == -2){
         cerr << "chiSquared: error -2 occurred";
         return -2;
      }
      
      //printing result to file
      output  << Chi << endl; 

   }

   return 0;
}
