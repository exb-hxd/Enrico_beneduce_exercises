#include "random.hpp"
#include "Metropolis.h"

#include <algorithm>
#include <exception>
#include <iostream>
#include <vector>

//start Metropolis definitions
void Metropolis :: Propose(Distribution which){
   if (which != IN_RANGE){
      if(which != GAUSS){
         throw "Metropolis: distribution is not GAUSS nor UNIFORM";
      }
   }
  
   for (size_t i = 0; i < ncoordinates; ++i){
      posnew [point] [i] = pos [point] [i] + Cumulate (which, 1);//(which == GAUSS) ? Gauss() : RannyuInRange(); // gives 1 value sampled from given distribution 
   }
   point = (point + 1) % npoints;
//   if (npoints > 1 ) std::cerr << point;
   if (point > npoints) throw std::length_error("npoints exceeded");
   ++ attempted;
}



Metropolis :: Metropolis(size_t nc, size_t np, double _sigma, double _delta, int nseed): ncoordinates(nc), 
npoints(np), point(0), pos(np, std::vector<double> (nc, 0)), posnew(pos), probold(-1) , accepted (0), attempted (0), saveStart (pos) { 

//   std::cerr << npoints;

   pos.shrink_to_fit();
   pos.shrink_to_fit();

   SetParameter(MEAN, 0);
   SetParameter(SIGMA, _sigma);
   SetParameter (MIN, -_delta);
   SetParameter (MAX, _delta);

   DoAll(nseed);
} 

//Metropolis :: Metropolis (size_t nc, size_t np, double _sigma, double _delta) : Metropolis (nc, np, _sigma, _delta, 1) {}

Metropolis :: ~Metropolis(){}
Metropolis :: Metropolis(const Metropolis& other) : Metropolis( other.ncoordinates, other.npoints, other.GetParameter(SIGMA), other.GetParameter(MAX)) {}

Metropolis& Metropolis :: operator= (const Metropolis& other) {
   ncoordinates = other.ncoordinates;
   npoints = other.npoints;
   point = other.point;
   pos = other.pos;
   posnew = other.posnew;
   probold = other.probold;
   attempted = other.attempted;
   accepted = other.accepted;
   saveStart = other.saveStart;
   SetParameter(MEAN, other.GetParameter(MEAN));
   SetParameter(SIGMA, other.GetParameter(SIGMA)); 
   SetParameter(MIN, other.GetParameter(MIN));
   SetParameter(MAX, other.GetParameter(MAX));
   return *this;
}



bool operator< (const Metropolis & lhs, const Metropolis & rhs) {
   if (lhs.npoints < rhs.npoints) return true;
   else return false;
}
 


void Metropolis :: Start(const std::vector<std::vector<double>>& pstart) {

   if(pstart.size() != npoints) {
      throw std::string("Metropolis :: Start() : wrong number of points");
      return;
   }
     
   for(auto vec : pstart){
      if( vec.size() !=  ncoordinates){
         throw "Metropolis :: Start() : wrong dimension"; 
         return;
      }
   }
   pos = pstart;
   saveStart = pstart;
   probold = -1; //probold is assigned value -1 in order to make Accept(functor) evaluate 
                 //it the first time without having to pass functor to Start()

   accepted = 0;
   attempted = 0; 
}

void Metropolis :: SetDelta ( double Value ){
   Value = std::abs( Value );

   SetParameter (MIN, -Value);
   SetParameter (MAX, Value);
}

  
//end Metropolis definition


