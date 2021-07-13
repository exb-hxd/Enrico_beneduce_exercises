#include "tools.hpp"
#include "Metropolis.h"
#include <fstream>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include <cstdlib>
#include <filesystem>


using namespace std;

double Module (const vector <double>& vec){ 
   return sqrt( accumulate (vec.begin(), vec.end(), 0.0, [](const double & t1, const double & t2) {return t1 + t2*t2;} ));
};

double CosTheta (const vector <double>& vec){
   if (vec.size() != 3) throw "theta: argument must be of size() = 3";
   return vec[2]/Module(vec);
}

int main(int argc, char** argv){



   auto orbital100Squared = [](const vector<vector<double>> & vec) -> double {
      return 1./M_PI * exp(- Module(vec[0])*2); 
   };
 
   auto orbital210Squared = [](const vector<vector<double>> & vec) -> double {
      return (1./(32. * M_PI)) * exp (- Module(vec[0])) * std::pow( Module (vec[0]) * CosTheta(vec[0]), 2);
   };
   
   string names[] = {"100uniform", "210uniform", "100Gauss", "210Gauss"};


//create directory for output and set parameters for output
  std::string dirname = "output/dir05";
  if (!filesystem::exists (dirname)) {
     cerr << "creating directory \"dir05\"";
     system (("mkdir " + dirname).c_str() );
  }
  unsigned int wd = 20;

//parameters for equilibration
   size_t nMaxEquilibrate = 50000;

   vector <vector <vector <double> > > startPoints {
      {1, vector <double>(3,0.)},
      {1, vector <double>(3,1.)},
      {1, vector <double>(3,2.)},
      {1, vector <double>(3,5.)},                      
   }; 


//creating equilibration data
for (auto const & pos : startPoints){ 

   std::map <string, BlockMean> fastRays {
      {names[0],1},
      {names[1],1},
      {names[2],1},
      {names[3],1},
   };


//creating output
   stringstream formatter;
   formatter << dirname << "/05-1.eq_x" << setprecision (1) << pos[0][0] << 'y' 
                                        << setprecision (1) << pos[0][1] << 'z'  
                                        << setprecision (1) << pos[0][2] << ".dat";
 
   ofstream output_eq (formatter.str()); 
   for (auto const & key : names){
      output_eq <<  setw(wd) << key;
   }
   output_eq << endl;

//initializing Metropolis instances
   map <string, class Metropolis> Samplers {
      {names[0],{3,1,1.,1.}},
      {names[1], {3,1,1.,1.}},
      {names[2],{3,1,1.,1.}},
      {names[3],{3,1,1.,1.}},
   };

   for (auto & [key, el ] : Samplers) {
      el.Start (pos);
   }

//perform equilibration 
   for (size_t n = 0; n < nMaxEquilibrate; ++n){
        Samplers.at(names[0]).Move(IN_RANGE, orbital100Squared);
        Samplers.at(names[1]).Move(IN_RANGE, orbital210Squared);
        Samplers.at(names[2]).Move(GAUSS, orbital100Squared);
        Samplers.at(names[3]).Move(GAUSS, orbital210Squared);
         
        for (auto & [key, el] : Samplers ){
 
           fastRays.at(key).AddValue (Module( (Samplers.at(key)).GetPos()->at(0) ));
           fastRays.at(key).AddMean();
           output_eq << setw (wd) << fastRays.at(key).GetMean();
           //uncomment lines above and comment following line if using with fastRays
//           output_eq << setw (wd) << Module ( (Samplers.at (key)).GetPos()->at(0));
        }
        output_eq << endl;
   }
}

//calculation of averages: initialization of Metropolis instances

//   string names[] = {"100uniform", "210uniform", "100Gauss", "210Gauss"};
     
   map <string, class Metropolis> Samplers {
      {names[0],{3,1,1.,1.}},     //__
      {names[1], {3,1,1.5,1.5}},  //__|
      {names[2],{3,1,0.8,0.8}},   //__|--"carefully reserached figures", chosen in order to keep acceptance about 50%
      {names[3],{3,1,1.,1.}},     //__|
   };
   
  
//   vector <vector <double>> startPoint (1, vector<double> (3, 1.)); 
   for (auto & [key,i] : Samplers){
      i.Start(startPoints [1] );
   }

   cout << "initial position:\n";
   for(auto & [key, i] : Samplers){
      cout << setw(10) << key; 
      auto vec = (i.GetPos())->at(0);
      std::for_each(vec.begin(), vec.end(), [](double val){cout << val << ' ';} );
      cout << endl;
   }


// variables for simulation
   size_t nBlocks = 100;
   size_t nsteps  =1000;


   std::map <std::string, size_t> nEquilibrate { //"carefully researched figures": chosen after looking at the graphs of instantaneous averages
      {names[0],2000},
      {names[1], 1000},
      {names[2], 3000},
      {names[3], 2000},
   };

   map <string, BlockMean> rays {
      {names[0], nsteps},
      {names[1], nsteps},
      {names[2], nsteps},
      {names[3], nsteps}
   };


  

//make out output file
   ofstream output(dirname + "/05-1.dat");
   for (auto const & [key , i] : rays){
      output <<  setw(wd) << key << setw(wd) << "error" << setw (wd) << "acceptance";
   }
   output << endl;


// discard out-of-equilibrium values
   for (size_t nDiscard = 0; nDiscard < nEquilibrate ["100uniform"]; ++nDiscard){
      Samplers.at("100uniform").Move(IN_RANGE, orbital100Squared);
   }
   for (size_t nDiscard = 0; nDiscard < nEquilibrate ["210uniform"]; ++nDiscard){
      Samplers.at("210uniform").Move(IN_RANGE, orbital210Squared);
   }
   for (size_t nDiscard = 0; nDiscard < nEquilibrate ["100Gauss"]; ++nDiscard){
      Samplers.at("100Gauss").Move(GAUSS, orbital100Squared);
   }
   for (size_t nDiscard = 0; nDiscard < nEquilibrate ["210Gauss"]; ++nDiscard){
      Samplers.at("210Gauss").Move(GAUSS, orbital210Squared);
   }

// compute averages
   for (size_t i = 0; i < nBlocks; ++i){
      for (size_t j = 0; j < nsteps; ++j){
         Samplers.at(names[0]).Move(IN_RANGE, orbital100Squared);
         Samplers.at(names[1]).Move(IN_RANGE, orbital210Squared);
         Samplers.at(names[2]).Move(GAUSS, orbital100Squared);
         Samplers.at(names[3]).Move(GAUSS, orbital210Squared);

         for (auto & [key, el] : rays){
            el.AddValue(Module( (Samplers.at(key)).GetPos()->at(0) ));
         }          
      }
      for(auto & [key, el] : rays){
         el.AddMean();
         output << setw(wd) << el.GetMean() << setw(wd) << statUncertainty(el.GetMean(), el.GetSqMean(), el.GetN()) << 
                   setw(wd) << Samplers.at(key).GetAcceptRate();
      }
      output << endl;
   }
  

return 0;
}
