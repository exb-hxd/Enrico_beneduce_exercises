#include "Metropolis.h"
#include "tools.hpp"
#include "armadillo"
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <sstream>
#include <map>


using namespace std;


inline double my_square(double x){ //substitutes std::pow( , 2)
   return x*x;
}

 
inline double trial_ground_state(double mu, double sigma, double x) {                                     
   return  ( std::exp(-my_square((x-mu)/sigma)/2.) +   //ground state wave function (its square gives the
             std::exp(-my_square((x+mu)/sigma)/2.)  ); //probability of measuring the particle in x)
}

inline double potential (double x){ //potential energy at position x
   return my_square(my_square(x)) - 5./2. * my_square(x);
}

inline double T_acts_on_psi(double mu, double sigma, double x){ //kinetic energy functional acting on psi
   return -0.5 * ( std::exp(-my_square((x-mu)/sigma)/2.) * ( my_square((x-mu)/sigma) - 1. ) + 
                  std::exp(-my_square((x+mu)/sigma)/2.) * ( my_square((x+mu)/sigma) - 1. )   ) /my_square(sigma);
}


struct Energy_computer { 
   std::string filename;

   Energy_computer ( const std::string& f) : filename (f) {}

   double operator () (const  std::vector<double>  & vec) {   

      double Mu = vec.at(0);
      double Sigma = vec.at(1);

      Metropolis met(1, 1, Sigma , Sigma); //initially, the steps are about the size of Sigma, which should be a reasonable value

      
      auto acceptance = [&] (const std::vector <std::vector <double>> & vec ) -> double {
                           return my_square (trial_ground_state (Mu, Sigma, vec[0][0]) );  
                        };

      try{
         met.Start({{Mu} } ); //start from a value that should have hig probability
         met.FixStep (GAUSS, acceptance, 1000, 10, 0.1);
      } 
      catch(std::string str){
         cout << str << endl;
      } 

      size_t nPerBlock = 10000;
      size_t nBlocks = 50;
      size_t wd = 20;
      enum quantities {ENERGY = 0, ACCEPTANCE  = 1};
      size_t nqt = 2;
      BlockMean qt ( nPerBlock, nqt);


      std::ofstream outputMet (filename);
      std::ofstream outputPoints ( "output/dir08/points.dat");


      if (!outputMet.good()) cout << "error opening output file";
      
      for (size_t i = 0; i < nBlocks; ++i){
         for(size_t j = 0; j < nPerBlock; ++j){
            met.Move(GAUSS, acceptance);
            double x = (*met.GetPos()) [0][0];
            double values [] = { [ENERGY] = ( potential (x) + T_acts_on_psi (Mu, Sigma,  x) ) / 
                                                  std::abs (trial_ground_state (Mu, Sigma, x) ),
                                 [ACCEPTANCE] = met.GetAcceptRate()};
            qt.AddValue( values, nqt );
            outputPoints << x << "  ";
         }      
         qt.AddMean();
         outputPoints << endl;
         outputMet << qt.GetMean(ENERGY) << setw (wd) << ' '  << statUncertainty (qt.GetMean(ENERGY), qt.GetSqMean(ENERGY), qt.GetN()) 
                   << setw(wd)  << ' ' << qt.GetMean(ACCEPTANCE)  << std::endl;
      }
//      outputMet << qt.GetMean(ENERGY) << setw (wd)  << statUncertainty (qt.GetMean(ENERGY), qt.GetSqMean(ENERGY), qt.GetN()) 
  //              << setw(wd) << qt.GetMean(ACCEPTANCE) << std::endl;
//      outputMet << "EndVal "<< endl;
      return qt.GetMean(ENERGY);
   }
};


int main (int argc, char **argv) {

   std::string dirname("output/dir08");
   std::string strNiterate = (argc >= 2 ? argv [1] : "100" );
   std::string strReduce = (argc >= 3 ? argv [2] : "20");

   if(filesystem::exists(dirname) == false){
      system(("mkdir " + dirname).c_str());
   }
   
   std::stringstream outname ;
   outname << dirname << "/Energy_and_acceptance.dat";
   string energyFilename (outname.str());

   if(filesystem::exists(energyFilename)) {
      system ( ( "rm " + energyFilename ).c_str() );
   }
   if (filesystem::exists (dirname + "/points.dat")){
      system (("rm " + dirname + "/points.dat").c_str()); 
   }


 
   double startMu = 0.5;
   double startSigma = 0.7;


   ofstream outCheck (dirname + "/trial_wave.dat");
   for (double i = -2; i < 2; i += 0.1){
      outCheck << trial_ground_state(startMu, startSigma, i) << ' ' 
               << T_acts_on_psi(startMu, startSigma,  i) << ' ' << potential(i) <<  endl; 
   }

   
  

   outname.str("");
   outname << dirname << "/Mu_Sigma.dat";
   
   ofstream output (outname.str());
   if (! output.good()){
      cout << "unable to open file" << outname.str();
      return 1;
   } 

 
   Energy_computer Energy (energyFilename);

   size_t wd = 2;

   double CenterMu = 1.;
   double SpanMu = 1.;

   double CenterSigma = 1.;
   double SpanSigma = 0.95; //sigma should not be allowed to go too close to 0

   arma::uword npoints = 10;


   arma::mat grid  (npoints, npoints);
   
   unsigned int  nIterate;
   double reduceFactor;
   try {
      nIterate = std::stoul (strNiterate);
   }
   catch (std::invalid_argument & e){
      std::cerr << e.what();
      nIterate = 10;
   }
   try {
      reduceFactor = std::stod (strReduce);
   }
   catch (std::invalid_argument & e){
      std::cerr << e.what();
      reduceFactor = 100;
   }
   double restrict = std::pow (static_cast<double> (reduceFactor) , 1./ static_cast <double> (nIterate) );


   for (unsigned int iIterate =0; iIterate < nIterate ; ++ iIterate){
      arma::vec mucoord = arma::linspace (CenterMu - SpanMu, CenterMu + SpanMu, npoints);
      arma::vec scoord = arma::linspace (CenterSigma - SpanSigma, CenterSigma + SpanSigma, npoints);

    
     
      for (arma::uword i = 0; i < npoints; ++i){
         for (arma::uword j = 0; j < npoints; ++j){
            grid (j,i) = Energy ({mucoord(j) , scoord (i)}) ;
         }
      }

//      std::cout << grid;
      
      arma::Col <arma::uword>  findvec = arma::find( grid == grid.min());
      arma::uword j = findvec[0] % npoints,
                  i = (findvec[0] / npoints) % npoints;

      double min = grid.min();

      std::cout << " i " << i << " j " << j << " min " << grid.min() << " foundmin " << grid (j,i) << std::endl;

      CenterMu = mucoord (j);
      CenterSigma = scoord (i);
      SpanMu /= restrict;
      SpanSigma /= restrict;      

      
      std::cout << "center energy " << min << " Mu " << CenterMu <<" Sigma " << CenterSigma << std::endl;      
      output << CenterMu << ' '  << CenterSigma << ' ' << Energy ({CenterMu, CenterSigma}) << std::endl ; //it is VERY important that Energy() is called 
   }                                                                                                      //in this line insetead of reusing "min" value. 
                                                                                                          // This way the energy and point files are updated so
                                                                                                          // that they refer to the optimum parameters and not to
                                                                                                          // to those in the rightmost - downmost cell in "grid" 

   return 0;
}
