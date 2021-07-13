#include "random.hpp"
#include "Genetic.h"
#include "circle.h"
#include "armadillo"

#include <algorithm>
#include <iostream>
#include <fstream>
//#include <cmath>
#include <string>
#include <vector>
#include <filesystem>


arma::mat makeCircle (Random & rnd, arma::uword ncities){ 
     PointOnCircle circle (1.);
   
     arma::mat cities = circle.randCoord (rnd); 
     for (arma::uword icity = 1; icity < ncities; ++icity ) {
   
        cities = join_rows (cities, circle.randCoord (rnd));
     }
     return cities;
}

arma::mat makeSquare (Random & rnd, arma::uword ncities){
     arma::mat cities  (2, ncities);
     cities.imbue ( [&] ()->double {return rnd.Rannyu();} );
     return cities;
}

namespace fs =  std::filesystem;

int main (int argc, char ** argv) {

  std::string args [3];
  std::string defaults [] = {"square", "10", "no"};

  for (int i = 1; i < argc; ++i ){
     args [i - 1] = argv[i];
  }
  for (int i = argc; i < 4; ++i){
     args [i - 1] = defaults [i - 1];
  }

  std::string name = (args [0] == "circle" ? "circle" : "square" );
  int nGen = std::stoi (args [1]);
  bool doSaveMeans = (args [2] == "yes"); 


//creating output directory
  std::string dirname ("output/dir09");
  if ( !std::filesystem::exists (dirname)){
     std::system (("mkdir " + dirname).c_str());
  }
  
  std::string filename (name + "_fitness.dat"); 

  { //limiting scope of init object, which is needed no more later
     std::ofstream init (dirname + '/' + filename );
     init << "p1 p2 p3 length sequence\n";
  }

// creating file for writing convergence avarages
  std::ofstream means (dirname + '/' + name + "_progressive_fitness.dat");
  means << (doSaveMeans ? "yes" : "no") << std::endl;


// reading input files with parameters for population size, p1, p2, p3, pCross, nGenerations
  std::map <std::string, double> Parameters {
     {"ncities", 32},
     {"nCouples",450},
     {"p1",0.1},
     {"p2",0},
     {"p3",0},
     {"pCross",0.5},
     {"maxNgenerations", 100} 
  };

   for (const auto & entry : fs::directory_iterator(std::filesystem::path (dirname) / "input_files")){
      std::cout << entry;
   
     std::ifstream input (entry.path().string());   
     
   
     {  // just a scope in order to make key and value disappear after use
           std::string key;
           double value;
           std::string readstr;
   
           while (std::getline (input, readstr)) {
              std::stringstream s (readstr);
              s >> key >> value;
              std::cout << readstr;
              if (input.eof()) break; 
              Parameters. at (key) = value;
        }
     }

     for (const auto & [key, el] : Parameters ){
        std::cout << key << ": " << el << std::endl;
     }
 
   
     Random rnd;
     rnd.DoAll(10);
   
     arma::uword ncities = Parameters ["ncities"];
     
     arma::mat cities = (name == "square" ? makeSquare (rnd, ncities) : makeCircle (rnd, ncities) ); 

     std::cout << cities;
   
   
   
   
   //file for saving point coordinates
     std::ofstream WriteCities (dirname + '/' + name + "_points.dat");
   
     WriteCities << cities;
   
     FitnessComputer environment ( cities );
   
     size_t nCouples = Parameters ["nCouples"];
   
     unsigned int maxNgenerations = Parameters ["maxNgenerations"];
     Population pop (nCouples, ncities, rnd );
     std::vector <double> mean_lengths;
     std::vector <CodePerm> BestIndividuals ;
   
    
   
     std::vector <double> probabilities  {Parameters ["p1"], 
                                          Parameters ["p2"], 
                                          Parameters ["p3"],
                                          Parameters ["pCross"]
                                         };
       
     
     std::ofstream output (dirname + '/' + filename, std::ios::app);
                                                                                             
     for (unsigned int iGen = 0; iGen < maxNgenerations; ++iGen){
        pop.NextGen (environment, rnd, probabilities);
   //      mean_lengths.push_back (pop.MeanOfBestHalf (environment));
         BestIndividuals.emplace_back (pop.GetBest (environment));
   //      BestIndividuals.back().Print (3);
          
         if (iGen % nGen == 0) {
            std::cout << "before sorting: ";
            BestIndividuals.front().Print (3); 
            std::sort (BestIndividuals.begin(), BestIndividuals.end(), environment);
            std::cout << "after sorting: ";
            
            if ( doSaveMeans ){
               means << pop.MeanOfBestHalf(environment) << ' ' << BestIndividuals.front().GetLength() << std::endl; 
            }

            BestIndividuals.front().Print (3);    
            BestIndividuals.erase ( BestIndividuals.begin() + 1, BestIndividuals.end());
         }
     }

     
     output << Parameters ["p1"] << ' ' 
            << Parameters ["p2"] << ' ' 
            << Parameters ["p3"] << ' '  
            << Parameters ["pCross"] << ' ' 
            << BestIndividuals [0].GetLength ();


     for (auto const & ind : BestIndividuals [0].GetFeno() ){
        output << ' ' << ind;
     }
     output << std::endl; 
   }
  return 0;
}
