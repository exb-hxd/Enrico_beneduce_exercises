//credits to https://hpc-tutorials.llnl.gov for the implementation of nearest neighbours exchange in ring topology


#include "random.hpp"
#include "Genetic.h"
#include "circle.h"
#include "armadillo"
#include "mpi.h"

#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;


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


int main (int argc, char ** argv) {


   //unpack arguments : circle or square, nMigrants, nGenerations, nSplit, 
  std::string args [4];
  std::string defaults [] = {"square", "4", "4", "30"};
  
  for (int i = 1; i < argc; ++i ){
     args [i - 1] = argv[i];
  }
  for (int i = argc; i <= 4; ++i){
     args [i - 1] = defaults [i - 1];
  }
  
  std::string name = (args [0] == "circle" ? "circle" : "square" );
 
  const int nMigrants = std::stoi (args [1]);
  unsigned int nGenerations = std::stoul (args [2]); //number of internal generations (= number of traditional GA moves in each loop)
  unsigned int nSplitLoop = std::stoul (args [3]); //maximum number of loops (= number of times individuals are exchanged between processees)

  const int buffer_size = 60;

  int nprocess, rank, next, prev, tag = 1;

  unsigned long int ** sendbuf ,  ** recbuf; //dynamic memory alloaction is necessary in order to 
                                        //control the parameter nMigrants from commandline (argv[1])
  sendbuf = new unsigned long int * [nMigrants];
  recbuf = new unsigned long int * [nMigrants];

  sendbuf [0] = new unsigned long int [nMigrants * buffer_size];
  recbuf [0] = new unsigned long int [nMigrants * buffer_size];

  for (int iMig = 1; iMig < nMigrants; ++iMig) {
     sendbuf [iMig] = new unsigned long int [buffer_size];
     recbuf [iMig] = new unsigned long int [buffer_size];
  }  //end dynamic allocation 


  unsigned long int  collectbuf [4] [buffer_size];

  bool endloop = false; //terminating condition   


  Random rnd_forPoints, rnd;
  //each process must use a different sequence of random numbers for mutation, but the same one for generating the set of points
  rnd_forPoints.DoAll(10);

  arma::uword ncities = buffer_size / 2 + 2;
  
  arma::mat cities = (name == "circle" ? makeCircle (rnd_forPoints, ncities) : makeSquare (rnd_forPoints, ncities));
  
  
  FitnessComputer environment ( cities );
  
  //FitnessComputer env1 (environment);  //if a bug in the executable seems 
                                         //linked to the bad behaviour of copy constructor 
                                         //of FitnessComputer (compiler defined), please try 
                                         //and compile after uncommenting this line

  MPI_Request reqs [2];  //required for non-blocking calls
  MPI_Status stats [2];  //required for Waitall routine

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocess);
  MPI_Comm_rank (MPI_COMM_WORLD, & rank);


  double tstart = MPI_Wtime ();
  unsigned int ntimes = 1;



  for (unsigned int itime = 0; itime < ntimes; ++ itime) {
     //determine left and right neighbours
     next = (rank + 1 ) % nprocess;
     prev = (rank - 1 + nprocess) % nprocess;
   
   
     rnd.DoAll (10 + rank);
   
     size_t nCouples = std::pow (ncities - 2, 2) / 2.;
   
   
     Population pop (nCouples, ncities, rnd );
     std::vector <double> mean_lengths;
     std::vector <CodePerm> BestIndividuals ;
   
     //initializing output to file 
     std::string dirname ("output/dir10");
     std::string filename (name +"_fitness_" + std::to_string(rank) + ".dat");
     if ( !std::filesystem::exists (dirname)){
        std::system (("mkdir " + dirname).c_str());
     }
     std::ofstream output (dirname + '/' + filename);
     std::ofstream collectBest;
     if (rank == 0) {
        collectBest.open(dirname + '/' + name + "_common_fitness.dat");
     }

     if (rank == 0){
        std::ofstream Points (dirname + '/' + name + "_points.dat");
        Points << cities;
     }
   
     //we directly assign probabilities the values that proved best in previous exercice, WITH SelectInd USING EXPONENT = 3 
     std::vector <double> probabilities = (name == "circle" ? std::vector <double>{ 0.125, 0.125, 0.125, 0.625} : std::vector <double>{0.125, 0., 0.25, 0.875} );  
  
     unsigned int nMaxEquals = 30; //number of generations we shall go on even if the "bestForNow" value does not improve 
     unsigned int nequals = 0;
     double lastBest = ncities; //this way it should be always bigger than actual path lengths

     for (unsigned int iLoop = 0; iLoop < nSplitLoop; ++ iLoop){
        

        std::cerr << "loop n " + std::to_string(iLoop) + '\n';
   
        //reading best individuals and preparing array for migration
        for (int iMig = 0; iMig < nMigrants; ++iMig){
           std::vector <arma::uword> bufvec1 (pop.GetBest (environment, iMig * 2).GetGen () ), 
                                           bufvec2 (pop.GetBest (environment, iMig * 2 + 1).GetGen () ); 
           bufvec1.insert (bufvec1.end(), bufvec2.begin(), bufvec2.end());
           for (unsigned int i = 0, end = bufvec1.size(); i < end; ++i){
              sendbuf [iMig] [i] = static_cast <unsigned long int> (bufvec1 [i]); 
           } 
        }
    
        //calling MPI non-blocking routines
        MPI_Irecv (recbuf [0], nMigrants * buffer_size, MPI_UNSIGNED_LONG_LONG, prev, tag, MPI_COMM_WORLD, & reqs[0] );
   //     std::cout << rank << " waiting for buffer from " << prev << std::endl;
        MPI_Isend (sendbuf [0] , nMigrants * buffer_size, MPI_UNSIGNED_LONG_LONG, next, tag, MPI_COMM_WORLD, & reqs[1] );
   //     std::cout << rank  << " sending buffer to " << next << std::endl;
  
        //internal loop (no exchanges)
        for (unsigned int iGen = 0; iGen < nGenerations; ++iGen){
           pop.NextGen (environment, rnd, probabilities);
           BestIndividuals.emplace_back (pop.GetBest ( environment));
           std::sort (BestIndividuals.begin(), BestIndividuals.end(), environment);
   //writing to output file
           output << pop.MeanOfBestHalf (environment) << ' ' << BestIndividuals.front().GetLength() << ' ';
           BestIndividuals.erase ( BestIndividuals.begin() + 1, BestIndividuals.end() );
   
           for (const auto & i : BestIndividuals.front().GetFeno() ) {
              output << i << ' ';
           }
           output << std::endl;
        }
          
   //     std::sort (BestIndividuals.begin(), BestIndividuals.end(), environment); should be unnecessary after the internal loop
        std::vector <arma::uword> BestGenes (BestIndividuals.front().GetGen ());
      
   
        MPI_Waitall (2, reqs, stats);
   /*     for (int i = 0; i < buffer_size; ++i ){
           std::cerr << recbuf [i] << ' ';
        }
        std::cerr << std::endl;
   */
        
      
        for (int iMig = 0; iMig < nMigrants; ++iMig) {
        std::vector <arma::uword> gen1 (recbuf [iMig], recbuf [iMig]  + ncities - 2      ),
                                  gen2 (recbuf [iMig] + ncities - 2, recbuf [iMig] + (ncities - 2) * 2);

        // update popolation with migrated elements

            pop.AddCouple (CodePerm (gen1), CodePerm(gen2));
            pop.DiscardCouples (environment, 1);    
        }

        //preparing buffer for  communication to process 0
        int size = BestGenes.size();
        for (int i  = 0; i < size; ++i) {
           sendbuf [0] [i] = static_cast <unsigned long int> (BestGenes [i]);
        }
        
        //Gahter-ing the best genes from each process
        MPI_Barrier (MPI_COMM_WORLD);
//        int whooseTurn  = static_cast <int> (iLoop) % nprocess;          
        MPI_Gather ( sendbuf [0], buffer_size , MPI_UNSIGNED_LONG_LONG, & collectbuf [rank] [0] , buffer_size, MPI_UNSIGNED_LONG_LONG, 0 , MPI_COMM_WORLD);
        
        if (rank == 0) { // process 0 performs last steps:
           std::vector <CodePerm> processVec;
           for (int i = 0; i < nprocess; ++i) {
              processVec.emplace_back (std::vector <arma::uword> (collectbuf [i], collectbuf [i] + size));
           } 
           std::sort (processVec.begin(), processVec.end(), environment);
           auto bestForNow = processVec.begin();
           environment.ComputeLength (*bestForNow);
   
           //writing to output file
           double newBest = bestForNow -> GetLength();
           collectBest << newBest <<  ' ';
   //        for (auto const & code : processVec ) collectBest << code.GetLength() << ' ';
           for (auto const & el : bestForNow -> GetFeno() ) collectBest << el << ' ';
           collectBest << std::endl; 

           //counting loops without improvement
           nequals += nGenerations;
           if (lastBest > newBest ) nequals = 0;
           lastBest = newBest;
           if (nequals >= nMaxEquals) {
              endloop = true;
           }
        }

        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast ( &endloop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);

      
        if (endloop){ //exit condition (all processes)
           if (rank == 0) {
              std::string par_file_name = name + "_parameters_and_effectiveness.dat";
              std::ofstream savePars; 
              if ( !fs::exists(dirname + "/" + par_file_name)) {
                 savePars.open (dirname + "/" + par_file_name);
                 savePars << "nprocess nMigrants nGenerations nTotGenerations fitness"; 
              }
              else savePars.open (dirname + "/" + par_file_name, std::ios::app);
              savePars << std::endl << nprocess << ' ' << nMigrants << ' ' << nGenerations << ' ' << (iLoop + 1)*nGenerations << ' ' << lastBest;
           }
           break;
        }
     } 
   }

  double tend = MPI_Wtime();

  if (rank == 0){
     std::ofstream Time ("output/dir10/time10" + name + ".dat", std::ios::app);
     Time << nprocess << ' ' << tend - tstart << std::endl;
  }

  
  MPI_Finalize();

  delete sendbuf [0];
  delete recbuf [0]; 
  delete sendbuf;
  delete recbuf;

  return 0;
}
