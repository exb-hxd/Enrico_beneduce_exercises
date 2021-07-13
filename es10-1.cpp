#include "tools.hpp"
#include "random.hpp"
#include "Genetic.h"
#include "GeneralAnnealer.h"
#include "armadillo"
#include "circle.h"

#include <iostream>
#include <filesystem>
#include <cstdlib>

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


int main (int argc, char **argv) {

   std::string name = "circle";
   if (argc > 1){
      std::string temp (argv[1]);
      if (temp == "square") {
         name = "square";
      }
   }

//creating output files  
   std::string dirname ("output/dir10");
   if (!std::filesystem::exists (dirname) ){
      std::cerr << "creating dir " << dirname;
      system ( ("mkdir " + dirname ).c_str() );
   }
   std::string filename ("annealing_progress_" + name + ".dat");
   std::ofstream writeOut (dirname + '/' + filename );

   std::ofstream writePoints (dirname + "/annealing_" + name + "_points.dat");

   Random rnd;
   rnd.DoAll(10);
   
   arma::uword ncities = 32;

   arma::mat cities = (name == "square" ? makeSquare (rnd, ncities) : makeCircle (rnd, ncities) ); 

   //saving cities coordinates
   writePoints << cities;

   FitnessComputer environment (cities);

   double Temp = 1.;

   CodePerm Sampler (ncities, rnd);
   double lastLength = environment.ComputeLength (Sampler);

   auto pr0Lambda = [&](CodePerm & c)->void {c.Mutate2 (rnd, 0.8);};
   auto pr1Lambda = [&](CodePerm & c)->void {c.MutateFeno1Rand (rnd, 1);  };
   auto pr2Lambda = [&](CodePerm & c)->void {c.MutateFeno2Rand (rnd, 1);  };
   auto pr3Lambda = [&](CodePerm & c)->void {c.MutateFeno3Rand (rnd, 1);  };
   auto accLambda = [&](CodePerm & c1)->double {return exp (- (environment.ComputeLength (c1))/ Temp); };
   auto calcLambda = [&](CodePerm & c1, CodePerm & c2)->double {return environment.ComputeLength (c1);};

   GeneralAnnealer <CodePerm> ann (ncities, rnd);
  
   arma::vec  TempVec = arma::logspace ( std::log10(ncities), std::log10 (0.01),  std::pow(ncities, 2));
   std:: cout << TempVec;
   unsigned int nsteps = 1000; 
   for (auto t : TempVec ){ 
      Temp = t;
      BlockMean meanLength (nsteps);
      for (unsigned int  i  = 0; i < nsteps; ++i ){
           lastLength = ann.Move (pr0Lambda, accLambda, calcLambda, rnd, lastLength);
           lastLength = ann.Move (pr1Lambda, accLambda, calcLambda, rnd, lastLength);
           lastLength = ann.Move (pr2Lambda, accLambda, calcLambda, rnd, lastLength);
           lastLength = ann.Move (pr3Lambda, accLambda, calcLambda, rnd, lastLength);

           meanLength.AddValue (lastLength);    
      }
      meanLength.AddMean(); //we only compute mean once for each temp, not block averages

      //writing to output file
      writeOut << Temp << ' ' << ann.GetAcceptRate () << ' ' << meanLength.GetMean() << ' ' << lastLength  << ' ';
      for (auto const & i : ann.GetPos().GetFeno() )   writeOut << i << ' ';
      writeOut << std::endl;

      //print progress to teminal
      ann.GetPos().Print({3});  
  //    double newnsteps = static_cast <unsigned int> (1000. /  (nsteps * ann.GetAcceptRate()) );
//      nsteps = (newnsteps > 100000 ? 100000 : (newnsteps < 100 ? 100 : newnsteps ));
   }
   return 0;

}
