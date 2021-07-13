#include "random.hpp"
#include "Genetic.h"
#include "armadillo"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

//start CodePerm definitions
CodePerm :: CodePerm (size_t nPoints, Random & rnd ): genes(nPoints -2) , fenotype (nPoints), IsTranslated (false), total_length (-1), IsEvaluated (false) {
   for (size_t i = 0, iEnd = genes.size(); i < iEnd; ++i){
      genes[i] = static_cast <arma::uword> (rnd.Rannyu() * static_cast <double> (nPoints - i - 1));
   }
}        

CodePerm :: CodePerm (const std::vector <arma::uword> & init_gen): genes(std::move(init_gen)), fenotype (genes.size() + 2) ,  
                                               IsTranslated (false), total_length (-1), IsEvaluated (false)   {
   for (size_t i = 0, nGenes = genes.size(); i < nGenes; ++i) {
      if ( genes [i] > nGenes - i) {
         throw std::invalid_argument ("CodePerm : init_gen [" + std::to_string(i) + "] = " + 
                                       std::to_string (genes[i]) + ", max is " + std::to_string (nGenes - i));
      }
   }
}

void CodePerm :: Print ( int iMember ) const {
   switch (iMember) {
      case 0:
         std::for_each ( genes.begin(), genes.end(), [](arma::uword n){std::cout << ' ' << n;});
         std::cout << std::endl;
         break;
      case 1: 
         std::for_each (fenotype.begin(), fenotype.end(), [](arma::uword n){std:: cout << ' ' << n;});
         std::cout << std::endl;
         break;
      case 2: 
         std::cout << ( IsTranslated  ? "IsTranslated: true\n" : "IsTranslated: false\n") ;
         break;
      case 3:
         std::cout << ' ' << total_length << std::endl;
         break;
      case 4:
         std::cout << ( IsEvaluated  ? "IsEvaluated: true\n" : "IsEvaluated: false\n") ; 
         break;
      default:
         std::cerr << "CodePerm::Print() : invalid index\n";
   }
}

void CodePerm :: Print ( std::vector <int> mustPrint ) const {  //prints members of CodePerm corresponding to positions of '1' in input string
//   std::array <int, 5> mustPrint (list); // CodePerm only has 5 members
   for (int i = 0, n = mustPrint.size(); i < n; ++i ){
      Print ( mustPrint[i] );
   }
}

void CodePerm :: Translate ( std::vector <arma::uword> &result) {  //can also be used to write an equivalent of fenotype to output vector

   std::vector<arma::uword> ordered (genes.size() + 1);
   for( arma::uword i = 0, iEnd = ordered.size() ; i < iEnd; ++i){
      ordered[i] = i + 1;
   }

   result.clear (); 
   result.push_back (0); //first city is always 0-th one
         
   for (const auto & index : genes) {
      result.push_back (ordered [index]);
      ordered.erase (ordered.begin() + index); 
   }

   result.push_back (ordered [0]);
}

void CodePerm :: TranslateBack ( std::vector <arma::uword> & result) {
   result.clear();
   
   std::vector<arma::uword> ordered (fenotype.size ());
   for( arma::uword i = 0, iEnd = ordered.size() ; i < iEnd; ++i){
      ordered[i] = i + 1;
   }
   
   for (size_t i = 1, size = fenotype.size() - 1; i < size; ++i) {
      auto  it =  std::find (ordered.begin(), ordered.end(), fenotype [i]);
      result.push_back ( it - ordered.begin() );
      ordered.erase (it);
   }
}

void CodePerm :: Mutate1 (Random & rnd, double probability){
   if (probability > 1 || probability < 0 ) {
      throw std::invalid_argument ("CodePerm :: Mutate1: probability must lie between 0 and 1");
   } 
   if (rnd.Rannyu() > probability) return;
   
   size_t size = genes.size();
   size_t ind  = static_cast <size_t> (std::floor (rnd.Rannyu() * static_cast <double> (size)));

   genes[ind] = static_cast <arma::uword> ( (size - ind) - genes.at(ind) );

   
   if (! IsGenes (genes) ) {
      throw std::domain_error ("CodePerm: value of genes went out of range");
   }
   
   IsTranslated = false;
   IsEvaluated = false;

}

void CodePerm :: Mutate2 (Random & rnd, double probability) {
   if (probability > 1 || probability < 0 ) {
      throw std::invalid_argument ("CodePerm :: Mutate1: probability must lie between 0 and 1");
   } 
   if (rnd.Rannyu() > probability) return;

   size_t size = genes.size();
   size_t start = static_cast <size_t> (std::floor (rnd.Rannyu() * static_cast <double> ( size - 1)));
   
   for (size_t i = start, end = start + static_cast <size_t> (std::floor (rnd.Rannyu() * static_cast <double> ( size - start)) + 1.); i < end; ++i){


      genes [i] = static_cast <arma::uword> ( (static_cast <int> (genes[i]) + static_cast <int> (std::floor (rnd.Rannyu() * static_cast <double> (i)))) % 
                                                               static_cast <int>(size - i + 1)   ) ;
   }

   if (! IsGenes (genes) ) {
      throw std::domain_error ("CodePerm: value of genes went out of range");
   }

   
   IsTranslated = false;
   IsEvaluated = false;
}

void CodePerm :: Mutate3 (size_t index, arma::uword increment) {
   size_t size = genes.size();
   try{
      genes[index] = static_cast <arma::uword> (static_cast <int> (genes.at (index) + increment) % static_cast<int>(size - index + 1) );
   }
   catch (std::out_of_range & e) {
      std::cout << "CodePerm: " << e.what();
      throw e;
   }
   if (! IsGenes (genes) ) {
      throw std::domain_error ("CodePerm: value of genes went out of range");
   }
   
   IsTranslated = false;
   IsEvaluated = false;
}

void CodePerm :: MutateFeno1 (size_t pos1, size_t pos2) {

   size_t nMovablePoints = genes.size() + 1;
   if (pos1 >= nMovablePoints ) throw std::invalid_argument ("MutateFeno1: pos1 = " + std::to_string (pos1) + ", max value is " + std::to_string (nMovablePoints) );
   if (pos2 >= nMovablePoints ) throw std::invalid_argument ("MutateFeno1: pos2 = " + std::to_string (pos2) + ", max value is " + std::to_string (nMovablePoints) );

   if (!IsTranslated) Translate();

   std::swap (fenotype [pos1 + 1], fenotype [pos2 + 1] );
//   IsTranslated = false;
   TranslateBack ();
   IsEvaluated = false;
}

void CodePerm :: MutateFeno2 (size_t start, size_t length, int nmove ) {
   
   size_t nMovablePoints = genes.size() + 1;
   if ( length < 2 ) throw std::invalid_argument ("MutateFeno2: length = " + std::to_string (length) + 
                                                       " < 2: mutation has no effect");
   if (start + length > nMovablePoints ) throw std::invalid_argument ("MutateFeno2: start + length = " + std::to_string (start + length) + 
                                                                ", max value is " + std::to_string (nMovablePoints ) );
   if (!IsTranslated) Translate();

   size_t cut = static_cast <size_t> (  start + ( (nmove % static_cast <int> (length) + length ) % static_cast <int> (length) ) );
   std::vector <arma::uword> tmp (fenotype.begin() + cut + 1, fenotype.begin()  + start + length + 1 );
   tmp.insert (tmp.end(), fenotype.begin() + start + 1, fenotype.begin() + cut +1);

   std::swap_ranges (tmp.begin(), tmp.end(), fenotype.begin() + start + 1);
//   IsTranslated = false;
   TranslateBack();
   IsEvaluated = false;
}

void CodePerm :: MutateFeno3 (size_t start, size_t length) {

   size_t nMovablePoints = genes.size() + 1;
   if (length < 2 ) throw std::invalid_argument ("MutateFeno3: length = " + std::to_string (length) + 
                                                       " < 2: mutation has no effect ");
   if (start + length > nMovablePoints ) throw std::invalid_argument ("MutateFeno3: start + length = " + std::to_string (start + length) + 
                                                                ", max value is " + std::to_string (nMovablePoints) );

   if (!IsTranslated) Translate();

   std::reverse (fenotype.begin() + start + 1, fenotype.begin() + start + 1 + length);

   TranslateBack();
   IsEvaluated = false;
}


void CodePerm :: MutateFeno1Rand (Random & rnd, double p) {
   if (p < 0 || p > 1) throw std::invalid_argument ("CodePerm :: MutateFeno1Rand: probability should lie betweeen 0 and 1" );
   if (rnd.Rannyu () > p ) return;    
    
   size_t ind1 = static_cast <size_t> (std::floor (genes.size()  * rnd.Rannyu() ));
   size_t ind2 = static_cast <size_t> (std::floor (genes.size()  * rnd.Rannyu() ));
   if (ind2 >= ind1) ++ ind2;
   MutateFeno1 (ind1, ind2);
}

void CodePerm :: MutateFeno2Rand (Random & rnd, double p) {
   if (p < 0 || p > 1) throw std::invalid_argument ("CodePerm :: MutateFeno2Rand: probability should lie betweeen 0 and 1" );
   if (rnd.Rannyu () > p ) return;    
        
   size_t start = static_cast <size_t> (std::floor (genes.size() * rnd.Rannyu()));
   size_t length = static_cast <size_t> (std::floor ((genes.size() - start) * rnd.Rannyu() + 2));
   int move = static_cast <int> (std::floor ( (length - 1) * rnd.Rannyu() + 1 ));
   MutateFeno2 (start, length, move);
}

void CodePerm :: MutateFeno3Rand (Random & rnd, double p) {
   if (p < 0 || p > 1) throw std::invalid_argument ("CodePerm :: MutateFeno3Rand: probability should lie betweeen 0 and 1" );
   if (rnd.Rannyu () > p ) return;    
        
   size_t start = static_cast <size_t> (std::floor (genes.size()  * rnd.Rannyu() ));
   int  length = static_cast <size_t> (std::floor( (genes.size() - start ) * rnd.Rannyu()) + 2);
    MutateFeno3 (start, length); 
}
  
void CodePerm :: MutateFeno2size (Random & rnd, double p , double scale_factor){

   if (p < 0 || p > 1) throw std::invalid_argument ("CodePerm :: MutateFeno2size: probability should lie betweeen 0 and 1" );
   if (scale_factor < 0 || scale_factor > 1) throw std::invalid_argument ("CodePerm :: MutateFeno1Rand: scale_factor should lie betweeen 0 and 1" );  
   if (rnd.Rannyu () > p ) return;    
     
   size_t start = static_cast <size_t> (std::floor (genes.size() * rnd.Rannyu() * scale_factor));
   size_t length = static_cast <size_t> (std::floor ((genes.size() - start) * rnd.Rannyu() * scale_factor + 2));
   int move = static_cast <int> (std::floor ( (length - 1) * rnd.Rannyu() + 1 ));
   MutateFeno2 (start, length, move);
}


void CodePerm :: MutateFeno3size (Random & rnd, double p, double scale_factor) {
   if (p < 0 || p > 1) throw std::invalid_argument ("CodePerm :: MutateFeno3size: probability should lie betweeen 0 and 1" );
   if (scale_factor < 0 || scale_factor > 1) throw std::invalid_argument ("CodePerm :: MutateFeno1Rand: scale_factor should lie betweeen 0 and 1" );  
   if (rnd.Rannyu () > p ) return;    
        
   size_t start = static_cast <size_t> (std::floor (genes.size()  * rnd.Rannyu() * scale_factor ));
   int  length = static_cast <size_t> (std::floor( (genes.size() - start ) * rnd.Rannyu() * scale_factor ) + 2);
    MutateFeno3 (start, length); 
}
/*
void CodePerm :: MutateFeno2RandSmall (Random rnd, double p, double ) {
   

} */   
// end CodePerm definitions

//start FitnessComputer definitions

FitnessComputer :: FitnessComputer (const arma::mat & points) {
   arma::uword nPoints = arma::size (points) [1] ;
   arma::uword nDim = arma::size (points) [0];
   std::cout << "FitnessComputer: constructor argument interpreted as " << nPoints
             <<  " points in " << nDim << " dimensions\n";

   distances_ptr = std::make_shared <arma::mat> (nPoints, nPoints);

   for (arma::uword i = 0; i < nPoints; ++i){
      for (arma::uword j = 0; j < nPoints; ++j){
        (*distances_ptr) (j,i) = arma::norm ( points.col(j) - points.col(i) , 2 );  // since we are only interested in the cities order, choosing 
      }									    // a different exponent for the norm computation makes no difference. 
   } 									    // "1" looks like being slightly computationally smarter

   std::cout << "FitnessComputer: distance matrix\n" <<  *distances_ptr;          
}

bool FitnessComputer :: operator ()  ( CodePerm & former,  CodePerm & latter)  {
   size_t nPoints_for = former.GetNpoints();
   if (nPoints_for != latter.GetNpoints()) {
      throw std::invalid_argument( "FitnessComputer :: operator () (CodePerm &, CodePerm &): unmatching sizes");
   }

   return ComputeLength (former) < ComputeLength (latter);   //best elements first
}

double FitnessComputer :: ComputeLength (CodePerm & individual) const { 
   if ( !individual.IsEvaluated ){
      size_t nPoints = individual.GetNpoints();
      if (nPoints > arma::size ( *distances_ptr )[0] ){
         throw std::length_error (  "FitnessComputer :: operator() (CodePerm & ): too many points (" + std::to_string (nPoints) 
                                  + " expected " + std::to_string (arma::size ( *distances_ptr ) [0])                            );
      }
      if( individual.IsTranslated == false) {
         individual.Translate();
      }
      
      arma::Col <arma::uword >  indices (individual.fenotype );
 //     std::cout << arma::trans (indices);
 //     std::cout << std::endl;

 //     std::cout << (*distances_ptr) ( indices, arma::join_cols (indices.subvec ( 1, nPoints - 1 ), indices.subvec (0, 0 ) )); 

      individual.total_length = arma::trace ( (*distances_ptr) ( indices, arma::join_cols (indices.subvec ( 1, nPoints - 1 ), indices.subvec (0, 1 ) )));
      individual.IsEvaluated = true;
   }
   return individual.total_length;
}
 
//end Fitness Computer definitions


//start Population definitions

Population :: Population (size_t nCouples, size_t nPoints, Random & rnd) {
   for (size_t iCouple = 0; iCouple < nCouples; ++iCouple ) {
      for (size_t iParent = 0, nParents = 2; iParent < nParents; ++iParent){ //it is assumed that nParents individuals a time
         individuals.emplace_back(nPoints, rnd);                              //mix up into a new individual
      }
   }
   IsSorted = false;
}


void Population :: Evaluate (FitnessComputer const & env) {
   std::for_each (individuals.begin(), individuals.end(), [=] (CodePerm & el) -> void {env.ComputeLength (el);} );
}

size_t Population :: SelectInd ( Random & rnd, size_t ExcludedInd ) {
   size_t ind = static_cast <size_t> (std::floor ((individuals.size() - 1) * std::pow (rnd.Rannyu(), 3)));
   if (ind >= ExcludedInd) ++ ind;
   return ind;
}


void Population :: Crossover (size_t ind1, size_t ind2, Random & rnd, double p ){  //does not check whether individuals[ind1].size() == individuals[ind2].size() !!
   //does not check wheter p is a true probability
   if (rnd.Rannyu() > p) {
      newPop.push_back ( individuals.at(ind1));
      newPop.push_back ( individuals.at(ind2));
   }
   else { 
      std::vector <arma::uword> cpy1 (individuals.at(ind1).GetGen()), cpy2 (individuals.at(ind2).GetGen());

      size_t startCross = static_cast <size_t> (std::floor (rnd.Rannyu() * static_cast <double> (cpy1.size() - 1)));
      size_t crossSize = static_cast <size_t> (std::ceil (rnd.Rannyu() * static_cast <double> (cpy1.size() - startCross - 1)));
      
      std::swap_ranges (cpy1.begin() + startCross, cpy1.begin() + startCross + crossSize, cpy2.begin() + startCross);        
      newPop.emplace_back ( cpy1 );
      newPop.emplace_back ( cpy2 );
   }
}

void Population :: Crossover (size_t ind1, size_t ind2, size_t cut, size_t length, Random & rnd, double p) { 
   //does not check wheter p is a true probability
   if (rnd.Rannyu() > p) {
      newPop.push_back ( individuals.at(ind1));
      newPop.push_back ( individuals.at(ind2));
   }
   else { 
      std::vector <arma::uword> cpy1 (individuals.at(ind1).GetGen()), cpy2 (individuals.at(ind2).GetGen());

      if (cut + length > cpy1.size()) {
         size_t checked_length = cpy1.size() - cut;
         std::cerr << "CrossoverDet : argument \"length\" = " << length << " is too big, reducing it to " << checked_length << std::endl;
         length = checked_length;
      }
      
      std::swap_ranges (cpy1.begin() + cut, cpy1.begin() + cut + length, cpy2.begin() + cut);        
      newPop.emplace_back ( cpy1 );
      newPop.emplace_back ( cpy2 );
   }
}

void Population :: Reproduce (Random & rnd, const FitnessComputer & env, double p) {
  if ( p < 0 || p > 1) throw std::invalid_argument ("Reproduce: given p is not a probability");
  if (!IsSorted) Sort (env);
  for (size_t i  = 0, size = individuals.size(); i < size; i += 2){
     size_t ind1 = SelectInd (rnd);
     size_t ind2 = SelectInd (rnd, ind1);
     Crossover (ind1, ind2, rnd, p);
  }

//  individuals.insert( individuals.end(), std::make_move_iterator (newPop.begin()) , std::make_move_iterator (newPop.end()) );
  individuals = std::move (newPop);  //less conservative
  newPop.clear();  //probably with no effect 
  IsSorted = false;
}

void Population :: Reproduce (Random & rnd, size_t cut, size_t length, const FitnessComputer & env, double p) {
   if ( p < 0 || p > 1) throw std::invalid_argument ("Reproduce: given p is not a probability");
   if (!IsSorted) Sort (env);
   for (size_t i = 0, size = individuals.size(); i < size; i += 2) { 
     size_t ind1 = SelectInd (rnd);
     size_t ind2 = SelectInd (rnd, ind1);
     Crossover (ind1, ind2, cut, length, rnd, p);
   }
//   individuals.insert( individuals.end(), std::make_move_iterator (newPop.begin()) , std::make_move_iterator (newPop.end()) ); //conservative
   individuals = std::move (newPop);  //less conservative
   newPop.clear();  //probably with no effect 
   IsSorted = false;
}

void Population :: MutateRand ( Random & rnd, double p1, double p2, double p3) {
   std::for_each ( individuals.begin(), individuals.end(), [&]( CodePerm & el )->void {el.MutateFeno1Rand(rnd, p1 ); 
                                                                                       el.MutateFeno2Rand(rnd, p2); 
                                                                                       el.MutateFeno3Rand(rnd, p3); 
                                                                                 });									     
   IsSorted = false;
}

void Population :: MutateDet (Random & rnd, std::vector <std::tuple <double, size_t, arma::uword> > & vec) {

   for (auto const & [p, ind, inc] : vec ) {
//      std::cout << "tuple: " << p << ' ' << ind << ' ' << inc << std::endl;
      if ( p < 0 || p > 1 ) throw std::invalid_argument ("Population :: MutateDet: " + std::to_string (p) + " not a probability"); 
   }
   
   std::sort (vec.begin(), vec.end()); //sorts for ascending values of p

//   for (auto const & [p, ind, inc] : vec )  
//      std::cout << "tuple: " << p << ' ' << ind << ' ' << inc << std::endl;
   
   std::for_each ( individuals.begin(), individuals.end(), [&]( CodePerm & el )->void {
      auto it = std::lower_bound (vec.begin(), vec.end(), std::tuple <double, size_t, arma::uword> (rnd.Rannyu(), 0, 0));
      if ( it == vec.end() ) { return;}
//      std::cout << ' ' << std::get <1> (*it) << ' ' ;
//      std::cout << std::get<0> (*it) << ' ' << std::get<1> (*it) << ' ' << std::get<2> (*it) << std::endl;
//      el.Print ({0}); 
      el.Mutate3 (std::get<size_t>(*it), std::get <arma::uword> (*it) );
//      el.Print ({0});
   });
   
   IsSorted = false;                                                                                        
}
/*
void Population :: MutateFeno (Random & rnd, )

//void Population :: MutateDet (Random & rnd, std::vector <double> p_vec, std::vector <size_t> index_vec, std::vector <arma::uword> increment_vec) 
*/
void Population :: NextGen (const FitnessComputer & env, Random & rnd, std::vector <double> & vec) {
//   size_t cut =  /*individuals [0].GetGen().size() / 2; */static_cast <size_t> (std::floor (rnd.Rannyu() * static_cast <double> (individuals[0].GetGen().size() - 1)));
//   size_t length = individuals[0].GetGen().size() - cut - 1;   //static_cast <size_t> (std::ceil (rnd.Rannyu() * static_cast <double> (individuals[0].GetGen().size() - cut - 1)));

      Reproduce (rnd, env, vec.at(3));
//   for (int itime = 0; itime < 1; ++itime){
      MutateRand (rnd, vec.at (0), vec.at(1), vec .at(2));
//   }

//   Cut (env);
}

void Population :: AddCouple (CodePerm && one, CodePerm && two) {
   individuals.emplace_back (std::move (one));
   individuals.emplace_back (std::move (two));
   IsSorted = false;
}

void Population :: DiscardCouples (const FitnessComputer & env, size_t nWorse) {
   if ( !IsSorted) Sort (env);
   individuals.erase (individuals.end() - nWorse * 2, individuals.end());
}


void Population :: PrintAll (std::vector <int> members, int iPop ) const {

   switch (iPop){
      case 0:
         for (const auto & el : individuals) {
            el.Print (members);
         }
         break;
      case 1: 
        for (const auto & el : individuals) {
           el.Print (members);
        }
        break;
      default:
         std::cerr << "Population :: PrintAll (iMember, iPop): invalid argument. iPop must be less than 2";
         break;
   }

}

CodePerm & Population :: GetBest (const FitnessComputer & env, size_t nBest) {
   if (!IsSorted) Sort(env);
   return individuals.at (nBest);
}

double Population :: MeanOfBestHalf (const FitnessComputer & env) {
  if (!IsSorted) Sort (env);
  double result = 0.;

  std::for_each (individuals.begin(), individuals.begin() + individuals.size() / 2, [&](const CodePerm & code )->void { result += code.GetLength(); });
  return result / static_cast <double> (individuals.size() / 2);
} 

 
// end Population definitions

//check functions
bool IsPermutation ( const std::vector <arma::uword>& one, const std::vector <arma::uword>& two, bool IsOrdered = false){
   if (one.size() != two.size()) {
      throw std::length_error( "IsPermutation(): wrong size" ) ; 
   }   

   auto onecpy(one), twocpy (two);

   if ( !IsOrdered ) std::sort ( onecpy.begin(), onecpy.end());
   std::sort (twocpy.begin(), twocpy.end());
   
   return onecpy == twocpy;
}


bool IsGenes (const std::vector <arma::uword> & sequence) {
   
   for (size_t i = 0, size = sequence.size(); i < size; ++i){
      if ( sequence [i] > size - i ) return false;
   } 
   return true;

}
//end check functions
