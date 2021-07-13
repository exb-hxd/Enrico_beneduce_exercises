#ifndef __Genetic__
#define __Genetic__

#include "random.hpp"
#include "armadillo"

#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <utility>


class CodePerm {  // short for "Permutation Codifier"
   std::vector<arma::uword> genes;
   std::vector<arma::uword> fenotype;
   bool IsTranslated;
   double total_length;
   bool IsEvaluated;   

   friend class FitnessComputer;

   public: 
      CodePerm (size_t , Random & rnd ); 
      CodePerm (const std::vector <arma::uword> &);
      void Print (int) const; 
      void Print ( std::vector<int> ) const;
      void Translate ( std::vector<arma::uword> &);
      void TranslateBack ( std::vector <arma::uword>&);
      void Mutate1 (Random &, double);
      void Mutate2 (Random &, double);      
      void Mutate3 (size_t, arma::uword);
      void MutateFeno1 (size_t, size_t); //swaps
      void MutateFeno2 (size_t, size_t, int); //moves cities in circle
      void MutateFeno3 (size_t, size_t); //inverts cities in range

      void MutateFeno1Rand (Random & rnd, double p);       
      void MutateFeno2Rand (Random & rnd, double p);       
      void MutateFeno3Rand (Random & rnd, double p);

      void MutateFeno2size (Random &, double, double);
      void MutateFeno3size (Random &, double, double);



      inline void Translate () {
         Translate (fenotype);
         IsTranslated = true;
      }  

      inline void TranslateBack () {
         TranslateBack (genes);
         IsTranslated = true;
      }

      inline const std::vector <arma::uword> & GetFeno () const { 
         return fenotype;
      }

      inline const std::vector <arma::uword> & GetGen () const {
         return genes;
      }
      
      inline bool IsTr () const {
         return IsTranslated;
      }

      inline size_t GetNpoints () const { 
         return  (genes.size() + 2); 
      } 

      inline double GetLength () const {
         return total_length;
      }
      
      inline bool IsEval () const {
         return IsEvaluated;
      }
};


class FitnessComputer {

   std::shared_ptr <arma::mat>  distances_ptr;

   public:
      FitnessComputer (const arma::mat &i );
      bool operator () ( CodePerm & , CodePerm & );
      double ComputeLength (CodePerm & ) const ;
};

class Population {
   std::vector <CodePerm> individuals;
   std::vector <CodePerm> newPop;
   bool IsSorted;

   
   void Crossover (size_t ind1, size_t ind2, Random & , double ); 
   void Crossover (size_t ind1, size_t ind2, size_t cut, size_t length, Random &, double);

   public: 
      Population (size_t nCouples, size_t nPoints, Random & rnd); //not sorted
      void Evaluate (FitnessComputer const & env); 
 
      void Sort (FitnessComputer const & env) {
         std::sort (individuals.begin(), individuals.end(), env);
         IsSorted = true;
      }

      void Cut (FitnessComputer const & env ) {
         if ( !IsSorted ) Sort (env);
         individuals.erase ( individuals.begin() + individuals.size() / 2, individuals.end() );
      }

      size_t SelectInd ( Random & , size_t ExcludedValue = INT_MAX);
      void Reproduce (Random &, const FitnessComputer & , double ); //sets IsSorted to false
      void Reproduce (Random &, size_t cut, size_t length, const FitnessComputer & , double); // sets IsSorted to false
      void MutateRand (Random &, double , double, double) ; //sets IsSorted to false
      void MutateDet (Random & , std::vector <std::tuple <double, size_t, arma::uword>> &); //sets IsSorted to false
      void NextGen (const FitnessComputer & env, Random & rnd, std::vector <double> &);
      void AddCouple (CodePerm && , CodePerm &&) ;
      void DiscardCouples (const FitnessComputer &, size_t );

/*      template <typename ... Args>
      void MutateDet (Args && ... args ){ //sets IsSorted to false
         MutateDet (, std::vector <std::tuple <double, size_t, arma::uword>> (std::forward <Args> (args) ...));
      }  
*/ 
      void Print (size_t i, int iMember) const {
         individuals.at(i).Print (iMember); 
      }

      void Print (size_t i, std::vector<int> members) const {
         individuals.at(i).Print(members);
      }
      void PrintAll (std::vector<int>, int iPop) const;
       
      CodePerm & GetBest (const FitnessComputer & env) {
         return GetBest (env, 0);
      };

      CodePerm & GetBest (const FitnessComputer & , size_t );
      double MeanOfBestHalf (const FitnessComputer & env);

};

bool IsPermutation (const std::vector <arma::uword>& , const std::vector <arma::uword>& );

bool IsGenes (const std::vector <arma::uword> &);


#endif
