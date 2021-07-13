#ifndef __GeneralAnnealer__
#define __GeneralAnnealer__

#include "random.hpp"

#include <algorithm>
#include <iostream>
#include <vector>



template <class Sampler> 
class GeneralAnnealer {

      Sampler posold;
      Sampler posnew;
      double pold;

      unsigned int accepted;
      unsigned int attempted;

      template <typename Functor>  // Functor should take an argument of type 
      void Propose ( Functor f) {  // Sampler by reference and modify its state in some way
         f (posnew);
         ++attempted;
      }
       

      template <typename Functor>   
      void Accept (Functor prob, Random & rnd){ //prob should be T(y|x)p(x) / T(x/y)p(y) in order to do real Metropolis sampling
         if (pold < 0 ) {
            pold = prob (posold);
         }         
         double pnew = prob (posnew); 
         double A = std::min(1., pnew / pold);
         if(rnd.Rannyu() <= A){
            ++ accepted;
            posold = posnew;
            pold = pnew;
            
         }
         else {
            posnew = posold;
         }
      }

      template <typename Functor, typename Calculator >   //Calculator may be used to compute quantities of interest
      auto Accept (Functor prob, Calculator calc, Random & rnd, decltype (calc(posnew, posold)) defaultResult ) -> decltype (calc (posnew, posold) ) {
         if (pold < 0 ) {
            pold = prob (posold);
         }
         
         using R = decltype ( calc (posnew, posold)); 
         R result; 
         double pnew = prob (posnew);

//         std::cout << "new " << pnew << " old " << pold;
         double A = std::min(1., pnew / pold);
         
/*         std::cout << "n: ";
         for (auto const & i : posnew) std::cout << i << ' ' ;
         std:: cout << "o: " ;
         for (auto const & i : posnew) std::cout << i << ' ' ;
         std::cout << std :: endl;
*/ 
         if(rnd.Rannyu() <= A){
            ++ accepted;
            result = calc (posnew, posold);
            posold = posnew;
            pold = pnew;
//            std:: cout << " acc\n" ; 
         }
         else {
//         std::cout << std::endl;
            posnew = posold;
            return defaultResult;
         }
         return result;
      }
 
   public:   
      GeneralAnnealer (const Sampler & init): posold (init), posnew (init), pold (-1), accepted (0), attempted (0) {}

      template<typename... Args>
      GeneralAnnealer (Args&&... args) : posold (std::forward<Args>(args)...), posnew (posold), accepted (0), attempted (0) {}  

      void Start(const Sampler & init) {
         posold = init;
         posnew = init;
   
         accepted = 0;
         attempted = 0;
      };      


      template <typename FunctorP, typename FunctorA>
      inline void Move (FunctorP p, FunctorA a, Random & rnd){
         Propose (p);
         Accept (a, rnd);
      }

      template <typename FunctorP, typename FunctorA, typename Calculator>
      inline auto Move (FunctorP p, FunctorA a, Calculator c, Random & rnd, decltype (c(posnew, posold) ) defres) -> decltype ( c (posnew, posold) ){
         Propose (p);
         return Accept (a, c, rnd, defres);
      }
      
      inline const Sampler & GetPos() const{
         return posold;
      }
      
      inline double const GetAcceptRate () {
         return double(accepted) / double(attempted);
      }

};


#endif 
