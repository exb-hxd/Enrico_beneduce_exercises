#ifndef __Metropolis__
#define __Metropolis__

#include "random.hpp"

#include <algorithm>
#include <iostream>
#include <vector>


class Metropolis : private Random{
      size_t ncoordinates;
      size_t npoints;
      size_t point;
      std::vector<std::vector <double>> pos;
      std::vector<std::vector <double>> posnew;
      double probold;
      size_t accepted;
      size_t attempted;
      std::vector <std::vector <double>> saveStart;
      

      void Propose(Distribution );

      
/*      template <typename Functor >
      void ProposeT(Functor f) {
         posnew = f (posold)
         
      }
*/
      template <typename Functor>
      void Accept (Functor func){
         if (probold < 0) {
            probold = func(pos);
         }

         double probnew = func(posnew);
         double A = std::min(1., probnew/probold);
         if(Rannyu() <= A){
            pos[point] = posnew[point];
            probold = probnew;
            ++ accepted;
         }
      }

   public:   
      Metropolis(size_t , size_t , double , double , int nseed = 1);
//      Metropolis(size_t , size_t, double, double);
      ~Metropolis();
      Metropolis(const Metropolis& other);
      Metropolis& operator= (const Metropolis& other);
      friend  bool operator< (const Metropolis &lhs, const Metropolis & rhs);


      void Start(const std::vector<std::vector<double>> &);      

      inline void SetP(Parameters par, double value){
         SetParameter(par, value);
      }

      template <typename Functor>
      inline void Move(Distribution d, Functor f){
         Propose (d);
         Accept (f);
      }
      
      inline std::vector<std::vector<double>>* const GetPos() {
         return &pos;
      }
      
      inline double const GetAcceptRate () {
         return attempted == 0 ? 0 : double(accepted) / double(attempted);
      }
      void SetDelta (double);
      inline void SetSigma (double Value){
         SetParameter (SIGMA, Value);
      }
      inline double GetSigma (){
         return GetParameter(SIGMA);
      }
      inline double GetDelta () {
         return GetParameter (MAX);
      }   
      
      template <typename Functor>      
      void FixStep ( Distribution d, Functor f , size_t nsteps = 1000, size_t nMaxTrials = 10 , double tolerance = 0.1) {
         for ( double acceptance = GetAcceptRate() ; std::abs (acceptance - 0.5) > tolerance;){ 

            Start (saveStart);
            for (size_t istep = 0; istep < nsteps; ++istep) {
               Move (d, f);
            } 
            acceptance = GetAcceptRate();
            switch (d){
               case IN_RANGE:
                  SetDelta (std::abs ( GetParameter (MAX) ) * (acceptance / 0.5)  );
                  break;
               case GAUSS:
                  SetSigma (std::abs ( GetParameter (SIGMA) ) * (acceptance / 0.5) );
                  break;
               default:
                  throw std::invalid_argument ("move with distribution n " + std::to_string (d) + " not implemented\n" );
                  break;
            }
       
            if (nMaxTrials-- ==  0 ) { //nMaxTrials is changed at every cycle!
               std::cerr << "not  able to achieve required ratio: got " << acceptance << std::endl;
               return; 
            }
            
         }
      }
};


#endif 
