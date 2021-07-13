#include <iostream>
#include <fstream>
#include "random.hpp"
#include <cmath>
#include <cstdlib>
#include "tools.hpp"
#include <string>

/*
class BlockMean {
      double val, mean, sqMean;
      size_t n, nPerBlock;

   public: 
      BlockMean( size_t m ):val(0), mean(0), sqMean(0), n(0), nPerBlock(m) {}
      ~BlockMean(){}

     void AddValue (double Value) {
        val += Value;
     }

     void AddMean() {
	val *= 1./nPerBlock;
        if (n == 0) {
	   mean += val;
	   sqMean += val*val;
	   ++n;
	}
	else{
	   mean = mean * (n/double(n+1)) + val * (1./(n+1));
	   sqMean = sqMean * (n/double(n+1)) + val*val * (1./(n+1));
	   ++n; 
	}
	val = 0;
     }  

     inline double GetMean(){return mean; }
     inline double GetSqMean() {return sqMean;}
     inline size_t GetN() {return n;}

};


double statUncertainty( double mean, double sqMean, size_t N){
	if (N <= 1) return 0.;

	return sqrt((sqMean-pow(mean,2))/(N-1));
}
*/
using namespace std;

int main (int argc, char ** argv) {

   //definining some constants:	
   const double InitialAssetPrice = 100;
   const double DeliveryTime = 1;
   const double StrikePrice = 100;
   const double RiskFreeRate = 0.1;
   const double Volatility = 0.25;

   //defininig functors to calculate option price
   auto CallPrice = [T = DeliveryTime, r = RiskFreeRate, K = StrikePrice] 
      (double S) -> double 
   {
      return exp(-r*T) * max(0., S-K);
   };

   auto PutPrice = [T = DeliveryTime, r = RiskFreeRate, K = StrikePrice] 
      (double S) -> double 
   {
      return exp(-r*T) * max(0., K-S);
   };

   auto GBMasset = [S0 = InitialAssetPrice, T = DeliveryTime, K = StrikePrice, r = RiskFreeRate, sigma = Volatility] 
      (double Z, double t = -1, double S = -1) -> double 
   {  
      if(S == -1 ){
         S = S0;
      }
      if(t == -1 ){
         t = T;
      }
      return S * exp((r - sigma*sigma/2) * t + sigma * Z * sqrt(t));
   };

   //setting sample size:
   const size_t Nthrows = 10000;
   const size_t Nblocks = 100;
   const size_t NthrowsPerBlock = Nthrows / Nblocks;
   
   //number of steps in discretization of GBM:
   const size_t Nsteps = 100;
   const double  dt = DeliveryTime / Nsteps; 

   //initializing variables for block avarages:
   BlockMean CallPrice_1leap(NthrowsPerBlock);
   BlockMean PutPrice_1leap(NthrowsPerBlock);
   BlockMean CallPrice_nLeaps(NthrowsPerBlock);
   BlockMean PutPrice_nLeaps(NthrowsPerBlock);

   //initializing random number generator:
   Random rnd;
   rnd.DoAll();

   rnd.SetParameter(MEAN, 0);
   rnd.SetParameter(SIGMA, 1);

   //creating output files; we first write information abuot the parameters used
   system ("mkdir output/dir03-1");
   ofstream output1 ("output/dir03-1/price1step");
   ofstream outputN ("output/dir03-1/priceNsteps");

   string info ("");
   string offset (" i i i");
   info += "InitialAssetPrice: "; info += to_string (InitialAssetPrice);
   info += "\nDeliveryTime: "; info += to_string (DeliveryTime); 
   info += "\nStrikePrice: "; info += to_string (StrikePrice);
   info += "\nRiskFreeRate: "; info += to_string (RiskFreeRate);
   info += "\nVolatility: "; info += to_string (Volatility); 
   info += "\nFirstCol: Call"; 
   info += "\nSecondCol: Put";


   outputN << info << "\nNsteps: " << to_string (Nsteps);// << offset;

   output1 << info << "\nEndPar\n ";
   outputN << "\nEndPar\n ";

   //off we go!
   for (size_t i = 0; i < Nblocks; ++i) {
      for (size_t j = 0; j < NthrowsPerBlock; ++j){
	 double val = rnd.Gauss();
         CallPrice_1leap.AddValue ( CallPrice( GBMasset( val ) ) );
	 PutPrice_1leap.AddValue ( PutPrice ( GBMasset( val ) ) );

	 //discretized path:
	 double S = InitialAssetPrice;
	 for (size_t k = 0; k < Nsteps; ++k){
	    S = GBMasset ( rnd.Gauss(), dt, S );
	 }
	 CallPrice_nLeaps.AddValue( CallPrice ( S ));
	 PutPrice_nLeaps.AddValue( PutPrice ( S ));
      }

      //updating block averages;
      CallPrice_1leap.AddMean();
      PutPrice_1leap.AddMean();
      CallPrice_nLeaps.AddMean();
      PutPrice_nLeaps.AddMean();

      output1 << i << ' ' << CallPrice_1leap.GetMean() << ' ' << statUncertainty(CallPrice_1leap.GetMean(), CallPrice_1leap.GetSqMean(), i) <<
	              ' ' << PutPrice_1leap.GetMean() << ' ' << statUncertainty(PutPrice_1leap.GetMean(), PutPrice_1leap.GetSqMean(), i) << endl;

      outputN << i << ' ' << CallPrice_nLeaps.GetMean() << ' ' << statUncertainty(CallPrice_nLeaps.GetMean(), CallPrice_nLeaps.GetSqMean(), i) <<
	              ' ' << PutPrice_nLeaps.GetMean() << ' ' << statUncertainty(PutPrice_nLeaps.GetMean(), PutPrice_nLeaps.GetSqMean(), i) << endl;

   
   }



return 0;
}
