#include <iostream>
#include <cmath>
#include "random.hpp"
#include <string>
#include <fstream>

using namespace std;

int main (int argc, char * argv[]){

unsigned int cumulatives [] = {1, 2, 5, 10, 100, 0};
size_t M = 10000;
//here should be instruction to create the folder for output files! -to be fixed later- 
string rootname = "output/dir01-2/01-2-";


Random rnd;
rnd.DoAll();

//set parameters
rnd.SetParameter(MEAN, 0);
rnd.SetParameter(SIGMA, 1);
rnd.SetParameter(LAMBDA, 1);
rnd.SetParameter(GAMMA, 1);

for (unsigned int* ptr = cumulatives; *ptr != 0; ++ptr ){
   string nameWithNum = rootname + to_string(*ptr);
   string filenameGauss = nameWithNum + "Gauss"; 
   string filenameExp = nameWithNum + "Exp";
   string filenameLorentz = nameWithNum + "Lorentz";

   //creating output files
   ofstream outputG ( filenameGauss.c_str() );
   ofstream outputE ( filenameExp.c_str() );
   ofstream outputL ( filenameLorentz.c_str() );

   if (!(outputG.is_open() && outputE.is_open() && outputL.is_open()) ) {
      cerr << "unable to open output files\n";
      return 1; 
   }
   for (size_t i = 0; i < M; ++i){

      outputG << rnd.Cumulate(GAUSS, *ptr )/ *ptr << endl;
      outputE << rnd.Cumulate(EXP, *ptr)/ *ptr << endl;
      outputL << rnd.Cumulate(LORENTZ, * ptr)/ *ptr << endl; 

   }
}

return 0;
}
