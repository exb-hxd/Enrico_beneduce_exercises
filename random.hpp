/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__


enum Distribution {IN_RANGE, GAUSS, EXP, LORENTZ};
enum Parameters {MIN, MAX, MEAN, SIGMA, GAMMA, LAMBDA};

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

  double min, max, mean, sigma, Gamma, lambda;
protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  void SetParameter (Parameters, double);
  double GetParameter  (Parameters) const;
  double Rannyu(void);
  double RannyuInRange(void);
  double RannyuAngle(void);
  double Exp (void);
  double Lorentz(void);
  double Cumulate(Distribution, unsigned int); 
  double Gauss(void);
  void DoAll(int);
  void DoAll(void);
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
