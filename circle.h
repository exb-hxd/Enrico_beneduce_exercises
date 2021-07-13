#include "random.hpp"
#include "armadillo"
#include <cmath>


class PointOnCircle {
   double radius;
   
   public:
      PointOnCircle (double rad) : radius (rad) {}

      arma::vec randCoord (Random & rnd) {
         double randAngle = rnd.Rannyu() * 2 * M_PI;
         return arma::vec {std::cos (randAngle), std::sin (randAngle)};   
      }

};
