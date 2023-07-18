#ifndef HANDY_FUNC
#define HANDY_FUNC
#include "structs.h"
inline bool eq(double a, double b) {return abs(a-b)<=CN_EPSILON;}
inline bool lt(double a, double b) {return a-b < -CN_EPSILON;}
inline bool gt(double a, double b) {return a-b > CN_EPSILON;}
inline bool le(double a, double b) {return a-b <= CN_EPSILON;}
inline bool ge(double a, double b) {return a-b >= -CN_EPSILON;}

inline double round_down (double f){
  return std::floor(f/CN_EPSILON)*CN_EPSILON;
}

inline pair<double,double> solveQuad(double a, double b, double c) {
  double delta(b*b-4*a*c);
  if (!(delta>-CN_EPSILON))
    cout<<"delta"<<delta<<endl;
  assert(delta>-CN_EPSILON);
  if (abs(delta)<CN_EPSILON) delta=0;
  double sqrtDelta(sqrt(delta));
  double t1( round_down( (-b-sqrtDelta)/(2*a) ));
  double t2( round_down( (-b+sqrtDelta)/(2*a) ));
  return make_pair(t1,t2);
}


#endif
