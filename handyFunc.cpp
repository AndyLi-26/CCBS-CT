#ifndef HANDY_FUNC
#define HANDY_FUNC
#include "structs.h"
inline bool eq(double a, double b) {return abs(a-b)<=CN_EPSILON;}
inline bool lt(double a, double b) {return a-b < -CN_EPSILON;}
inline bool gt(double a, double b) {return a-b > CN_EPSILON;}
inline bool le(double a, double b) {return a-b <= CN_EPSILON;}
inline bool ge(double a, double b) {return a-b >= -CN_EPSILON;}

inline bool eq_raw(double a, double b) {return a==b;}
inline bool lt_raw(double a, double b) {return a<b;}
inline bool gt_raw(double a, double b) {return a>b;}
inline bool le_raw(double a, double b) {return a<=b;}
inline bool ge_raw(double a, double b) {return a>=b;}

inline void prt_double(double a)
{
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&a);
    for (int i =8 - 1; i >= 0; --i) {
        for (int j = 7; j >= 0; --j) {
            std::cout << ((bytes[i] >> j) & 1);
        }
        std::cout << ' '; // Separate bytes with a spacs
    }
}

inline double round_down (double f){
  return std::floor(f/CN_EPSILON)*CN_EPSILON;
}

inline Interval solveQuad(double a, double b, double c) {
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

inline unsigned int floorSingle(float number){
    unsigned int* ptr = (unsigned int*)&number;
    int s = *ptr >> 31;
    int e = *ptr & 0x7f800000;
    e >>= 23;
    e-=127;
    int m = *ptr & 0x007fffff;

    unsigned long temp=1;
    temp<<=23;
    unsigned int retval= m|temp;
    retval>>=(23-e-8);

}

inline unsigned long floorDouble(double number){
  unsigned long* ptr = (unsigned long*)&number;
  long e = *ptr & 0x7ff0000000000000;
  e >>= 52;
  e-=1023;
  long m = *ptr & 0x00fffffffffffff;

  unsigned long temp=1;
  temp<<=52;
  unsigned long retval= m|temp;
  retval>>=(52-e-8);
  return retval;
}

#endif
