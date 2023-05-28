#include "handyFunc.h"
bool eq(double a, double b) {return abs(a-b)<=CN_EPSILON;}
bool lt(double a, double b) {return a-b < -CN_EPSILON;}
bool gt(double a, double b) {return a-b > CN_EPSILON;}
bool le(double a, double b) {return a-b <= CN_EPSILON;}
bool ge(double a, double b) {return a-b >= -CN_EPSILON;}
