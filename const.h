#ifndef CONST_H
#define CONST_H

#define CN_USE_EDGE_SPLIT 1//1 -true, 0 - false, use edge spliting
#define CN_USE_CARDINAL  0 //1 - true, 0 - false
#define CN_HLH_TYPE      0 // 0 - no hlh, 1 - solve lpp by simplex, 2 - greedly take disjoint conflicts
#define CN_USE_DS        0 //1 - true, 0 - false
#define CN_TIMELIMIT     30 // in seconds
#define CN_AGENT_SIZE    sqrt(2.0)/4.0 //radius; only values in range (0; 0.5] are supported
#define CN_CONNECTEDNESS 2 // possible variants 2,3,4,5
#define CN_PRECISION     1e-4
#define CN_RESOLUTION    0.1
#define CN_FOCAL_WEIGHT  1.0 // experimental function, focal is supported only on the high-level
#define CN_OBSTL         1
#define CN_EPSILON       1e-4
//#define CN_MAP_RESOLUTION sqrt(2.0)/4.0
#define CN_INFINITY		 1e+8
#define CN_LOG           "_log"
#define CN_AGENT_NUM	 20

#define FAIL(message) std::cerr<<message << std::endl<<std::flush; exit(1);

#endif // CONST_H
