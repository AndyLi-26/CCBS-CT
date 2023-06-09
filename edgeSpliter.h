#ifndef EDGESPLIT_H
#define EDGESPLIT_H
#include "const.h"
#include "structs.h"
#include "map.h"
#include "heuristic.h"
#include <math.h>
#include <string>
using namespace std;
class edgeSpliter
{
  public:
    Config config;
    double r;
    double d;
    ofstream output;
    bool writeNode;
    edgeSpliter(Config _config):config(_config), r(_config.agent_size), d(_config.agent_size*2) {
      writeNode= (config.F_debug_info!="" && config.debug>0);
    }
    ~edgeSpliter() {}
    void find_deltas(const Conflict &conf, Map_deltas &deltasR, Map_deltas &deltasL, Map &map, Heuristic &h_values);

  private:
    Vector2D case2(Vector2D P0,Vector2D v, Vector2D P2);
    void waiting(Move m1, Move m2,Map_deltas &deltas, Map &map, Heuristic &h_values, int a);
    void moving(Move m1, Move m2,Map_deltas &deltas, Map &map, Heuristic &h_values, int a);
    double round_down(double f);
    bool validNewNode(Vector2D node1, Vector2D node2, Vector2D New);
    double solveQuad(double a, double b, double c);
    Vector2D fitPoint(Vector2D P, int n1, int n2, Map &map);
};

#endif
