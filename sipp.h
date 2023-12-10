#ifndef SIPP_H
#define SIPP_H
#include "structs.h"
#include "map.h"
#include "heuristic.h"
#include "handyFunc.cpp"
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <map>
#include <set>
class SIPP
{
public:

    Config config;
    SIPP()  {}
    ~SIPP() {}
    Path find_path(Agent agent, const Map &map, std::list<Constraint> cons, Heuristic &h_values,bool p=false);
    Path find_path_aux(Agent agent, const Map &map, std::list<Constraint> cons, Heuristic &h_values,bool p=false);

private:
    Agent agent;
    bool p;
    std::vector<Path> find_partial_path(std::vector<Node> starts, std::vector<Node> goals, const Map &map, Heuristic &h_values, double max_f = CN_INFINITY);
    Path add_part(Path result, Path part);
    void find_successors(Node curNode, const Map &map, std::list<Node> &succs, Heuristic &h_values, Node goal);
    void add_open(Node newNode);
    Node find_min();
    double dist(const Node &a, const Node &b);
    std::vector<Node> reconstruct_path(Node curNode);
    void make_constraints(std::list<Constraint> &cons,const Map &map);
    void clear();
    void add_collision_interval(int id, std::pair<double, double> interval);
    void add_move_constraint(Move move);
    std::vector<Node> get_endpoints(int node_id, double node_i, double node_j, double t1, double t2);
    pair<double,double> check_endpoint(Node start, Node goal,const Map &map);

    typedef std::tuple<int,int,bool> node_idx;
    //std::unordered_map<int, Node> close;
    boost::unordered_map<node_idx, Node> close;
    std::list<Node> open;
    //std::unordered_map<int, std::pair<double, bool>> visited;
    boost::unordered_map<node_idx, std::pair<double, bool>> visited;
    std::map<std::pair<int, int>, vector<Interval>> constraints;//stores sets of constraints associated with moves
    std::map<int, vector<Interval>> collision_intervals;//stores sets of collision intervals associated with cells
    std::vector<Move> landmarks;
    Path path;
    int call_find_pp=0;
    int node_open=0;
    int node_exp=0;
    void prt_landmarks();
    void prt_cons();
    void prt_node(Node n);
    void prt_nodes(vector<Node> nodes);
void prt_intervals();
void prt_constraints(list<Constraint> constraints);
void prt_info(){
cout<<"called find_partial_path: "<<call_find_pp<<endl;
cout<<"node opened: "<<node_open<<endl;
cout<<"node expanded: "<<node_exp<<endl;

}
};

#endif // SIPP_H
