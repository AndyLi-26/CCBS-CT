#ifndef SIPP_H
#define SIPP_H
#include "structs.h"
#include "map.h"
//#include "config.h"
#include "heuristic.h"
#include "const.h"
#include <unordered_map>
#include <map>
#include <set>
class SIPP
{
public:

    Config config;
    SIPP(){}
    ~SIPP() {}
    Path find_path(Agent agent, const Map &map, std::list<Constraint> cons, Heuristic &h_values, bool p=false);

private:
    bool p=false;
    Agent agent;
    std::vector<Path> find_partial_path(std::vector<Node> starts, std::vector<Node> goals, const Map &map, Heuristic &h_values, double max_f = CN_INFINITY);
    Path add_part(Path result, Path part);
    void find_successors(Node curNode, const Map &map, std::list<Node> &succs, Heuristic &h_values, Node goal);
    void add_open(Node newNode);
    Node find_min();
    double dist(const Node &a, const Node &b);
    std::vector<Node> reconstruct_path(Node curNode);
    void make_constraints(std::list<Constraint> &cons, const Map &map);
    void clear();
    void add_collision_interval(int id, std::pair<double, double> interval);
    void add_move_constraint(Move move);
    std::vector<Node> get_endpoints(int node_id, double node_i, double node_j, double t1, double t2);
    double check_endpoint(Node start, Node goal,int exit_id);
    void prt_constraint(Constraint c);
    void prt_constraints(std::list<Constraint> constraints);
    void prt_cons();
    void prt_intervals();
    void prt_path(Path P);
    void prt_paths(std::vector<Path> paths);
    void prt_nodes(std::vector<Node> starts);

    std::unordered_map<int, Node> close;
    std::list<Node> open;
    std::unordered_map<int, std::pair<double, bool>> visited;
    std::map<std::pair<int, int>, std::vector<Move>> constraints;//stores sets of constraints associated with moves
    std::unordered_map<int, std::vector<std::pair<double, double>>> collision_intervals;//stores sets of collision intervals associated with cells
    std::vector<Move> landmarks;
    Path path;
};

#endif // SIPP_H
