#ifndef SIPP_H
#define SIPP_H
#include "structs.h"
#include "map.h"
#include "heuristic.h"
#include "handyFunc.cpp"
#include <unordered_map>
//#include <boost/unordered_map.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/unordered_set.hpp>
#include <map>
#include <set>
#include "node_pool.h"
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
    std::vector<Path> search(std::vector<Node*> starts, std::vector<Node*> goals, const Map &map, Heuristic &h_values, double max_f = CN_INFINITY);
    Path add_part(Path result, Path part);
    void find_successors(Node* curNode, const Map &map, std::list<Node*> &succs, Heuristic &h_values, Node goal);
    //void add_open(Node newNode);
    Node* find_min();
    //double dist(const Node &a, const Node &b);
    std::vector<Node> reconstruct_path(Node* curNode);
    void make_constraints(std::list<Constraint> &cons,const Map &map);
    void clear();
    void add_collision_interval(int id, std::pair<double, double> interval);
    void add_move_constraint(Move move);
    vector<Node*> get_endpoints(const int node_id, const double node_i, const double node_j, const double t1, const double t2);
    vector<Node*> endpoints_from_prev(vector<Path> results);
    pair<double,double> check_endpoint(Node* start, Node* goal,const Map &map);

	size_t nodePoolSize= (500000);
	node_pool* nodes = new node_pool(nodePoolSize);
    //typedef std::tuple<int,int,bool> node_idx; //node id, interval id, from landmark
    typedef boost::heap::pairing_heap< Node*, boost::heap::compare<Node::compare_node> > heap_open_t;
    typedef boost::unordered_set<Node*, Node::NodeHasher, Node::eqnode> hashtable_t;
    heap_open_t open;
    hashtable_t visited;
    //std::unordered_map<int, Node> close;
    //boost::unordered_map<node_idx, Node> close;
    //std::list<Node> open;
    //std::unordered_map<int, std::pair<double, bool>> visited;
    //boost::unordered_map<node_idx, std::pair<double, bool>> visited; //g, visited
    std::map<std::pair<int, int>, vector<Interval>> constraints;//stores sets of constraints associated with moves
    std::map<int, vector<Interval>> collision_intervals;//stores sets of collision intervals associated with cells
    std::vector<Move> landmarks;
    Path path;
    int call_find_pp=0;
    int node_open=0;
    int node_exp=0;
    int node_idx=1;
    void prt_landmarks();
    void prt_cons();
    void prt_node(Node n);
    void prt_nodes(vector<Node*> nodes);
    void prt_solution(std::vector<Node> nodes);
void prt_intervals();
void prt_constraints(list<Constraint> constraints);
void prt_info(){
cout<<"called find_partial_path: "<<call_find_pp<<endl;
cout<<"node opened: "<<node_open<<endl;
cout<<"node expanded: "<<node_exp<<endl;

}
};

#endif // SIPP_H
