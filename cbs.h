#ifndef CBS_H
#define CBS_H
#include <chrono>
#include <math.h>
#include <stack>
#include "structs.h"
#include "map.h"
#include "task.h"
#include "edgeSpliter.h"
//#include "config.h"
#include "sipp.h"
#include "heuristic.h"
#include "simplex/simplex.h"
#include "simplex/pilal.h"

class CBS
{
  public:
    CBS() {}
    Solution find_solution(Map &map, const Task &task, const Config &cfg);
    bool init_root(Map &map, const Task &task);
    std::list<Constraint> get_constraints(CBS_Node *node, int agent_id = -1);
    //std::list<Constraint> merge_constraints(std::list<Constraint> constraints);
    bool validate_constraints(std::list<Constraint> constraints, int agent);
    bool check_positive_constraints(std::list<Constraint> constraints, Constraint constraint);
    Conflict check_paths(const sPath &pathA, const sPath &pathB);
    Conflict check_paths_ori(const sPath &pathA, const sPath &pathB);
    bool check_conflict(Move move1, Move move2);
    bool check_conflict_ori(Move move1, Move move2);
    double get_hl_heuristic(const std::list<Conflict> &conflicts);
    std::vector<Conflict> get_all_conflicts(const std::vector<sPath> &paths, int id);
    list<Constraint> get_constraint(int agent, Move move1, Move move2);
    list<Constraint> get_wait_constraint(int agent, Move move1, Move move2);
    Interval binary_wait_search_constraint(int agent, Move move1, Move move2);
    void find_new_conflicts(Map &map, const Task &task, CBS_Node &node, std::vector<sPath> &paths, const sPath &path,
        const std::list<Conflict> &conflicts, const std::list<Conflict> &semicard_conflicts, const std::list<Conflict> &cardinal_conflicts,
        int &low_level_searches, int &low_level_expanded);
    double get_cost(CBS_Node node, int agent_id);
    std::vector<sPath> get_paths(CBS_Node *node, unsigned int agents_size);
    Conflict get_conflict(std::list<Conflict> &conflicts);
    Move find_sub_conflict(Move m1,Move m2,CBS_Node *node);
    Move split_conf_move(Move m1,Move m2, CBS_Node *node);
    Interval binary_search_constraint(int agent, Move move1, Move move2);
    Conflict modify_conflict(Conflict conflict,CBS_Node *node);
    Move modify_move(Move move,int new_id);
    //bool validNewNode(Vector2D X1,Vector2D X2,Vector2D New);
    //typedef std::pair<Map_delta,Map_delta> Map_delta_pair;
    typedef std::pair<Vector2D,Vector2D> node_pair;
    typedef std::pair<int,int> edge;

    //std::pair<Vector2D,double> findAngle(edge e1, edge e2);
    void split_edge(Conflict conflict, std::vector<sPath> paths, Map_deltas &deltasR, Map_deltas &deltasL);
    //node_pair newNodeMoving(Conflict conflict);
    Constraint get_split_constraint(int agent, Move move1, Move move2);
    void gen_new_map(CBS_Node *node);
    void gen_original_map(CBS_Node *node);
    int id2ind(int v1, int v2, int agent);
    void prt_move(Move m1);
    void prt_constraint(Constraint c);
    void prt_constraints(std::list<Constraint> constraints);
    void prt_conflict(Conflict conflict);
    void prt_conflicts(list<Conflict> conflicts);
    void prt_path(sPath p);
    void prt_paths(std::vector<sPath> paths);
    void prt_map_deltas(Map_deltas R,Map_deltas L);
    void prt_map_deltas_aux(Map_deltas md);
    void saveCT(const string &fileName,CBS_Node *goal_node,unsigned int agent_num);
    void printBT_aux();
    void printBT(const string& prefix, const int node_id, bool isLeft);
    void check_collison(CBS_Node *node);
    void prt_history(CBS_Node *node);
    bool validate_path(list<Constraint> constraints,sPath p);

    void post_check(vector<sPath> Paths);
    Vector2D ind2Vec(int nodeId);
    CBS_Tree tree;
    SIPP planner;
    Solution solution;
    Heuristic h_values;
    Config config;
    Map* map;
    Map* original;
    typedef boost::unordered_map<int,CBS_Node_aux*> tree_aux;
    tree_aux tree_info;

};

#endif // CBS_H
