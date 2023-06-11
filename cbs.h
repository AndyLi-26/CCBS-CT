#ifndef CBS_H
#define CBS_H
#include <chrono>
#include <math.h>
#include <stack> 
#include "structs.h"
#include "map.h"
#include "task.h"
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
    bool validate_constraints(std::list<Constraint> constraints, int agent);
    bool check_positive_constraints(std::list<Constraint> constraints, Constraint constraint);
    Conflict check_paths(const sPath &pathA, const sPath &pathB);
    Conflict check_paths_ori(const sPath &pathA, const sPath &pathB);
    bool check_conflict(Move move1, Move move2);
    bool check_conflict_ori(Move move1, Move move2);
    double get_hl_heuristic(const std::list<Conflict> &conflicts);
    std::vector<Conflict> get_all_conflicts(const std::vector<sPath> &paths, int id);
    typedef std::pair<std::list<Constraint>,bool> con_return;
    con_return get_constraint(int agent, Move move1, Move move2, CBS_Node *node);
    CBS::con_return get_wait_constraint(int agent, Move move1, Move move2);
    void find_new_conflicts(Map &map, const Task &task, CBS_Node &node, std::vector<sPath> &paths, const sPath &path,
        const std::list<Conflict> &conflicts, const std::list<Conflict> &semicard_conflicts, const std::list<Conflict> &cardinal_conflicts,
        int &low_level_searches, int &low_level_expanded);
    double get_cost(CBS_Node node, int agent_id);
    std::vector<sPath> get_paths(CBS_Node *node, unsigned int agents_size);
    Conflict get_conflict(std::list<Conflict> &conflicts);
    Move find_sub_conflict(Move m1,Move m2,CBS_Node *node);
    Move split_conf_move(Move m1,Move m2, CBS_Node *node);
    Conflict modify_conflict(Conflict conflict,CBS_Node *node);
    Move modify_move(Move move,int new_id);
    bool validNewNode(Vector2D X1,Vector2D X2,Vector2D New);
    //typedef std::pair<Map_delta,Map_delta> Map_delta_pair;
    typedef std::pair<Vector2D,Vector2D> node_pair;
    typedef std::pair<int,int> edge;

    int id2ind(int v1, int v2, int agent);

    void post_check(vector<sPath> Paths);
    Vector2D ind2Vec(int nodeId);
    CBS_Tree tree;
    SIPP planner;
    Solution solution;
    Heuristic h_values;
    Config config;
    Map* map;
    Map* original;

};

#endif // CBS_H
