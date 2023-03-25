#ifndef MAP_H
#define MAP_H

#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include "tinyxml2.h"
#include "const.h"
#include "structs.h"
#include <boost/unordered_map.hpp>
#include <cstring>
#include <string>
#include <boost/tokenizer.hpp>

class Map
{
  private:
    std::vector<std::vector<int>> grid;
    std::vector<gNode> nodes;
    std::vector<std::vector<Node>> valid_moves;
    int  height, width, size;
    int  connectedness;
    double agent_size;
    bool map_is_roadmap;
    bool check_line(int x1, int y1, int x2, int y2);
    void get_grid(string FileName);
    void get_roadmap(string FileName);
    int nodes_num; //public node limit
    int init_node_num;
    void prt_nodes();
    typedef std::tuple<double,double,int,int> new_node_ind;  //<x,y,n1,n2> -> locate at (x,y), n1->new_node->n2
    typedef std::pair<double, double> ori_node_ind;
    typedef boost::unordered_map<new_node_ind,int> new_table;
    typedef boost::unordered_map<ori_node_ind,int> ori_table;
    new_table new_node_table;
    ori_table ori_node_table;

  public:
    Map(double size, int k){ agent_size = size; connectedness = k; }
    Map(Map *m);
    ~Map(){}
    std::vector<std::vector<double>> min_clear_time; //[main node][enter node]
    double min_min_clearT(int node);
    int  get_size() const { return size; }
    int get_new_node_num() const {return size-init_node_num;}
    int get_init_node_num() const {return init_node_num;}
    void get_map(string FileName);
    bool is_roadmap() const {return map_is_roadmap;}
    bool cell_is_obstacle(int i, int j) const;
    int  get_width() const {return width;}
    gNode get_gNode(int id) const {if(id < int(nodes.size())) return nodes[id]; return gNode();}
    int  get_id(int i, int j) const;
    double get_i (int id) const;
    double get_j (int id) const;
    Vector2D get_coord(int id) const;
    double get_dist(int id1, int id2) const;
    int add_node(double i, double j, int node1, int node2);
    std::vector<Node> get_valid_moves(int id) const;
    double fit2grid(double val){return round(val/CN_PRECISION)*CN_PRECISION;}
    Vector2D fit2line(Vector2D precise_pos, int node1, int node2);
    void print_map();
    void printPPM();
    //void prt_ind(node_index n);
    void prt_set(std::set<int> s) const;
    void prt_validmoves() const;
    void alter(Map_deltas map_delta);
    void alter_back(Map_deltas map_delta);
    void pre_process();
    bool equal(Map *m);
};

#endif // MAP_H
