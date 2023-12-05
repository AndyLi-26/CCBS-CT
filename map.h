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
#include "handyFunc.cpp"
#include <boost/unordered_map.hpp>
#include <cstring>
#include <string>
#include <boost/tokenizer.hpp>

class Map
{
    private:
        vector<std::vector<int>> grid;
        vector<gNode> nodes;
        vector<bool> activated;
        vector<std::vector<Node>> valid_moves;
        int  height, width, size;
        int  connectedness;
        double agent_size;
        double d;
        bool map_is_roadmap;
        bool check_line(int x1, int y1, int x2, int y2);
        void get_grid(string FileName);
        void get_roadmap(string FileName);
        int nodes_num; //public node limit
        int init_node_num;
        void prt_nodes();
        Config config;
        typedef std::tuple<double,double,int,int> new_node_ind;  //<x,y,n1,n2> -> locate at (x,y), n1->new_node->n2
                                                                 //typedef std::pair<double, double> ori_node_ind;
        typedef boost::unordered_map<new_node_ind,int> new_table;
        //typedef boost::unordered_map<ori_node_ind,int> ori_table;
        typedef pair<pair<int,int>,int> clearT_ind;
        typedef boost::unordered_map<clearT_ind,double> cacheTable; //cache for min_clear T
        double get_min_clear_t_sameDestination(int main_n, int enter_n);
        Vector2D find_intersection(Line l1, Line l2);
        double dis_P_l(Vector2D P, Line l);
        bool in_path(Vector2D p, Line l);
        cacheTable min_clearT;
        new_table new_node_table;
        std::vector<std::vector<double>> min_clear_time; //[main node][enter node]
                                                         //ori_table ori_node_table;

    public:
        Map(const Config cfg){
            agent_size = cfg.agent_size;
            d=2*agent_size;
            connectedness = cfg.connectedness;
            config=cfg;
        }
        Map(Map *m);
        ~Map(){}
        int  get_size() const { return size; }
        int get_new_node_num() const {return size-init_node_num;}
        int get_init_node_num() const {return init_node_num;}
        bool isNewNode(int n) const {return n>=init_node_num;}
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
        //cr part
        double get_min_clear_t(Move m1, int s2);
        double fit2grid(double val){return round(val/CN_PRECISION)*CN_PRECISION;}
        Vector2D fit2line(Vector2D precise_pos, int node1, int node2);
        void print_map();
        void printPPM();
        double cos2min_t(double cos_theta);
        int id2ind(int v1,int v2);
        bool node_in_path(int p, pair<int, int> l);
        void prt_set(std::set<int> s) const;
        void prt_validmoves() const;
        void alter(Map_deltas map_delta);
        void alter_back(Map_deltas map_delta);
        bool equal(Map *m);
        void pre_process();

};

#endif // MAP_H
