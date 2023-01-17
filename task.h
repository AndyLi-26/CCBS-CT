#ifndef TASK_H
#define TASK_H

#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>
#include "tinyxml2.h"
#include "structs.h"
#include "const.h"
#include <fstream>
#include "map.h"
#include <cstring>
#include <string>
#include <boost/tokenizer.hpp>
class Task
{
  private:
    std::vector<Agent> agents;
    bool isGrid;
  public:
    void get_task(string FileName);
    unsigned int get_agents_size() const { return agents.size(); }
    void make_ids(int width);
    void make_ij(const Map &map);
    int agent_num;
    Agent get_agent(int id) const;
    void print_task()
    {
       for(auto agent:agents)
            std::cout<<"<agent start_i=\""<<agent.start_i<<"\" start_j=\""<<agent.start_j<<"\" goal_i=\""<<agent.goal_i<<"\" goal_j=\""<<agent.goal_j<<"\"/>\n";
    }
    Task(int _agent_num, bool isGrid);
	int get_agent_num() {return agents.size();}
	void prt_agents();
};

#endif // TASK_H
