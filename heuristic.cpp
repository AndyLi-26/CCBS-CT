#include "heuristic.h"

void Heuristic::init(unsigned int size, unsigned int agents)
{
    h_values.clear();
    h_values.resize(size);
    for(unsigned int i = 0; i < size; i++)
        h_values[i].resize(agents, -1);
}

void Heuristic::count(const Map &map, Agent agent)
{
    Node curNode(agent.goal_id, 0, 0, agent.goal_i, agent.goal_j), newNode;
    open.clear();
    open.insert(curNode);
    while(!open.empty())
    {
        curNode = find_min();
        h_values[curNode.id][agent.id] = curNode.g;
        std::vector<Node> valid_moves = map.get_valid_moves(curNode.id);
        for(auto move: valid_moves)
        {
            newNode.i = move.i;
            newNode.j = move.j;
            newNode.id = move.id;
            newNode.g = curNode.g + dist(curNode, newNode);
            if(h_values[newNode.id][agent.id] < 0)
            {
                auto it = open.get<1>().find(newNode.id);
                if(it != open.get<1>().end())
                {
                    if(it->g > newNode.g)
                        open.get<1>().erase(it);
                    else
                        continue;
                }
                open.insert(newNode);
            }
        }
    }
}

void Heuristic::add_node(const Map &map, Map_delta m_del)
{
  if (map.get_size()==h_values.size()) return;
  int new_id(m_del.add_node);
  int n1(m_del.del_edge.first),n2(m_del.del_edge.second);
  int agents=h_values[0].size();
  vector<double> temp_vec;
  double d1(map.get_dist(new_id,n1)),d2(map.get_dist(new_id,n2));
  double h1,h2;
  

  for (int i=0; i<agents;i++)
  {
    h1=h_values[n1][i]+d1;
    h2=h_values[n2][i]+d2;
    temp_vec.push_back(min(h1,h2));
  }
  h_values.push_back(temp_vec);
}

Node Heuristic::find_min()
{
    Node min = *open.begin();
    open.erase(open.begin());
    return min;
}

void Heuristic::prt()
{
  for(int i=0;i<h_values.size();i++)
  {
    cout<<"| "<<i<<":"; 
    for (int j=0;j<h_values[i].size();j++)
    {
      cout<<" "<<h_values[i][j];
    }
    cout<<endl;
  } 
}
