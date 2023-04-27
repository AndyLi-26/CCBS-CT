#include "logger.h"
void logger::write_to_log_summary()
{
  ofstream fw(config.F_solution, std::ios::out | std::ios::app);
  if (!fw.is_open()){
    FAIL("did not open solution file properly");
  }
  fw<<solution.time.count()<<","
    <<solution.flowtime<<","
    <<solution.makespan<<","
    <<endl;
  
  fw.close();
}

void logger::write_to_log_path(const Map &map)
{
  ofstream fw(config.F_solution, std::ios::out);
  if (!fw.is_open()){
    FAIL("did not open result file properly");
  }
  for (int i=0;i<int(solution.paths.size());i++)
  {
    fw<<i<<","<<solution.paths[i].cost<<"\n";
    auto iter = solution.paths[i].nodes.begin();
    auto it = solution.paths[i].nodes.begin();

    while(iter != std::prev(solution.paths[i].nodes.end()))
    {
      fw<<map.get_i(it->id)<<","<<map.get_j(it->id)<<",";
      iter++;
      fw<<map.get_i(iter->id)<<","<<map.get_j(iter->id)<<",";
      fw<<iter->g-it->g<<","<<endl;
      it++;
    }

  }
  fw.close();
}
void logger::write_exp_result(int task_ind)
{
  ofstream fw(config.F_exp, std::ios::out | std::ios::app);
  if (!fw.is_open()){
    FAIL("did not open result file properly");
  }
  string ES=config.use_edge_split ? "ES":"0";
  string CR=config.cons_reason ? "CR":"0";
  string DS=config.use_disjoint_splitting? "DS":"0";
  int f=solution.found? 1:0;
  fw<<config.agent_num<<","
    <<config.agent_size<<","
    <<task_ind<<","
    <<solution.time.count()<<","
    <<f<<","
    <<solution.init_cost<<","
    <<solution.flowtime<<","
    <<solution.makespan<<","
    <<solution.high_level_expanded<<","
    <<solution.new_node<<","
    <<endl;

}
void logger::write_nodes(const Map &map)
{

}
