#include "logger.h"
void logger::write_to_log_summary(const Solution &solution)
{
  ofstream fw(config.F_result, std::ios::out | std::ios::app);
  if (!fw.is_open()){
    FAIL("did not open result file properly");
  }
  fw<<solution.time.count()<<","
    <<solution.flowtime<<","
    <<solution.makespan<<","
    <<endl;
  
  fw.close();
}

void logger::write_to_log_path(const Solution &solution, const Map &map)
{
  ofstream fw(config.F_result, std::ios::out);
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

void logger::write_nodes(const Map &map)
{

}
