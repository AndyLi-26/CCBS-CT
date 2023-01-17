#include <iostream>
#include <boost/program_options.hpp>
#include <fstream>
#include "map.h"
#include "task.h"
#include "cbs.h"
#include "logger.h"
#include "structs.h"
int main(int argc, const char *argv[])
{
  if(argc > 2)
  {
    Config config;

    namespace po=boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("map,m", po::value<std::string>()->required(), "input file for map")
      ("tasks,t", po::value<std::string>()->required(), "input file for all the tasks")
      ("debug_info",po::value<std::string>()->default_value(""),"output file for debug info, such as new node, path etc")
      ("result_file,o",po::value<std::string>()->default_value(""),"result data file")
      ("HI_h,h",po::value<int>()->default_value(0),"HI level heutistic, 0: none, 1:simplex model, 2: count" )
      ("focal_weigth",po::value<float>()->default_value(1.0),"focal weright")
      ("precision",po::value<float>()->default_value(0.00001),"the precision used in math calculation and number comparison")
      ("debug",po::value<int>()->default_value(0),"debug information")
      ("connectdness",po::value<int>()->default_value(2))
      ("agent_size,a",po::value<float>()->default_value(4.5))
      ("agent_num",po::value<int>()->required(),"number of agent")
      ("timelimit",po::value<int>()->default_value(30))
      ("Cardinal","use cardinal for choose conflict")
      ("Disjoint_splitting,DS","use Disjoint_splitting")
      ("Edge_split,ES","use edge split");
      
    po::variables_map temp;
    po::store(po::parse_command_line(argc,argv,desc),temp);
    const po::variables_map vm=temp;

    config.agent_size = vm["agent_size"].as<float>();
    config.connectdness = vm["connectdness"].as<int>();
    string fmap=vm["map"].as<string>();
    Map map = Map(config.agent_size,  config.connectdness);
    map.get_map(fmap);
    cout<<"read map success"<<endl;
    
    config.agent_num=vm["agent_number"].as<int>();
    string ftask=vm["tasks"].as<string>();
    Task task(config.agent_num,fmap[fmap.size()-1]=='d');
    task.get_task(ftask);
    if(map.is_roadmap())
      task.make_ij(map);
    else
      task.make_ids(map.get_width());
    std::cout<<"read task "<<endl;
    
    config.F_debug_info=vm["debug_info"].as<string>();
    config.F_result=vm["result_file"].as<string>();
    config.hlh_type=vm["HI_h"].as<int>();
    config.focal_weight=vm["focal_weigth"].as<float>();
    config.precision=vm["precision"].as<float>();
    config.timelimit=vm["timelimit"].as<int>();
    config.debug = vm["debug"].as<int>();

    bool card=vm.count("Cardinal");
    bool DS=vm.count("Disjoint_splitting");
    bool ES=vm.count("Edge_split");
    config.use_cardinal=card;
    config.use_disjoint_splitting=DS;
    config.use_edge_split=ES;
    
    cout<<"finishing reading PO"<<endl<<flush;

    cout<<"agents_num="<<task.get_agent_num()<<endl;
    task.prt_agents();
    CBS cbs;
    Solution solution = cbs.find_solution(map, task, config);
    logger log(config);
    auto found = solution.found?"true":"false";
    auto Use_edge = config.use_edge_split?"true":"false";

    std::cout<< "Soulution found: " << found <<"\nUse Edge Splitting: "<< Use_edge <<
      "\nRuntime: "<<solution.time.count() << "\nMakespan: " << solution.makespan << "\nFlowtime: " << solution.flowtime<< "\nInitial Cost: "<<solution.init_cost<< "\nCollision Checking Time: " << solution.check_time
      << "\nHL expanded: " << solution.high_level_expanded << "\nLL searches: " << solution.low_level_expansions << "\nLL expanded(avg): " << solution.low_level_expanded << std::endl;
    std::cout<<"agent_size: "<<config.agent_size<<std::endl;
    std::cout<<"introduce new node: "<<map.get_new_node_num()<<std::endl;
    log.write_to_log_summary(solution);
    log.write_to_log_path(solution, map);
  }
  else
  {
    std::cout<<"Error! Not enough input parameters are specified!\n";
  }
  return 0;
}
