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
      ("solution_file",po::value<std::string>()->default_value(""),"solution file")
      ("result_file,o",po::value<std::string>()->default_value(""),"result file")
      ("HI_h,h",po::value<int>()->default_value(0),"HI level heutistic, 0: none, 1:simplex model, 2: count" )
      ("focal_weigth",po::value<float>()->default_value(1.0),"focal weright")
      //("precision",po::value<float>()->default_value(0.00001),"the precision used in math calculation and number comparison")
      //("resolution",po::value<float>()->default_value(0.01),"the resolution of the nodes on the edge")
      ("debug",po::value<int>()->default_value(0),"debug information")
      ("connectdness",po::value<int>()->default_value(2))
      ("agent_size,a",po::value<float>()->default_value(4.5))
      ("agent_num",po::value<int>()->required(),"number of agent")
      ("timelimit",po::value<int>()->default_value(30))
      ("Cardinal","use cardinal for choose conflict")
      ("DS","use Disjoint_splitting")
      ("cons_reason","use minimum clearance time reasoning")
      ("ES","use edge split")
      ("extra_info",po::value<int>()->default_value(-1),"task index");
      
    po::variables_map temp;
    po::store(po::parse_command_line(argc,argv,desc),temp);
    const po::variables_map vm=temp;

    config.agent_size = vm["agent_size"].as<float>();
    config.connectdness = vm["connectdness"].as<int>();
    string fmap=vm["map"].as<string>();
    Map map = Map(config.agent_size,  config.connectdness);
    map.get_map(fmap);
    cout<<"read map success"<<endl;
    map.pre_process();
    
    config.agent_num=vm["agent_num"].as<int>();
    string ftask=vm["tasks"].as<string>();
    Task task(config.agent_num,fmap[fmap.size()-1]=='d');
    task.get_task(ftask);
    if(map.is_roadmap())
      task.make_ij(map);
    else
      task.make_ids(map.get_width());
    cout<<"read task success"<<endl;
    config.F_debug_info=vm["debug_info"].as<string>();
    config.F_solution=vm["solution_file"].as<string>();
    config.F_exp=vm["result_file"].as<string>();
    config.hlh_type=vm["HI_h"].as<int>();
    config.focal_weight=vm["focal_weigth"].as<float>();
    //precision=vm["precision"].as<float>();
    //resolution=vm["resolution"].as<float>();
    cout<<"agent_size"<<config.agent_size<<endl<<flush;
    cout<<"precision "<<CN_PRECISION<<"resolution "<<CN_RESOLUTION<<endl<<flush;


    config.timelimit=vm["timelimit"].as<int>();
    config.debug = vm["debug"].as<int>();

    bool card=vm.count("Cardinal");
    bool DS=vm.count("DS");
    bool ES=vm.count("ES");
    bool CR=vm.count("cons_reason");
    config.use_cardinal=card;
    config.use_disjoint_splitting=DS;
    config.use_edge_split=ES;
    config.cons_reason=CR;
    
    cout<<"agents_num="<<task.get_agent_num()<<endl;
    task.prt_agents();
    CBS cbs;
    Solution solution = cbs.find_solution(map, task, config);
    logger log(config,solution);
    auto found = solution.found?"true":"false";
    auto Use_edge = config.use_edge_split?"true":"false";
    auto Cons_reason =config.cons_reason?"true":"false";
    int task_ind=vm["extra_info"].as<int>();

    std::cout<< "Soulution found: " << found <<"\nUse Edge Splitting: "<< Use_edge <<
      "\nconstraint reasoning: "<<Cons_reason<<
      "\nRuntime: "<<solution.time.count() << "\nMakespan: " << solution.makespan << "\nFlowtime: " << solution.flowtime<< "\nInitial Cost: "<<solution.init_cost<< "\nCollision Checking Time: " << solution.check_time
      << "\nHL expanded: " << solution.high_level_expanded << "\nLL searches: " << solution.low_level_expansions << "\nLL expanded(avg): " << solution.low_level_expanded << std::endl;
    std::cout<<"agent_size: "<<config.agent_size<<std::endl;
    std::cout<<"introduce new node: "<<map.get_new_node_num()<<std::endl;
    if (config.F_solution!="")
    {
      log.write_to_log_summary();
      log.write_to_log_path(map);
    }
    log.write_exp_result(task_ind);
  }
  else
  {
    std::cout<<"Error! Not enough input parameters are specified!\n";
  }
  return 0;
}
