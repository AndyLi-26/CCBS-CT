#include "task.h"
Task::Task(int _agent_num)
{
    agents.clear();
	agent_num=_agent_num;
}

void Task::get_task(string FileName, bool isGrid)
{

  string line;
  ifstream myfile(FileName.c_str());
  if (myfile.is_open())
  {
    getline (myfile,line);
    int tasks = atoi ( (line).c_str() );
    for (int i=0;i<tasks;i++)
    {
      getline (myfile,line);
      char_separator<char> sep(",");
      tokenizer< char_separator<char> > tok(line, sep);
      tokenizer< char_separator<char> >::iterator beg=tok.begin();
      Agent a;
      if (isGrid)
      {
        int x1 = atoi ((*(beg++)).c_str()); // read id1
        int y1 = atoi ((*(beg++)).c_str()); // read id2
        int x2 = atoi ((*(beg++)).c_str());
        int y2 = atoi ((*(beg++)).c_str());
        a.start_i=x1;
        a.start_j=y1;
        a.goal_i=x2;
        a.goal_j=y2;
      }
      else{
        int id1 = atoi ( (*beg).c_str() ); // read id1
        beg++;
        int id2 = atoi ( (*beg).c_str() ); // read id2
        a.start_id = id1;
        a.goal_id=id2;
        a.id = int(agents.size());
        agents.push_back(a);
        if(int(agents.size()) == agent_num)
            break;
      }
    }
  }
  else
  {
    cerr<<"map file:("<<FileName <<") not found"<<endl;
    exit(10);
  }
}
/*
bool Task::get_task(string FileName, int k)
{
    tinyxml2::XMLElement *root = 0, *agent = 0;
    tinyxml2::XMLDocument doc;

    // Load XML File
    if (doc.LoadFile(FileName) != tinyxml2::XMLError::XML_SUCCESS)
    {
        std::cout << "Error opening TASK XML file!" << std::endl;
        return false;
    }

    // Get ROOT element
    root = doc.FirstChildElement(CNS_TAG_ROOT);
    if (!root)
    {
        std::cout << "Error! No '" << CNS_TAG_ROOT << "' tag found in XML file!" << std::endl;
        return false;
    }

    for (agent = root->FirstChildElement(); agent; agent = agent->NextSiblingElement())
    {
        Agent a;
        a.start_i = agent->DoubleAttribute(CNS_TAG_START_I);
        a.start_j = agent->DoubleAttribute(CNS_TAG_START_J);
        a.start_id = agent->IntAttribute(CNS_TAG_START_ID);
        a.goal_i = agent->DoubleAttribute(CNS_TAG_GOAL_I);
        a.goal_j = agent->DoubleAttribute(CNS_TAG_GOAL_J);
        a.goal_id = agent->IntAttribute(CNS_TAG_GOAL_ID);
        a.id = int(agents.size());
        agents.push_back(a);
        if(int(agents.size()) == agent_num)
            break;
    }
    return true;
}
*/
void Task::make_ids(int width)
{
    for(size_t i = 0; i < agents.size(); i++)
    {
        agents[i].start_id = int(agents[i].start_i)*width + int(agents[i].start_j);
        agents[i].goal_id = int(agents[i].goal_i)*width + int(agents[i].goal_j);
        //std::cout<<agents[i].start_i<<" "<<agents[i].start_j<<"  "<<agents[i].goal_i<<" "<<agents[i].goal_j<<"\n";
    }
}

void Task::make_ij(const Map& map)
{
    for(unsigned int i = 0; i < agents.size(); i++)
    {
        gNode start = map.get_gNode(agents[i].start_id), goal = map.get_gNode(agents[i].goal_id);
        agents[i].start_i = start.i;
        agents[i].start_j = start.j;
        agents[i].goal_i = goal.i;
        agents[i].goal_j = goal.j;
    }

}

Agent Task::get_agent(int id) const
{
    if(id >= 0 && id < int(agents.size()))
        return agents[id];
    else
        return Agent();
}

void Task::prt_agents()
{
	for(Agent a:agents){
		std::cout<<a.id<<":"<<a.start_id<<"->"<<a.goal_id<<std::endl;
	}
}
