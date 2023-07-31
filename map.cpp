#include "map.h"
Map::Map(Map *m) 
{ 
	grid=m->grid;
	nodes=m->nodes;
	height=m->height;
	width=m->width;
	size=m->size;
	init_node_num=m->init_node_num;
	for (int i=0;i<m->valid_moves.size();++i){
		std::vector<Node> neighbors;
		for (int j=0;j<m->valid_moves[i].size();++j){
			Node new_node;
			new_node.id=m->valid_moves[i][j].id;
			new_node.i=m->valid_moves[i][j].i;
			new_node.j=m->valid_moves[i][j].j;
			neighbors.push_back(new_node);
		}
		valid_moves.push_back(neighbors);
	}
	
	for (gNode original:m->nodes){
		gNode temp;
		temp.i=original.i;
		temp.j=original.j;
		nodes.push_back(temp);
	}
	
}

bool Map::equal(Map *m)
{
	if (grid!=m->grid)
		return false;
	if (height!=m->height)
		return false;
	if (width!=m->width)
		return false;
	if (init_node_num!=m->init_node_num)
		return false;
	for (int i=0;i<m->valid_moves.size();++i){
		for (int j=0;j<m->valid_moves[i].size();++j){
			if (valid_moves[i][j].id!=m->valid_moves[i][j].id ||
			!eq(valid_moves[i][j].i, m->valid_moves[i][j].i) ||
			!eq(valid_moves[i][j].j, m->valid_moves[i][j].j) ){
				std::cout<<"this:"<<std::endl;
				prt_validmoves();
				std::cout<<"ori"<<std::endl;
				m->prt_validmoves();
				return false;
			}
		}
	}
	return true;
}

void Map::get_map(string FileName)
{
  string line;
  ifstream myfile (FileName.c_str());
  if (FileName[FileName.size()-1]=='d'){
    map_is_roadmap=false;
    get_grid(FileName);
  }
  else if(FileName[FileName.size()-1]=='h'){
    map_is_roadmap=true;
    get_roadmap(FileName);
  }
  else{
    FAIL("File name has unrecognisable extension, please ensure the file end with .grid for grid or .graph for graph(roadmap)");
  }
}

double Map::get_i(int id) const
{
  if(!map_is_roadmap)
    return int(id/width);
  else
    return nodes.at(id).i;
}

Vector2D Map::get_coord(int id) const
{
  if (!map_is_roadmap)
    return Vector2D(int(id/width),int(id%width));
  else
    return Vector2D(nodes.at(id).i,nodes.at(id).j);
}

double Map::get_dist(int id1, int id2) const
{
  if (!map_is_roadmap) return -1;
  double dx(pow((nodes[id1].i-nodes[id2].i),2));
  double dy(pow((nodes[id1].j-nodes[id2].j),2));
  return sqrt(dx+dy);
}

double Map::get_min_clear_t(Move m1, int s2)
{
  cacheTable::iterator it;
  clearT_ind tempInd=make_pair(make_pair(m1.id1,m1.id2),s2);
  it=min_clearT.find(tempInd);
  if (it!=min_clearT.end()) //check cache
  {
    return it->second;
  }
  if (m1.id2==s2) // special case
  {
    double result = get_min_clear_t_sameDestination(m1.id2,m1.id1);
    //assert(false);
    min_clearT[tempInd]=result;
    return result;
  }

  if (s2==m1.id1 || gt(dis_P_l(s2,make_pair(get_coord(m1.id1),get_coord(m1.id2))),d))  // if the node has a larger distance than 2r between line and node, then no need to do futher calc
  {
    min_clearT[tempInd]=-1;
    return -1;
  }

  if (!in_path(get_coord(s2),make_pair(get_coord(m1.id1),get_coord(m1.id2))))  // if the node is in 2r distance, than check if it is in path
  {
    min_clearT[tempInd]=-1;
    return -1;
  }

  //general cases
  double min_t(-1),t_pre,t,offset,t1,t2;
  for (Node n:get_valid_moves(s2))
  {
    Vector2D A(get_coord(s2)),C(get_coord(n.id));
    Vector2D B(get_coord(m1.id1)), D(get_coord(m1.id2));
    if(in_path(B,make_pair(A,C)))
    {
      min_t=CN_INFINITY;
      continue;
    }
    Line l1(make_pair(B,D)), l2(make_pair(A,C));
    Vector2D A_(find_intersection(l1,l2));
    double dA_A_(A_.dis(A)),dD_A_(D.dis(A_));

    Vector2D v(B-A_);
    Vector2D unitV(v/v.mod());
    double quad_a(unitV.i*unitV.i+unitV.j*unitV.j);
    double dx(C.i-A_.i),dy(C.j-A_.j);
    double quad_b(-2*dx*unitV.i-2*dy*unitV.j);
    double quad_c(dx*dx+dy*dy-d*d);
    double delta=quad_b*quad_b-4*quad_a*quad_c;
    if (ge(delta,0))
    {
      
      if (eq(delta,0)) delta=0;
      t1=(-quad_b-sqrt(delta))/(2*quad_a);
      t2=(-quad_b+sqrt(delta))/(2*quad_a);
      if (gt(t1,B.dis(A_)))
        t_pre=CN_INFINITY;
      else if (gt(t2,0))
        t_pre=t2+C.dis(A_);
      else
      {
        double a(B.dis(C)),b(C.dis(A_)),c(B.dis(A_));
        double cos_theta((b*b+c*c-a*a)/(2*b*c));
        if (eq(cos_theta,-1)){ //zero
          min_t=0;
          break;
        }
        t_pre=cos2min_t(cos_theta);
      }
    }
    else{
      double a(B.dis(C)),b(C.dis(A_)),c(B.dis(A_));
      double cos_theta((b*b+c*c-a*a)/(2*b*c));
      if (eq(cos_theta,-1)){ //zero
        min_t=0;
        break;
      }
      t_pre=cos2min_t(cos_theta);
    }
    t=t_pre+dA_A_-dD_A_;
    if (le(t,0))
    {
      min_t=0;
      break;
    }
    if (min_t==-1 || t<min_t){
      min_t=t;
    }
  }
  min_clearT[tempInd]=min_t;
  return min_t;
}

bool Map::in_path(Vector2D p, Line l)
{
  //change from two point form to vector form
  Vector2D A(l.first);
  Vector2D V(l.second-l.first);
  double dis(V.mod());
  Vector2D unitv(V/dis); 

  double vx(unitv.i), vy(unitv.j);
  double dx(p.i-l.first.i), dy(p.j-l.first.j);

  double a(vx*vx+vy*vy), b(-2*dx*vx-2*dy*vy), c(dx*dx + dy*dy - d*d);
  double delta(b*b-4*a*c);
  if(le(delta,0))
    return false;
  double t((-b-sqrt(delta))/2*a);

  if (gt(t,0) && lt(t,dis))
    return true;
  else 
    return false;
}

bool Map::node_in_path(int p, pair<int,int> l)
{
  Vector2D P(get_coord(p));
  Vector2D A(get_coord(l.first)), B(get_coord(l.second));
  return in_path(P,make_pair(A,B));
}

double Map::get_min_clear_t_sameDestination(int main_n, int enter_n)
{
  //cout<<"main_n: "<<main_n<<"-->"<<enter_n<<endl;
  if (main_n>=init_node_num) return -1;
  double min_t(-1),t,t1,t2;
  for (Node exit:valid_moves[main_n]){
    //if (valid_moves[main_n][enter_n].id==exit.id) continue;
    if (enter_n==exit.id) continue;
    int nA(main_n),nB(enter_n),nC(exit.id);
    Vector2D A(get_coord(nA)),B(get_coord(nB)),C(get_coord(nC));
    //check if C is in the path
    Vector2D v(B-A);
    Vector2D unitV(v/v.mod());
    double quad_a(unitV.i*unitV.i+unitV.j*unitV.j);
    double dx(C.i-A.i),dy(C.j-A.j);
    double quad_b(-2*dx*unitV.i-2*dy*unitV.j);
    double quad_c(dx*dx+dy*dy-d*d);
    double delta=quad_b*quad_b-4*quad_a*quad_c;
    if (ge(delta,0))
    {
      
      if (eq(delta,0)) delta=0;
      t1=(-quad_b-sqrt(delta))/(2*quad_a);
      t2=(-quad_b+sqrt(delta))/(2*quad_a);
      if (gt(t1,get_dist(nB,nA)))
        t=CN_INFINITY;
      else if (gt(t2,0))
        t=t2+get_dist(nA,nC);
      else
      {
        double a(get_dist(nB,nC)),b(get_dist(nA,nC)),c(get_dist(nA,nB));
        //cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<endl;
        double cos_theta((b*b+c*c-a*a)/(2*b*c));
        if (eq(cos_theta,-1)) //zero
          return 0;
        t=cos2min_t(cos_theta);
      }
    }
    else
    {
      double a(get_dist(nB,nC)),b(get_dist(nA,nC)),c(get_dist(nA,nB));
      //cout<<"nA: "<<nA<<" nB: "<<nB<<" nC: "<<nC<<endl;
      //cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<endl;
      double cos_theta((b*b+c*c-a*a)/(2*b*c));
      if (eq(cos_theta,-1)){ //zero
        return 0;
      }
      t=cos2min_t(cos_theta);
    }
    if (min_t==-1 || t<min_t){
      min_t=t;
    }
  }
  return min_t;
}

double Map::cos2min_t(double cos_theta)
{
  //double temp((8*agent_size*agent_size)/(1-cos_theta));

  //double t=(sqrt(temp));
  //cout<<"cos theta: "<<cos_theta<<" retval: "<<sqrt((8*agent_size*agent_size)/(1-cos_theta))<<endl;
  return sqrt((8*agent_size*agent_size)/(1-cos_theta));
}

double Map::get_j(int id) const
{
  if(!map_is_roadmap)
    return int(id%width);
  else
    return nodes.at(id).j;
}

void Map::get_grid(string FileName)
{


}

/*
   bool Map::get_grid(string FileName)
   {

   tinyxml2::XMLElement *root = nullptr, *map = nullptr, *element = nullptr, *mapnode = nullptr;

   std::string value;
   std::stringstream stream;
   bool hasGridMem(false), hasGrid(false), hasHeight(false), hasWidth(false);

   tinyxml2::XMLDocument doc;
   if (doc.LoadFile(FileName) != tinyxml2::XMLError::XML_SUCCESS)
   {
   std::cout << "Error opening XML file!" << std::endl;
   return false;
   }
   root = doc.FirstChildElement(CNS_TAG_ROOT);
   if (!root)
   {
   std::cout << "Error! No '" << CNS_TAG_ROOT << "' tag found in XML file!" << std::endl;
   return false;
   }
   map = root->FirstChildElement(CNS_TAG_MAP);
   if (!map)
   {
   std::cout << "Error! No '" << CNS_TAG_MAP << "' tag found in XML file!" << std::endl;
   return false;
   }

   for (mapnode = map->FirstChildElement(); mapnode; mapnode = mapnode->NextSiblingElement())
   {
   element = mapnode->ToElement();
   value = mapnode->Value();
   std::transform(value.begin(), value.end(), value.begin(), ::tolower);

   stream.str("");
   stream.clear();
   stream << element->GetText();

   if (!hasGridMem && hasHeight && hasWidth)
   {
   grid.resize(height);
   for (int i = 0; i < height; ++i)
   grid[i].resize(width);
   hasGridMem = true;
   }

   if (value == CNS_TAG_HEIGHT)
   {
   if (hasHeight)
   {
   std::cout << "Warning! Duplicate '" << CNS_TAG_HEIGHT << "' encountered." << std::endl;
   std::cout << "Only first value of '" << CNS_TAG_HEIGHT << "' =" << height << "will be used."
   << std::endl;
   }
   else
   {
   if (!((stream >> height) && (height > 0)))
   {
   std::cout << "Warning! Invalid value of '" << CNS_TAG_HEIGHT
   << "' tag encountered (or could not convert to integer)." << std::endl;
   std::cout << "Value of '" << CNS_TAG_HEIGHT << "' tag should be an integer >=0" << std::endl;
   std::cout << "Continue reading XML and hope correct value of '" << CNS_TAG_HEIGHT
   << "' tag will be encountered later..." << std::endl;
   }
   else
   hasHeight = true;
   }
   }
   else if (value == CNS_TAG_WIDTH)
   {
if (hasWidth)
{
  std::cout << "Warning! Duplicate '" << CNS_TAG_WIDTH << "' encountered." << std::endl;
  std::cout << "Only first value of '" << CNS_TAG_WIDTH << "' =" << width << "will be used." << std::endl;
}
else
{
  if (!((stream >> width) && (width > 0)))
  {
    std::cout << "Warning! Invalid value of '" << CNS_TAG_WIDTH
      << "' tag encountered (or could not convert to integer)." << std::endl;
    std::cout << "Value of '" << CNS_TAG_WIDTH << "' tag should be an integer AND >0" << std::endl;
    std::cout << "Continue reading XML and hope correct value of '" << CNS_TAG_WIDTH
      << "' tag will be encountered later..." << std::endl;

  }
  else
    hasWidth = true;
}
}
else if (value == CNS_TAG_GRID)
{
  int grid_i(0), grid_j(0);
  hasGrid = true;
  if (!(hasHeight && hasWidth))
  {
    std::cout << "Error! No '" << CNS_TAG_WIDTH << "' tag or '" << CNS_TAG_HEIGHT << "' tag before '"
      << CNS_TAG_GRID << "'tag encountered!" << std::endl;
    return false;
  }
  element = mapnode->FirstChildElement();
  while (grid_i < height)
  {
    if (!element)
    {
      std::cout << "Error! Not enough '" << CNS_TAG_ROW << "' tags inside '" << CNS_TAG_GRID << "' tag."
        << std::endl;
      std::cout << "Number of '" << CNS_TAG_ROW
        << "' tags should be equal (or greater) than the value of '" << CNS_TAG_HEIGHT
        << "' tag which is " << height << std::endl;
      return false;
    }
    std::string str = element->GetText();
    std::vector<std::string> elems;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, ' '))
      elems.push_back(item);
    grid_j = 0;
    int val;
    if (elems.size() > 0)
      for (grid_j = 0; grid_j < width; ++grid_j)
      {
        if (grid_j == int(elems.size()))
          break;
        stream.str("");
        stream.clear();
        stream << elems[grid_j];
        stream >> val;
        grid[grid_i][grid_j] = val;
      }

    if (grid_j != width)
    {
      std::cout << "Invalid value on " << CNS_TAG_GRID << " in the " << grid_i + 1 << " " << CNS_TAG_ROW
        << std::endl;
      return false;
    }
    ++grid_i;
    element = element->NextSiblingElement();
  }
}
}
if (!hasGrid) {
  std::cout << "Error! There is no tag 'grid' in xml-file!\n";
  return false;
}
size = width*height;
std::vector<Step> moves;
valid_moves.resize(height*width);
if(connectedness == 2)
  moves = {{0,1}, {1,0}, {-1,0},  {0,-1}};
else if(connectedness == 3)
  moves = {{0,1}, {1,1}, {1,0},  {1,-1},  {0,-1},  {-1,-1}, {-1,0}, {-1,1}};
else if(connectedness == 4)
  moves = {{0,1}, {1,1}, {1,0},  {1,-1},  {0,-1},  {-1,-1}, {-1,0}, {-1,1},
    {1,2}, {2,1}, {2,-1}, {1,-2}, {-1,-2}, {-2,-1}, {-2,1},  {-1,2}};
else
moves = {{0,1},   {1,1},   {1,0},   {1,-1},  {0,-1},  {-1,-1}, {-1,0}, {-1,1},
  {1,2},   {2,1},   {2,-1},  {1,-2},  {-1,-2}, {-2,-1}, {-2,1}, {-1,2},
  {1,3},   {2,3},   {3,2},   {3,1},   {3,-1},  {3,-2},  {2,-3}, {1,-3},
  {-1,-3}, {-2,-3}, {-3,-2}, {-3,-1}, {-3,1},  {-3,2},  {-2,3}, {-1,3}};
for(int i = 0; i < height; i++)
for(int j = 0; j < width; j++)
{
  std::vector<bool> valid(moves.size(), true);
  for(unsigned int k = 0; k < moves.size(); k++)
    if((i + moves[k].i) < 0 || (i + moves[k].i) >= height || (j + moves[k].j) < 0 || (j + moves[k].j) >= width
        || cell_is_obstacle(i + moves[k].i, j + moves[k].j)
        || !check_line(i, j, i + moves[k].i, j + moves[k].j))
      valid[k] = false;
  std::vector<Node> v_moves = {};
  for(unsigned int k = 0; k < valid.size(); k++)
    if(valid[k])
      v_moves.push_back(Node((i + moves[k].i)*width + moves[k].j + j, 0, 0, i + moves[k].i, j + moves[k].j));
  valid_moves[i*width+j] = v_moves;
}
return true;
}
*/

void Map::get_roadmap(string FileName)
{
  ifstream myfile(FileName);
  int n_num;
  if(!(myfile>>n_num)){
    FAIL("Invalid node number, expecting an int");
  }

  double x,y;
  for (int i=0;i<n_num;++i)
  {
    if (!(myfile>>x>>y)){
      FAIL("Invalid coord, expecting 2 floats");
    }
    gNode node(x,y);
    //node.agent.insert(-1);
    nodes.push_back(node);
    //ori_node_ind ind(std::make_pair(x,y));
    //ori_node_table[ind]=nodes.size()-1;
  }
  nodes_num=n_num;
  init_node_num=n_num;

  int edges,id1,id2;
  if(!(myfile>>edges)){
    FAIL("Invalid edge number, expecting an int");
  }
  for (int i=0;i<edges;i++)
  {
    if (!(myfile>>id1>>id2)){
      FAIL("Invalid edge, expecting 2 ints");
    }
    nodes[id1].neighbors.push_back(id2);
  }
  for(gNode cur:nodes)
  {
    Node node;
    std::vector<Node> neighbors;
    neighbors.clear();
    for(unsigned int i = 0; i < cur.neighbors.size(); i++)
    {
      node.i = nodes[cur.neighbors[i]].i;
      node.j = nodes[cur.neighbors[i]].j;
      node.id = cur.neighbors[i];
      neighbors.push_back(node);
    }
    valid_moves.push_back(neighbors);
  }
  size = int(nodes.size());
  activated=vector<bool>(size,true);
}

/*
   double Map::min_min_clearT(int node)
   {
   vector<double> v=min_clear_time.at(node);
   return *min_element(v.begin(),v.end());
   }
   */

/*
   void Map::get_roadmap(string FileName)
   {
   string line;
   ifstream myfile(FileName.c_str());
   if (myfile.is_open())
   {
   getline (myfile,line);
   int nodes = atoi ( (line).c_str() );
   for(int i=0;i<nodes;++i)
   {
   getline (myfile,line);
   char_separator<char> sep(",");
   tokenizer< char_separator<char> > tok(line, sep);
   tokenizer< char_separator<char> >::iterator beg=tok.begin();
   double x = stod ( (*beg).c_str() ); // read x
   beg++;
   double y = stod ( (*beg).c_str() ); // read y
   gNode node(x,y);
   node.agent.insert(-1);
   nodes.push_back(node);
   ori_node_ind ind(std::make_pair(x,y));
   ori_node_table[ind]=nodes.size()-1;
   }
   nodes_num=nodes;
   init_node_num=nodes;

   getline(myfile,line);
   int edges = atoi ((*line).c_str());
   for (int i=0;i<edges;i++)
   {
   getline (myfile,line);
   char_separator<char> sep(",");
   tokenizer< char_separator<char> > tok(line, sep);
   tokenizer< char_separator<char> >::iterator beg=tok.begin();
   int id1 = atoi ( (*beg).c_str() ); // read id1
   beg++;
   int id2 = atoi ( (*beg).c_str() ); // read id2
   nodes[id1].neighbors.push_back(id2);
   }

   for(gNode cur:nodes)
   {
   Node node;
   std::vector<Node> neighbors;
   neighbors.clear();
   for(unsigned int i = 0; i < cur.neighbors.size(); i++)
   {
   node.i = nodes[cur.neighbors[i]].i;
   node.j = nodes[cur.neighbors[i]].j;
   node.id = cur.neighbors[i];
   neighbors.push_back(node);
   }
   valid_moves.push_back(neighbors);
   }
   size = int(nodes.size());
   }
   else
   {
   cerr<<"map file:("<<FileName <<") not found"<<endl;
   exit(10);
   }
   }
   */
/*
   bool Map::get_roadmap(string FileName)
   {
   tinyxml2::XMLDocument doc;
   if (doc.LoadFile(FileName) != tinyxml2::XMLError::XML_SUCCESS)
   {
   std::cout << "Error opening XML file!" << std::endl;
   return false;
   }
   tinyxml2::XMLElement *root = 0, *element = 0, *data;
   std::string value;
   std::stringstream stream;
   root = doc.FirstChildElement("graphml")->FirstChildElement("graph");
   for(element = root->FirstChildElement("node"); element; element = element->NextSiblingElement("node"))
   {
   data = element->FirstChildElement();

   stream.str("");
   stream.clear();
   stream << data->GetText();
   stream >> value;
   auto it = value.find_first_of(",");
   stream.str("");
   stream.clear();
   stream << value.substr(0, it);
   double i;
   stream >> i;
   i=fit2grid(i);
   stream.str("");
   stream.clear();
   value.erase(0, ++it);
   stream << value;
   double j;
   stream >> j;
   j=fit2grid(j);
   gNode node(i,j);
   node.agent.insert(-1);
   nodes.push_back(node);
//nodes_table.insert(n);
ori_node_ind ind(std::make_pair(i,j));
ori_node_table[ind]=nodes.size()-1;
}
//prt_nodes();
nodes_num=nodes.size();
init_node_num=nodes_num;

for(element = root->FirstChildElement("edge"); element; element = element->NextSiblingElement("edge"))
{
std::string source = std::string(element->Attribute("source"));
std::string target = std::string(element->Attribute("target"));
source.erase(source.begin(),++source.begin());
target.erase(target.begin(),++target.begin());
int id1, id2;
stream.str("");
stream.clear();
stream << source;
stream >> id1;
stream.str("");
stream.clear();
stream << target;
stream >> id2;
nodes[id1].neighbors.push_back(id2);
}
for(gNode cur:nodes)
{
Node node;
std::vector<Node> neighbors;
neighbors.clear();
for(unsigned int i = 0; i < cur.neighbors.size(); i++)
{
node.i = nodes[cur.neighbors[i]].i;
node.j = nodes[cur.neighbors[i]].j;
node.id = cur.neighbors[i];
neighbors.push_back(node);
}
valid_moves.push_back(neighbors);
}
size = int(nodes.size());
return true;
}
*/
/*
   Vector2D Map::fit2line(Vector2D precise_pos, int node1, int node2)
   {
   cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
   cout<<"n1: "<<node1<<" n2: "<<node2<<endl;
   Vector2D P1,P2;
   int dir(node1>node2);
   cout<<"dir: "<<dir<<endl;
   if (dir)
   {
   P1=get_coord(node2);
   P2=get_coord(node1);
   }
   else
   {
   P1=get_coord(node1);
   P2=get_coord(node2);
   }
   cout<<"P1:"<<P1<<" P2:"<<P2<<endl<<flush;
   Vector2D v=P2-P1;
   Vector2D unitv=v/sqrt(v*v);
   cout<<"unitv: "<<unitv<<endl;
   Vector2D temp = precise_pos-P1;//get_coord(node1);
   double dis=sqrt(temp*temp);

   if (dir)
   {
   double edge_length=sqrt(v*v);
   dis=edge_length-dis;
   }

   cout<<"dis: "<<dis<<endl;
   int rounded=(int) (dis/config.min_dis);
   double temp2=(rounded+dir)*config.min_dis;
   Vector2D newP(P1+unitv*temp2);
   cout<<"================================================================"<<endl;
   cout<<get_coord(node1)<<"--"<<precise_pos<<"--"<<get_coord(node2)<<endl<<"this turned into: "<<newP<<endl<<flush;
   assert(abs(newP.i-precise_pos.i)<config.min_dis && abs(newP.j-precise_pos.j)<config.min_dis);
   return newP;
   }
   */
int Map::add_node(double i, double j, int node1, int node2)
{
  //std::cout<<"new agent:"<<agent<<std::endl;
  new_table::iterator it;

  int n1(node1), n2(node2);
  int prev(node2), swap;
  while (isNewNode(n1))
  {
    auto temp=get_valid_moves(n1);
    swap=n1;
    n1 = temp.at(0).id==prev ? temp.at(1).id : temp.at(0).id;
    prev=swap;
  }

  prev=node1;
  while (isNewNode(n2))
  {
    auto temp=get_valid_moves(n2);
    swap=n2;
    n2 = temp.at(0).id==prev ? temp.at(1).id : temp.at(0).id;
    prev=swap;
  }

  int s_id= (n1<n2) ? n1 : n2;
  int l_id= (n1>n2) ? n1 : n2;

  new_node_ind ind(std::make_tuple(i,j,s_id,l_id));
  //check if exist
  it = new_node_table.find(ind);
  //std::cout<<std::endl<<std::endl<<"info"<<std::endl;
  //prt_ind(ind);
  //std::cout<<"("<<ind.first.first<<","<<ind.first.second<<") a:"<<ind.second<<std::endl;
  //prt_nodes();
  int node_id;
  if (it != new_node_table.end())
  {
    node_id=it->second;
    if (activated.at(node_id))
      return -1;
    else
      return node_id;
  }
  /*
     for (int idx=init_node_num;idx<size;idx++)
     {
     if (!activated.at(idx)) continue;
     if (eq(nodes.at(idx).i,i) && eq(nodes.at(idx).j,j))
     return -1;
     }
     */
  {
    node_id=nodes.size();
    gNode tempnode(i,j);
    nodes.push_back(tempnode);
    this->size+=1;
    //nodes[node_id].agent=agent;
    new_node_table[ind]=node_id;

    std::vector<Node> neighbors;
    Node valid1;
    valid1.i = nodes[node1].i;
    valid1.j = nodes[node1].j;
    valid1.id = node1;
    //valid1.agent.insert(-1);
    neighbors.push_back(valid1);
    Node valid2;
    valid2.i = nodes[node2].i;
    valid2.j = nodes[node2].j;
    //valid2.agent.insert(-1);
    valid2.id = node2;
    neighbors.push_back(valid2);
    valid_moves.push_back(neighbors);
  }
  //add nodes to node list
  nodes[node_id].neighbors.push_back(node1);
  nodes[node_id].neighbors.push_back(node2);
  nodes[node1].neighbors.push_back(node_id);
  nodes[node2].neighbors.push_back(node_id);
  activated.push_back(false);
  return node_id;
}

void Map::print_map()
{
  std::cout<<height<<"x"<<width<<std::endl;
  for(int i = 0; i < height; i++)
  {
    std::cout<<"<row>";
    for(int j = 0; j < width; j++)
      std::cout<<grid[i][j]<<" ";
    std::cout<<"</row>"<<std::endl;
  }
}

void Map::printPPM()
{
  std::cout<<"P3\n"<<width<<" "<<height<<"\n255\n";
  for(int i = 0; i < height; i++)
    for(int j = 0; j < width; j++)
    {
      if(grid[i][j]==1)
        std::cout<<"0 0 0\n";
      else
        std::cout<<"255 255 255\n";
    }
}

bool Map::cell_is_obstacle(int i, int j) const
{
  return (grid[i][j] == CN_OBSTL);
}
/*
   std::vector<Node> Map::get_valid_moves(int id,int agent) const
   {
   std::vector<Node> retval=valid_moves.at(id);
   if (agent==-1)
   return retval;

   for (int i=0;i<retval.size();++i){
   if(retval.at(i).agent!=-1 && retval.at(i).agent!=agent) {
//if (retval.at(i).agent.find(-1)==retval.at(i).agent.end() && retval.at(i).agent.find(agent)==retval.at(i).agent.end()){
//if (retval[i].agent!=-1 && retval[i].agent!=agent)
//if (valid_moves[retval[i].id].size()!=2){
//	std::cout<<"id: "<<id<<", agent: "<<agent<<", i:"<<i<<std::endl;
//	prt_validmoves();
//	assert(false);
//}
if (valid_moves.at(retval[i].id)[0].id==id)
retval[i]=valid_moves[retval[i].id][1];
else
retval[i]=valid_moves[retval[i].id][0];
}
}
return retval;
}
*/
std::vector<Node> Map::get_valid_moves(int id) const
{
  return valid_moves.at(id);
}
bool Map::check_line(int x1, int y1, int x2, int y2)
{
  int delta_x(std::abs(x1 - x2));
  int delta_y(std::abs(y1 - y2));
  if((delta_x > delta_y && x1 > x2) || (delta_y >= delta_x && y1 > y2))
  {
    std::swap(x1, x2);
    std::swap(y1, y2);
  }
  int step_x(x1 < x2 ? 1 : -1);
  int step_y(y1 < y2 ? 1 : -1);
  int error(0), x(x1), y(y1);
  int gap = int(agent_size*sqrt(pow(delta_x, 2) + pow(delta_y, 2)) + double(delta_x + delta_y)/2 - CN_PRECISION);
  int k, num;

  if(delta_x > delta_y)
  {
    int extraCheck = int(agent_size*delta_y/sqrt(pow(delta_x, 2) + pow(delta_y, 2)) + 0.5 - CN_PRECISION);
    for(int n = 1; n <= extraCheck; n++)
    {
      error += delta_y;
      num = (gap - error)/delta_x;
      for(k = 1; k <= num; k++)
        if(cell_is_obstacle(x1 - n*step_x, y1 + k*step_y))
          return false;
      for(k = 1; k <= num; k++)
        if(cell_is_obstacle(x2 + n*step_x, y2 - k*step_y))
          return false;
    }
    error = 0;
    for(x = x1; x != x2 + step_x; x++)
    {
      if(cell_is_obstacle(x, y))
        return false;
      if(x < x2 - extraCheck)
      {
        num = (gap + error)/delta_x;
        for(k = 1; k <= num; k++)
          if(cell_is_obstacle(x, y + k*step_y))
            return false;
      }
      if(x > x1 + extraCheck)
      {
        num = (gap - error)/delta_x;
        for(k = 1; k <= num; k++)
          if(cell_is_obstacle(x, y - k*step_y))
            return false;
      }
      error += delta_y;
      if((error<<1) > delta_x)
      {
        y += step_y;
        error -= delta_x;
      }
    }
  }
  else
  {
    int extraCheck = int(agent_size*delta_x/sqrt(pow(delta_x, 2) + pow(delta_y, 2)) + 0.5 - CN_PRECISION);
    for(int n = 1; n <= extraCheck; n++)
    {
      error += delta_x;
      num = (gap - error)/delta_y;
      for(k = 1; k <= num; k++)
        if(cell_is_obstacle(x1 + k*step_x, y1 - n*step_y))
          return false;
      for(k = 1; k <= num; k++)
        if(cell_is_obstacle(x2 - k*step_x, y2 + n*step_y))
          return false;
    }
    error = 0;
    for(y = y1; y != y2 + step_y; y += step_y)
    {
      if(cell_is_obstacle(x, y))
        return false;
      if(y < y2 - extraCheck)
      {
        num = (gap + error)/delta_y;
        for(k = 1; k <= num; k++)
          if(cell_is_obstacle(x + k*step_x, y))
            return false;
      }
      if(y > y1 + extraCheck)
      {
        num = (gap - error)/delta_y;
        for(k = 1; k <= num; k++)
          if(cell_is_obstacle(x - k*step_x, y))
            return false;
      }
      error += delta_x;
      if((error<<1) > delta_y)
      {
        x += step_x;
        error -= delta_y;
      }
    }
  }
  return true;
}
/*
   void Map::prt_nodes(){
   ori_node_table::iterator it_ori;
   new_node_table::iterator it_new;
   std::cout<<"ori_table"<<std::endl;
   for(auto it_ori=ori_table.begin(); it_ori!=ori_table.end(); ++it_ori){
   std::cout<<"id:"<<it_ori->second;
   std::cout<<"("<<it_ori->first->first<<","<<it_ori->first->second<<")";
   }

   }
   */
void Map::prt_validmoves() const
{
  //std::cout<<"Map:"<<std::endl;
  for (int i=0;i<valid_moves.size();++i){
    string act= activated.at(i) ? "+" : "-";
    auto vecNode=valid_moves[i];
    std::cout<<"|"<<i<<" "<<act<<
      " ("<<nodes[i].i<<","<<nodes[i].j<<") : ";
    for (Node n:vecNode){
      std::cout<<n.id<<", ";
    }
    std::cout<<std::endl;
  }
  std::cout<<"&&&&&&&"<<std::endl<<std::endl;
}

void Map::prt_set(std::set<int> s) const
{
  for (int i:s)
    std::cout<<i<<",";
}

void Map::alter(Map_deltas map_deltas)
{
  for (Map_delta delta:map_deltas){
    int node_id(delta.add_node);
    assert(!activated.at(node_id));
    activated.at(node_id)=true;
    if (node_id==-1)
      return;

    int v1(delta.del_edge.first),v2(delta.del_edge.second);
    valid_moves.at(node_id).at(0).id=v1;
    valid_moves.at(node_id).at(0).i=get_i(v1);
    valid_moves.at(node_id).at(0).j=get_j(v1);
    valid_moves.at(node_id).at(1).id=v2;
    valid_moves.at(node_id).at(1).i=get_i(v2);
    valid_moves.at(node_id).at(1).j=get_j(v2);


    for (int i=0;i<valid_moves[v1].size();++i){
      if (valid_moves[v1][i].id==v2){
        valid_moves[v1][i].i=get_i(node_id);
        valid_moves[v1][i].j=get_j(node_id);
        valid_moves[v1][i].id=node_id;
        //valid_moves[v1][i].agent=nodes[node_id].agent;
        break;
      }
    }

    for (int i=0;i<valid_moves[v2].size();++i){
      if (valid_moves[v2][i].id==v1){
        valid_moves[v2][i].i=get_i(node_id);
        valid_moves[v2][i].j=get_j(node_id);
        valid_moves[v2][i].id=node_id;
        //valid_moves[v2][i].agent=nodes[node_id].agent;
        break;
      }
    }
  }
}

void Map::alter_back(Map_deltas map_deltas)
{
  for (Map_delta delta: map_deltas){
    int node_id(delta.add_node);
    assert(activated.at(node_id));
    activated.at(node_id)=false;
    if (node_id==-1)
      return;

    int v1(delta.del_edge.first),v2(delta.del_edge.second);

    for (int i=0;i<valid_moves[v1].size();++i){
      if (valid_moves[v1][i].id==node_id){
        valid_moves[v1][i].i=get_i(v2);
        valid_moves[v1][i].j=get_j(v2);
        valid_moves[v1][i].id=v2;
        //valid_moves[v1][i].agent=nodes[i].agent;
      }
    }

    for (int i=0;i<valid_moves[v2].size();++i){
      if (valid_moves[v2][i].id==node_id){
        valid_moves[v2][i].i=get_i(v1);
        valid_moves[v2][i].j=get_j(v1);
        valid_moves[v2][i].id=v1;
        //valid_moves[v2][i].agent=nodes[i].agent;
      }
    }
  }
}

int Map::id2ind(int v1, int v2)
{
  std::vector<Node> temp_moves=valid_moves.at(v1);

  for (int i=0;i<temp_moves.size();++i){
    if (temp_moves.at(i).id==v2)
      return i;
  }
  assert(false);
  return -2;
}

Vector2D Map::find_intersection(Line l1, Line l2)
{
  double x1(l1.first.i), y1(l1.first.j), x2(l1.second.i),y2(l1.second.j);
  double x3(l2.first.i), y3(l2.first.j), x4(l2.second.i),y4(l2.second.j);
  double c1(x1*y2-y1*x2), c2(x3*y4-y3*x4);
  double denominator= (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
  double px=( c1*(x3-x4) - (x1-x2)*c2 ) / denominator;
  double py=( c1*(y3-y4) - (y1-y2)*c2 ) / denominator;
  return Vector2D(px,py);

}

double Map::dis_P_l(Vector2D P, Line l)
{
  assert(!(l.first==l.second));
  double x0(P.i),y0(P.j);
  double x1(l.first.i),y1(l.first.j);
  double x2(l.second.i),y2(l.second.j);
  double numerator = abs( (x2-x1)*(y1-y0)-(x1-x0)*(y2-y1) );
  double denominator = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
  return numerator/denominator;
}

void Map::pre_process()
{
  auto dist=[](Vector2D A, Vector2D B) -> double {return sqrt((A-B)*(A-B));};
  for (int i=0; i<size;++i){
    vector<double> min_clearV;
    min_clearV.clear();
    for (int j=0;j<valid_moves[i].size();++j){
      double min_t(-1);
      for (Node exit:valid_moves[i]){
        if (valid_moves[i][j].id==exit.id) continue;
        //calc theta
        Vector2D A(get_coord(i)),B(get_coord(valid_moves[i][j].id)),C(get_coord(exit.id));
        //cout<<"node A: "<<i<<A<<", node B: "<<valid_moves[i][j].id<<B<<" , node C:"<<exit.id<<C<<endl;

        double a(dist(B,C)),b(dist(A,C)),c(dist(A,B));
        double cos_theta((b*b+c*c-a*a)/(2*b*c));
        //cout<<"a: "<<a<<"b: "<<b<<"c: "<<c<<endl;
        //cout<<"cos: "<<cos_theta<<endl;
        if (cos_theta==-1){ //zero
          min_t=0;
          break;
        }
        //if (cos_theta==-1) continue;//inf large
        //find min_t
        double temp((8*agent_size*agent_size)/(1-cos_theta));
        double t=(sqrt(temp));
        //double t(sqrt((2*agent_size*agent_size)/(1-cos_theta)));
        //cout<<"t: "<<t<<endl;
        if (min_t==-1 || t<min_t){
          min_t=t;
        }
      }
      //cout<<"min t:"<<min_t<<endl;
      min_clearV.push_back(min_t);
    }
    min_clear_time.push_back(min_clearV);
  }
  for(int i=0;i<min_clear_time.size();i++)
  {
    cout<<i<<": "<<endl;
    for (double j:min_clear_time.at(i))
      cout<<j<<", ";
    cout<<endl;
  }
}
