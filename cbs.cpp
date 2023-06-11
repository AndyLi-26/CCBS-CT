#include "cbs.h"
bool CBS::init_root(Map &map, const Task &task)
{
  CBS_Node root;
  tree.set_focal_weight(config.focal_weight);
  sPath path;
  for(int i = 0; i < int(task.get_agents_size()); i++)
  {
    Agent agent = task.get_agent(i);
    path = planner.find_path(agent, map, {}, h_values);
    if(path.cost < 0)
      return false;
    root.paths.push_back(path);
    root.cost += path.cost;
  }
  root.low_level_expanded = 0;
  root.parent = nullptr;
  root.id = 1;
  root.id_str = "1";
  auto conflicts = get_all_conflicts(root.paths, -1);
  root.conflicts_num = conflicts.size();
  
  for(auto conflict: conflicts)
    if(!config.use_cardinal)
      root.conflicts.push_back(conflict);
    else
    {
     con_return tempA=get_constraint(conflict.agent1, conflict.move1, conflict.move2, &root),tempB=get_constraint(conflict.agent2, conflict.move2, conflict.move1, &root);

      auto pathA = planner.find_path(task.get_agent(conflict.agent1), map, {tempA.first}, h_values);
      auto pathB = planner.find_path(task.get_agent(conflict.agent2), map, {tempB.first}, h_values);
      //conflict.path1 = pathA;
      //conflict.path2 = pathB;
      if(pathA.cost > root.paths[conflict.agent1].cost && pathB.cost > root.paths[conflict.agent2].cost)
      {
        conflict.overcost = std::min(pathA.cost - root.paths[conflict.agent1].cost, pathB.cost - root.paths[conflict.agent2].cost);
        root.cardinal_conflicts.push_back(conflict);
      }
      else if(pathA.cost > root.paths[conflict.agent1].cost || pathB.cost > root.paths[conflict.agent2].cost)
        root.semicard_conflicts.push_back(conflict);
      else
        root.conflicts.push_back(conflict);
    }
  solution.init_cost = root.cost;
  tree.add_node(root);
  return true;
}

bool CBS::check_conflict(Move move1, Move move2)
{
  double r(2*config.agent_size);
  if (move1.t2+r<move2.t1-CN_EPSILON || move2.t2+r<move1.t1-CN_EPSILON)
    return false;
  double startTimeA(move1.t1), endTimeA(move1.t2), startTimeB(move2.t1), endTimeB(move2.t2);
  double m1i1(map->get_i(move1.id1)), m1i2(map->get_i(move1.id2)), m1j1(map->get_j(move1.id1)), m1j2(map->get_j(move1.id2));
  double m2i1(map->get_i(move2.id1)), m2i2(map->get_i(move2.id2)), m2j1(map->get_j(move2.id1)), m2j2(map->get_j(move2.id2));
  Vector2D A(m1i1, m1j1);
  Vector2D B(m2i1, m2j1);
  Vector2D VA((m1i2 - m1i1)/(move1.t2 - move1.t1), (m1j2 - m1j1)/(move1.t2 - move1.t1));
  Vector2D VB((m2i2 - m2i1)/(move2.t2 - move2.t1), (m2j2 - m2j1)/(move2.t2 - move2.t1));
  if(startTimeB > startTimeA)
  {
    A += VA*(startTimeB-startTimeA);
    startTimeA = startTimeB;
  }
  else if(startTimeB < startTimeA)
  {
    B += VB*(startTimeA - startTimeB);
    startTimeB = startTimeA;
  }
  Vector2D w(B - A);
  double c(w*w - r*r);
  if(c < 0)
    return true;

  Vector2D v(VA - VB);
  double a(v*v);
  double b(w*v);
  double dscr(b*b - a*c);
  if(dscr - CN_EPSILON < 0)
    return false;
  double ctime = (b - sqrt(dscr))/a;
  if(ctime >  -CN_EPSILON  && ctime < std::min(endTimeB,endTimeA) - startTimeA +  CN_EPSILON )
    return true;
  return false;
}

//original version
CBS::con_return CBS::get_wait_constraint(int agent, Move move1, Move move2)
{
  list<Constraint> constraint(0);
  Constraint waitCon;
  if (move2.id1==move2.id2) //colliding while waiting
  {
    waitCon=Constraint(agent, move2.t1, move2.t2, move1.id1, -1,move1.id2);
    assert(waitCon.t2-waitCon.t1>CN_EPSILON);
    constraint.push_back(waitCon);
    return make_pair(constraint,false);
  }
  double radius = 2*config.agent_size;
  std::pair<double,double> interval;
  double i0(map->get_i(move2.id1)), j0(map->get_j(move2.id1)), i1(map->get_i(move2.id2)), j1(map->get_j(move2.id2)), i2(map->get_i(move1.id1)), j2(map->get_j(move1.id1));
  Point point(i2,j2), p0(i0,j0), p1(i1,j1);
  int cls = point.classify(p0, p1);
  double dist = fabs((i0 - i1)*j2 + (j1 - j0)*i2 + (j0*i1 - i0*j1))/sqrt(pow(i0 - i1, 2) + pow(j0 - j1, 2));
  double da = (i0 - i2)*(i0 - i2) + (j0 - j2)*(j0 - j2);
  double db = (i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2);
  double ha = sqrt(da - dist*dist);
  double temp=radius*radius - dist*dist;
  double size = temp > CN_EPSILON ? sqrt(temp) : 0;
  if(cls == 3)
  {

    interval.first = move2.t1;
    interval.second = move2.t1 + (sqrt(radius*radius - dist*dist) - ha);
  }
  else if(cls == 4)
  {
    interval.first = move2.t2 - sqrt(radius*radius - dist*dist) + sqrt(db - dist*dist);
    interval.second = move2.t2;
  }
  else if(da < radius*radius - CN_EPSILON)
  {
    if(db < radius*radius - CN_EPSILON)
    {
      interval.first = move2.t1;
      interval.second = move2.t2;
    }
    else
    {
      double hb = sqrt(db - dist*dist);
      interval.first = move2.t1;
      interval.second = move2.t2 - hb + size;
    }
  }
  else
  {
    if(db < radius*radius - CN_EPSILON)
    {
      interval.first = move2.t1 + ha - size;
      interval.second = move2.t2;
      //assert(abs(interval.first-move1.t2)>CN_EPSILON);
    }
    else
    {
      interval.first = move2.t1 + ha - size;
      interval.second = move2.t1 + ha + size;
    }
  }

  interval.second=fmax(interval.second,move1.t1+CN_EPSILON);
  interval.first=fmin(interval.first,move1.t2);
  waitCon=Constraint(agent, interval.first, interval.second, move1.id1, -1,move1.id2);
  assert(waitCon.t2-waitCon.t1>CN_EPSILON);
  constraint.push_back(waitCon);
  if (config.cons_reason)
  {
    if (move2.id2==move1.id1)
    {
      vector<Node> nodes=map->get_valid_moves(move1.id1);
      for (Node n :nodes){
        Vector2D toNode(map->get_coord(move1.id1)),fromNode(map->get_coord(n.id));
        Vector2D temp(toNode-fromNode);
        double dis=sqrt(temp*temp);
        int enter_index=id2ind(move1.id1,n.id,agent);
        assert(enter_index!=-2);
        double min_clear=map->get_min_clear_t_ori(move1.id1,enter_index);
        Move tempMove(-1,-1,n.id,move1.id1);
        double min_clear2=map->get_min_clear_t(tempMove,move2.id2);
        assert(eq(min_clear,min_clear2));
        //double min_clear=map->get_min_clear_t(move1,n.id);
        if (min_clear<-CN_EPSILON) //min_clear==-1
          continue;
        //double startT(interval.second-dis),endT(startT+min_clear);
        double startT(interval.second-dis),endT(move2.t2-dis+min_clear);
        if (gt(startT,0) && gt(min_clear,0))
        {
          assert(gt(endT,startT));
          int exit_id=id2ind(n.id,move1.id1,agent);
          Constraint edgeCon(agent,startT,endT,n.id,exit_id,move1.id1); 
          assert(edgeCon.t2-edgeCon.t1>CN_EPSILON);
          constraint.push_back(edgeCon);
        }
      }
    }
    return make_pair(constraint,true);
  }
  else{
    return make_pair(constraint,false);
  }
  //return Constraint(agent, interval.first-CN_PRECISION, interval.second+CN_PRECISION, move1.id1, -1,move1.id2); // exit index -1 means wait
}

double CBS::get_hl_heuristic(const std::list<Conflict> &conflicts)
{
  if(conflicts.empty() || config.hlh_type == 0)
    return 0;
  else if (config.hlh_type == 1)
  {
    optimization::Simplex simplex("simplex");
    std::map<int, int> colliding_agents;
    for(auto c: conflicts)
    {
      colliding_agents.insert({c.agent1, colliding_agents.size()});
      colliding_agents.insert({c.agent2, colliding_agents.size()});
    }

    pilal::Matrix coefficients(conflicts.size(), colliding_agents.size(), 0);
    std::vector<double> overcosts(conflicts.size());
    int i(0);
    for(auto c:conflicts)
    {
      coefficients.at(i, colliding_agents.at(c.agent1)) = 1;
      coefficients.at(i, colliding_agents.at(c.agent2)) = 1;
      overcosts[i] = c.overcost;
      i++;
    }
    simplex.set_problem(coefficients, overcosts);
    simplex.solve();
    return simplex.get_solution();
  }
  else
  {
    double h_value(0);
    std::vector<std::tuple<double, int, int>> values;
    values.reserve(conflicts.size());
    std::set<int> used;
    for(auto c:conflicts)
      values.push_back(std::make_tuple(c.overcost, c.agent1, c.agent2));
    std::sort(values.begin(), values.end(), std::greater<std::tuple<double, int, int>>());
    for(auto v: values)
    {
      if(used.find(get<1>(v)) != used.end() || used.find(get<2>(v)) != used.end())
        continue;
      h_value += get<0>(v);
      used.insert(get<1>(v));
      used.insert(get<2>(v));
    }
    return h_value;
  }
}

CBS::con_return CBS::get_constraint(int agent, Move move1, Move move2, CBS_Node *node)
{
  Move BKUP1=move1,BKUP2=move2;
  double startT=0,endT=0;

  if(move1.id1 == move1.id2)
  {
    //return make_pair(get_wait_constraint(agent,move1,move2),false);
    CBS::con_return cons=get_wait_constraint(agent, move1, move2);
    return cons;
  }

  move1=split_conf_move(move1,move2,node);
  double startTimeA(move1.t1), endTimeA(move1.t2);
  Vector2D A(map->get_i(move1.id1), map->get_j(move1.id1)), A2(map->get_i(move1.id2), map->get_j(move1.id2)),
           B(map->get_i(move2.id1), map->get_j(move2.id1)), B2(map->get_i(move2.id2), map->get_j(move2.id2));
  if(move2.t2 == CN_INFINITY){
    int exit_index=id2ind(move1.id1, move1.id2,agent);
    list<Constraint> c(0);
    Constraint cons(agent, move1.t1, CN_INFINITY, move1.id1, exit_index, move1.id2);
    c.push_back(cons);
    return make_pair(c,false);
  }

  if (config.cons_reason){
    if (move2.id2<map->get_init_node_num() && move1.id2==move2.id2){
      Vector2D A(map->get_coord(move1.id2)),B(map->get_coord(move1.id1));
      double dist=sqrt((A-B)*(A-B));
      int enter_index=id2ind(move1.id2, move1.id1,agent);
      assert(enter_index!=-2);
      double min_clear=map->get_min_clear_t_ori(move1.id2,enter_index);
      double min_clear2=map->get_min_clear_t(move1,move2.id2);
      cout<<"m1: "<<move1<<"m2: "<<move2<<endl;
      cout<<"min_clear: "<<min_clear<<" min_clear2: "<<min_clear2<<endl;
      assert(eq(min_clear,min_clear2));
      startT=move1.t1;
      if (min_clear>-CN_EPSILON)
      {
        endT=move2.t2-dist+min_clear;
      }
      else
      {
        endT=move2.t2+dist;
      }
    }
  }
  double delta = move2.t2 - move1.t1;
  while(delta > CN_PRECISION/2.0)
    //while(delta > CN_PRECISION)
  {
    if(check_conflict(move1, move2))
    {
      move1.t1 += delta;
      move1.t2 += delta;
    }
    else
    {
      move1.t1 -= delta;
      move1.t2 -= delta;
    }
    if(move1.t1 > move2.t2 + CN_EPSILON)
    {
      move1.t1 = move2.t2;
      move1.t2 = move1.t1 + endTimeA - startTimeA;
      break;
    }
    delta /= 2.0;
  }
  //if(delta < CN_PRECISION/2.0 + CN_PRECISION && check_conflict(move1, move2))
  if(delta < CN_PRECISION/2.0 + CN_EPSILON && check_conflict(move1, move2))
  {
    move1.t1 = fmin(move1.t1 + delta*2, move2.t2);
    move1.t1 = fmax(move1.t1,startTimeA+CN_PRECISION);
    move1.t2 = move1.t1 + endTimeA - startTimeA;
  }
  //consider min_clear time for node conflict
  int exit_index=id2ind(move1.id1, move1.id2,agent);
  if (config.cons_reason){
    /*
       if ( endT<startT || startTimeA-CN_PRECISION>move1.t1+CN_PRECISION)
       {
       check_conflict(BKUP1,BKUP2);
       assert(false);
       }
       */
    if (move1.t1>endT)
    {
      list<Constraint> c(0);
      Constraint cons(agent, startTimeA, move1.t1, move1.id1, exit_index, move1.id2);
      c.push_back(cons);
      return make_pair(c,false);
    }
    else
    {
      list<Constraint> c(0);
      Constraint cons(agent, startT, endT,move1.id1, exit_index, move1.id2);
      c.push_back(cons);
      return make_pair(c,true);
    }
  }
  else{
    list<Constraint> c(0);
    Constraint cons(agent, startTimeA, move1.t1, move1.id1, exit_index, move1.id2);
    c.push_back(cons);
    return make_pair(c,false);
  }
}

int CBS::id2ind(int v1, int v2,int agent)
{
  std::vector<Node> temp_moves=map->get_valid_moves(v1);
  for (int i=0;i<temp_moves.size();++i){
    if (temp_moves[i].id==v2)
      return i;
  }
  printBT_aux();
  assert(false);
  return -2;
}

Conflict CBS::get_conflict(std::list<Conflict> &conflicts)
{
  auto best_it = conflicts.begin();
  for(auto it = conflicts.begin(); it != conflicts.end(); it++)
  {
    if(it->overcost > 0)
    {
      if(best_it->overcost < it->overcost || (fabs(best_it->overcost - it->overcost) < CN_EPSILON && best_it->t < it->t))
        best_it = it;
    }
    else if(best_it->t < it->t)
      best_it = it;
  }

  Conflict conflict = *best_it;
  conflicts.erase(best_it);
  return conflict;
}

Move CBS::find_sub_conflict(Move m1,Move m2,CBS_Node *node)
{//split m1 into a list of moves and return the one is causing conflict with m2
  stack<CBS_Node> infos;
  CBS_Node* curNode = node;
  while(curNode->parent != nullptr){
    infos.push(*curNode);
    curNode = curNode->parent;
  }
  list<int> subM1;
  list<int>::iterator it;

  subM1.push_back(m1.id1);
  subM1.push_back(m1.id2);

  CBS_Node tempNode;
  while (!infos.empty())
  {
    for (auto md:infos.top().deltas)
    {
      //split into a list of moves
      it=find(subM1.begin(),subM1.end(),md.del_edge.first);
      if (it!=subM1.end() )
      {
        if (it!=(--subM1.end()))
        {
          if (*(++it)==md.del_edge.second)
          {
            subM1.insert(it,md.add_node);
          }
          it--;
        }
        if (it!=subM1.begin())
        {
          if (*(--it)==md.del_edge.second)
          {
            it++;
            subM1.insert(it,md.add_node);
          }
        }
      }
    }
    infos.pop();
  }
  double curT(m1.t1);
  for(it=subM1.begin();it!=(--subM1.end());it++)
  {
    int id1((*it++)),id2((*it--));
    double dt(map->get_dist(id1,id2));

    Move temp(curT,curT+dt,id1,id2);
    if (check_conflict(temp,m2) || curT+dt>m2.t2-CN_EPSILON)
    {
      return temp;
    }
    curT+=dt;
  }
  assert(false);
}

Move CBS::split_conf_move(Move m1,Move m2, CBS_Node *node)
{
  if (m1.id1==m1.id2) return m1;

  std::vector<Node> temp_moves=map->get_valid_moves(m1.id1);
  for (int i=0;i<temp_moves.size();++i)
    if (temp_moves[i].id==m1.id2)
      return m1;
  //check conflict for all move1 with move2, return the first one
  Move temp;
  temp=find_sub_conflict(m1,m2,node);
  return temp;
}

Conflict CBS::modify_conflict(Conflict conflict, CBS_Node *node)
{
  Move new_m1,new_m2;
  if (conflict.move1.id1 != conflict.move1.id2)
  {
    bool trigger=true;
    std::vector<Node> temp_moves=map->get_valid_moves(conflict.move1.id1);
    for (int i=0;i<temp_moves.size();++i)
    {
      if (temp_moves[i].id==conflict.move1.id2)
      {
        trigger=false;
        new_m1=conflict.move1;
        break;
      }
    }
    if (trigger)
    {
      new_m1=find_sub_conflict(conflict.move1,conflict.move2,node);
      //check conflict for all move1 with move2, return the first one
    }
  }
  else
  {
    new_m1=conflict.move1;
  }

  if (conflict.move2.id1 != conflict.move2.id2)
  {
    bool trigger=true;
    std::vector<Node> temp_moves=map->get_valid_moves(conflict.move2.id1);
    for (int i=0;i<temp_moves.size();++i)
    {
      if (temp_moves[i].id==conflict.move2.id2)
      {
        trigger=false;
        new_m2=conflict.move2;
        break;
      }
    }
    if (trigger)
    {
      new_m2=find_sub_conflict(conflict.move2,new_m1,node);
      //check conflict for all move1 with move2, return the first one
    }
  }
  else
  {
    new_m2=conflict.move2;
  }

  return Conflict(conflict.agent1, conflict.agent2,new_m1,new_m2);
}

Solution CBS::find_solution(Map &map, const Task &task, const Config &cfg)
{
  config = cfg;
  planner.config=cfg;
  this->map = &map;
  Map temp(map);
  this->original = &temp;
  h_values.init(map.get_size(), task.get_agents_size());
  for(int i = 0; i < int(task.get_agents_size()); i++)
  {
    Agent agent = task.get_agent(i);
    h_values.count(map, agent);
  }
  auto t = std::chrono::high_resolution_clock::now();
  int cardinal_solved = 0, semicardinal_solved = 0;
  if(!this->init_root(map, task))
    return solution;
  solution.init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - t);
  solution.found = true;
  CBS_Node node;
  std::chrono::duration<double> time_spent;
  int expanded(1);
  double time(0);
  std::list<Conflict> conflicts;
  Conflict conflict;
  std::vector<int> conflicting_agents;
  std::vector<std::pair<int, int>> conflicting_pairs;
  int low_level_searches(0);
  int low_level_expanded(0);
  int id = 2;
  bool BREAK=false;
  do
  {
    auto parent = tree.get_front();
    node = *parent;
    node.cost -= node.h;
    parent->conflicts.clear();
    parent->cardinal_conflicts.clear();
    parent->semicard_conflicts.clear();
    //if (node.id==4)
    //  p=false;

    auto paths = get_paths(&node, task.get_agents_size());
    auto time_now = std::chrono::high_resolution_clock::now();
    conflicts = node.conflicts;
    auto cardinal_conflicts = node.cardinal_conflicts;
    auto semicard_conflicts = node.semicard_conflicts;
    if(conflicts.empty() && semicard_conflicts.empty() && cardinal_conflicts.empty())
    {
      break; //i.e. no conflicts => solution found
    }
    if(!cardinal_conflicts.empty())
    {
      conflict = get_conflict(cardinal_conflicts);
      cardinal_solved++;
    }
    else if(!semicard_conflicts.empty())
    {
      conflict = get_conflict(semicard_conflicts);
      semicardinal_solved++;
    }
    else
      conflict = get_conflict(conflicts);

    if (config.use_edge_split)
      gen_new_map(&node);
    check_conflict(conflict.move1,conflict.move2);
    conflict=modify_conflict(conflict,&node);
    parent->cur_conflict=conflict; //del when experiment
                                   //Map_delta_pair info;
    Map_deltas deltasL,deltasR;
    deltasL.clear(); deltasR.clear();
    if(config.use_edge_split)
    {
      ESer.find_deltas(conflict, deltasR, deltasL, map, h_values);
      //split_edge(conflict, paths,deltasR,deltasL);
    }

    time_spent = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - time_now);
    time += time_spent.count();
    expanded++;
    std::list<Constraint> constraintsA = get_constraints(&node, conflict.agent1);
    //std::list<Constraint> constraintsA_New;
    list<Constraint> constraintA;
    pair<list<Constraint>,bool> consBoolPairA;
    if (debug>0)
    {
      cout<<"+"<<endl;
    }

    if (config.use_edge_split && !deltasR.empty() && conflict.move1.id1!=conflict.move1.id2)
    {
      if (debug>1)
      {
        cout<<"move1: ";
      }
      Move temp_move=modify_move(conflict.move1,(deltasR.begin())->add_node);
      if (debug>1)
      {
        cout<<"new_m: ";
      }
      map.alter(deltasR);
      consBoolPairA = get_constraint(conflict.agent1, temp_move, conflict.move2,&node);
      map.alter_back(deltasR);
      constraintA=consBoolPairA.first;  
      if (debug>0)
      {
        cout<<"new constraintA:  ";
      }
      assert(std::find(constraintsA.begin(),constraintsA.end(),*constraintA.begin())==constraintsA.end());
      constraintsA.insert(constraintsA.end(),constraintA.begin(),constraintA.end());
    }
    else
    {
      consBoolPairA=get_constraint(conflict.agent1, conflict.move1, conflict.move2,&node);
      constraintA = consBoolPairA.first; 
      if (debug>0)
      {
        cout<<"new constraintA:  ";
      }
      if (std::find(constraintsA.begin(),constraintsA.end(),*constraintA.begin())!=constraintsA.end()){
        cout<<"breaking B"<<endl;
        BREAK=true;
      }
      constraintsA.insert(constraintsA.end(),constraintA.begin(),constraintA.end());
    }
    if (debug>0)
    {
      cout<<"||"<<endl;
      cout<<"----"<<endl;
    }
    //cout<<"original"<<endl;



    //cout<<"curNode"<<endl;
    sPath pathA;
    if (config.use_edge_split && !deltasR.empty()){
      map.alter(deltasR);
      if (debug>1){
        cout<<"modA"<<endl;
      }
      pathA = planner.find_path(task.get_agent(conflict.agent1), map, constraintsA, h_values);
      map.alter_back(deltasR);
      if (debug>1){
        cout<<"prev: node"<<endl;
      }
    }
    else{
      if(config.debug>1)
        cout<<"start planning"<<endl;
      pathA = planner.find_path(task.get_agent(conflict.agent1), map, constraintsA, h_values);
    }
    if (debug>0){
      cout<<"new_path:";
      cout<<endl;
    }
    low_level_searches++;
    low_level_expanded += pathA.expanded;
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if (debug>0)
      cout<<"-----------------------------------------"<<endl;
    std::list<Constraint> constraintsB = get_constraints(&node, conflict.agent2);
    list<Constraint> constraintB;
    pair<list<Constraint>,bool> consBoolPairB;
    if (debug>0){
      cout<<"+"<<endl;
    }
    if (config.use_edge_split && !deltasL.empty() && conflict.move2.id1!=conflict.move2.id2){
      if (debug>1){
        cout<<"move2: ";
      }
      Move temp_move=modify_move(conflict.move2,(deltasL.begin())->add_node);
      if (debug>1){
        cout<<"new_m: ";
      }
      map.alter(deltasL);
      consBoolPairB= get_constraint(conflict.agent2, temp_move, conflict.move1,&node);
      map.alter_back(deltasL);
      constraintB= consBoolPairB.first;
      if (debug>0)
      {
        cout<<"new constraintB:  ";
      }
      assert(std::find(constraintsB.begin(),constraintsB.end(),*constraintB.begin())==constraintsB.end());
      constraintsB.insert(constraintsB.end(),constraintB.begin(),constraintB.end());
    }
    else{
      consBoolPairB = get_constraint(conflict.agent2, conflict.move2, conflict.move1,&node);
      constraintB= consBoolPairB.first;
      if (debug>0)
      {
        cout<<"new constraintB:  ";
      }
      if (std::find(constraintsB.begin(),constraintsB.end(),*constraintB.begin())!=constraintsB.end()){
        cout<<"breaking B"<<endl;
        BREAK=true;
      }
      if (debug>1){
        cout<<"moveB1:";
        cout<<", moveB2: ";
      }
      constraintsB.insert(constraintsB.end(),constraintB.begin(),constraintB.end());
    }
    if (debug>0){
      cout<<"||"<<endl;
      cout<<"----"<<endl;
      cout<<endl;
      cout<<flush;
    }
    sPath pathB;

    if (config.use_edge_split && !deltasL.empty()){
      map.alter(deltasL);
      if (debug>1){
        cout<<"modB"<<endl;
      }
      pathB = planner.find_path(task.get_agent(conflict.agent2), map, constraintsB, h_values);
      map.alter_back(deltasL);
      if (debug>2){
        cout<<"prev: node"<<endl;
      }
    }
    else{
      pathB = planner.find_path(task.get_agent(conflict.agent2), map, constraintsB, h_values);
    }
    if (debug>0){
      cout<<"new_path:";
      cout<<endl<<flush;
      //assert(false);
    }
      assert(!BREAK); 

    low_level_searches++;
    low_level_expanded += pathB.expanded;

    CBS_Node right({pathA}, parent, constraintA,deltasR, node.cost + pathA.cost - get_cost(node, conflict.agent1), 0, node.total_cons + 1);
    CBS_Node left ({pathB}, parent, constraintB,deltasL, node.cost + pathB.cost - get_cost(node, conflict.agent2), 0, node.total_cons + 1);
    Constraint positive;
    bool inserted = false;
    bool left_ok = true, right_ok = true;
    if(config.use_disjoint_splitting)
    {
      Constraint tempConsA(*constraintA.begin()),tempConsB(*constraintB.begin());
      int agent1positives(0), agent2positives(0);
      for(auto c: constraintsA)
        if(c.positive)
          agent1positives++;
      for(auto c: constraintsB)
        if(c.positive)
          agent2positives++;

      /*
         Constraint cons2changeA,cons2changeB;
         if (config.use_edge_split && info.first.add_node!=-1)
         cons2changeA=tempA;
         else
         cons2changeA=constraintA;

         if (config.use_edge_split && info.second.add_node!=-1)
         cons2changeB=tempB;
         else
         cons2changeB=constraintB;
         */

      if(conflict.move1.id1 != conflict.move1.id2 && agent2positives > agent1positives && pathA.cost > 0)
      //if(!consBoolPairB.second && conflict.move1.id1 != conflict.move1.id2 && agent2positives > agent1positives && pathA.cost > 0)
      {
        int exit_index=id2ind(conflict.move1.id1, conflict.move1.id2,conflict.agent1);
        positive = Constraint(conflict.agent1, tempConsA.t1, tempConsA.t2, tempConsA.id1, exit_index, tempConsA.to_id, true);
        //positive = Constraint(conflict.agent1, constraintA.t1, constraintA.t2, conflict.move1.id1, exit_index, conflict.move1.id2, true);
        if(check_positive_constraints(constraintsA, positive))
        {
          left.positive_constraint = positive;
          left.total_cons++;
          constraintsB.push_back(left.positive_constraint);
          if (debug>0)
          {
            cout<<"new positive_constraint:"<<endl;
          }
          inserted = true;
        }
      }
      if(conflict.move2.id1 != conflict.move2.id2 && !inserted && pathB.cost > 0)
      //if(!consBoolPairA.second && conflict.move2.id1 != conflict.move2.id2 && !inserted && pathB.cost > 0)
      {
        int exit_index=id2ind(conflict.move2.id1, conflict.move2.id2,conflict.agent2);
        positive = Constraint(conflict.agent2, tempConsB.t1, tempConsB.t2, conflict.move2.id1, exit_index, conflict.move2.id2, true);
        if(check_positive_constraints(constraintsB, positive))
        {
          right.positive_constraint = positive;
          right.total_cons++;
          constraintsA.push_back(right.positive_constraint);
          inserted = true;
          if (debug>0)
          {
            cout<<"new positive_constraint:"<<endl;
          }
        }
      }
      if(conflict.move1.id1 != conflict.move1.id2 && !inserted && pathA.cost > 0)
      //if(!consBoolPairA.second && conflict.move1.id1 != conflict.move1.id2 && !inserted && pathA.cost > 0)
      {
        int exit_index=id2ind(conflict.move1.id1, conflict.move1.id2,conflict.agent1);
        positive = Constraint(conflict.agent1, tempConsA.t1, tempConsA.t2, conflict.move1.id1, exit_index, conflict.move1.id2, true);
        if(check_positive_constraints(constraintsA, positive))
        {
          inserted = true;
          left.positive_constraint = positive;
          left.total_cons++;
          constraintsB.push_back(left.positive_constraint);
        }
      }
    }
    if (config.use_edge_split){
      gen_original_map(&node);
      assert(map.equal(original));
    }
    right.id_str = node.id_str + "0";
    left.id_str = node.id_str + "1";
    right.id = id++;
    left.id = id++;

    if(right_ok && pathA.cost > 0 && validate_constraints(constraintsA, pathA.agentID))
    {
      time_now = std::chrono::high_resolution_clock::now();
      if (config.use_edge_split)
        gen_new_map(&right);

      find_new_conflicts(map, task, right, paths, pathA, conflicts, semicard_conflicts, cardinal_conflicts, low_level_searches, low_level_expanded);
      if (config.use_edge_split)
        gen_original_map(&right);
      time_spent = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - time_now);
      time += time_spent.count();
      if(right.cost > 0)
      {
        right.h = get_hl_heuristic(right.cardinal_conflicts);
        right.cost += right.h;
        tree.add_node(right);
        /*
           if (node.cost+node.h-0.1>right.cost+CN_PRECISION){
           string file="CT_tree_no_sol.dot";
           saveCT(file,&right,task.get_agents_size());
           assert(false);
           }
           */
      }
    }
    if(left_ok && pathB.cost > 0 && validate_constraints(constraintsB, pathB.agentID))
    {
      time_now = std::chrono::high_resolution_clock::now();
      if (config.use_edge_split)
        gen_new_map(&left);
      find_new_conflicts(map, task, left, paths, pathB, conflicts, semicard_conflicts, cardinal_conflicts, low_level_searches, low_level_expanded);
      if (config.use_edge_split)
        gen_original_map(&left);
      time_spent = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - time_now);
      time += time_spent.count();           
      if(left.cost > 0)
      {
        left.h = get_hl_heuristic(left.cardinal_conflicts);
        left.cost += left.h;
        tree.add_node(left);

      }
    }
    time_spent = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - t);
    if(time_spent.count() > config.timelimit )
    {
      solution.found = false;
      if (debug>0){
        printBT_aux();
      }
      break;
    }
  }
  while(tree.get_open_size() > 0);
  solution.paths = get_paths(&node, task.get_agents_size());
  solution.flowtime = node.cost;
  solution.low_level_expansions = low_level_searches;
  solution.low_level_expanded = double(low_level_expanded)/std::max(low_level_searches, 1);
  solution.high_level_expanded = expanded;
  solution.high_level_generated = int(tree.get_size());
  for(auto path:solution.paths)
    solution.makespan = (solution.makespan > path.cost) ? solution.makespan : path.cost;
  solution.time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - t);
  solution.check_time = time;
  solution.cardinal_solved = cardinal_solved;
  solution.semicardinal_solved = semicardinal_solved;
  solution.new_node = map.get_new_node_num();

  if (solution.found)
  {
    post_check(solution.paths);
  }
  return solution;
}

bool CBS::check_positive_constraints(std::list<Constraint> constraints, Constraint constraint)
{
  std::list<Constraint> positives;
  for(auto c: constraints)
    if(c.positive && c.agent == constraint.agent)
      positives.push_back(c);

  for(auto p: positives)
  {
    if(p.id1 == constraint.id1 && p.id2 == constraint.id2 && p.t1 - CN_EPSILON < constraint.t1 && p.t2 + CN_EPSILON > constraint.t2) // agent needs to perform two equal actions simultaneously => it's impossible
      return false;
    if(p.id1 == constraint.id1 && p.id2 == constraint.id2 && constraint.t1 - CN_EPSILON < p.t1 && constraint.t2 + CN_EPSILON > p.t2)
      return false;
  }
  return true;
}

bool CBS::validate_constraints(std::list<Constraint> constraints, int agent)
{
  std::list<Constraint> positives;
  for(auto c: constraints)
    if(c.positive && c.agent == agent)
      positives.push_back(c);
  if (config.debug>1)
  {
    cout<<"validate_constraints:"<<endl;
    cout<<"---start checking--"<<endl;
  }
  for(auto p: positives)
    for(auto c: constraints)
    {
      if(c.positive)
        continue;
      if (config.debug>1)
      {
        cout<<"---"<<endl;
      }
      if(p.agent == c.agent && p.id1 == c.id1 && p.id2 == c.id2) //if the same action
      {
        if(config.debug>1)
        {
          cout<<"inner checking p.t1:"<<p.t1<<" c.t1:"<<c.t1<<" and p.t2:"<<p.t2<<" c.t2"<<c.t2<<endl;
          cout<<": ="<< (p.t1 - c.t1) << "and "<< (p.t2 - c.t2)<<endl; 
          cout<<": ="<< (p.t1  > c.t1 - CN_EPSILON) << "and "<< (p.t2 < c.t2 + CN_EPSILON)<<endl; 
        }
        if(p.t1 > c.t1 - CN_EPSILON && p.t2 < c.t2 + CN_EPSILON) //if the whole positive interval is inside collision interval
          return false;
      }
      if(config.debug>1)
      {
        cout<<"pass"<<endl;
        cout<<"---"<<endl;
      }
    }
  return true;
}

void CBS::find_new_conflicts(Map &map, const Task &task, CBS_Node &node, std::vector<sPath> &paths, const sPath &path,
    const std::list<Conflict> &conflicts, const std::list<Conflict> &semicard_conflicts, const std::list<Conflict> &cardinal_conflicts,
    int &low_level_searches, int &low_level_expanded)
{
  auto oldpath = paths[path.agentID];
  paths[path.agentID] = path;
  auto new_conflicts = get_all_conflicts(paths, path.agentID);
  paths[path.agentID] = oldpath;
  std::list<Conflict> conflictsA({}), semicard_conflictsA({}), cardinal_conflictsA({});
  for(auto c: conflicts)
    if(c.agent1 != path.agentID && c.agent2 != path.agentID)
      conflictsA.push_back(c);
  for(auto c: semicard_conflicts)
    if(c.agent1 != path.agentID && c.agent2 != path.agentID)
      semicard_conflictsA.push_back(c);
  for(auto c: cardinal_conflicts)
    if(c.agent1 != path.agentID && c.agent2 != path.agentID)
      cardinal_conflictsA.push_back(c);
  if(!config.use_cardinal)
  {
    node.conflicts = conflictsA;
    for(auto n:new_conflicts)
      node.conflicts.push_back(n);
    node.cardinal_conflicts.clear();
    node.semicard_conflicts.clear();
    node.conflicts_num = node.conflicts.size();
    return;
  }
  for(auto c: new_conflicts)
  {
    std::list<Constraint> constraintsA, constraintsB;
    if(path.agentID == c.agent1)
    {
      constraintsA = get_constraints(&node, c.agent1);
      con_return temp=get_constraint(c.agent1, c.move1, c.move2, &node);
      constraintsA.insert(constraintsA.end(),temp.first.begin(),temp.first.end());
      auto new_pathA = planner.find_path(task.get_agent(c.agent1), map, constraintsA, h_values);
      constraintsB = get_constraints(&node, c.agent2);
      temp=get_constraint(c.agent2, c.move2, c.move1,&node);
      constraintsB.insert(constraintsB.end(),temp.first.begin(),temp.first.end());
      auto new_pathB = planner.find_path(task.get_agent(c.agent2), map, constraintsB, h_values);
      double old_cost = get_cost(node, c.agent2);
      if(new_pathA.cost < 0 && new_pathB.cost < 0)
      {
        node.cost = -1;
        return;
      }
      else if (new_pathA.cost < 0)
      {
        c.overcost = new_pathB.cost - old_cost;
        cardinal_conflictsA.push_back(c);
      }
      else if (new_pathB.cost < 0)
      {
        c.overcost = new_pathA.cost - path.cost;
        cardinal_conflictsA.push_back(c);
      }
      else if(new_pathA.cost > path.cost && new_pathB.cost > old_cost)
      {
        c.overcost = std::min(new_pathA.cost - path.cost, new_pathB.cost - old_cost);
        cardinal_conflictsA.push_back(c);
      }
      else if(new_pathA.cost > path.cost || new_pathB.cost > old_cost)
        semicard_conflictsA.push_back(c);
      else
        conflictsA.push_back(c);
      low_level_searches += 2;
      low_level_expanded += (new_pathA.expanded + new_pathB.expanded);
    }
    else
    {
      constraintsA = get_constraints(&node, c.agent2);
      con_return temp=get_constraint(c.agent2, c.move2, c.move1,&node);
      constraintsA.insert(constraintsA.end(),temp.first.begin(),temp.first.end());
      auto new_pathA = planner.find_path(task.get_agent(c.agent2), map, constraintsA, h_values);
      constraintsB = get_constraints(&node, c.agent1);
      temp=get_constraint(c.agent1, c.move1, c.move2,&node);
      constraintsB.insert(constraintsB.end(),temp.first.begin(),temp.first.end());
      auto new_pathB = planner.find_path(task.get_agent(c.agent1), map, constraintsB, h_values);
      double old_cost = get_cost(node, c.agent1);
      if(new_pathA.cost < 0 && new_pathB.cost < 0)
      {
        node.cost = -1;
        return;
      }
      else if (new_pathA.cost < 0)
      {
        c.overcost = new_pathB.cost - old_cost;
        cardinal_conflictsA.push_back(c);
      }
      else if (new_pathB.cost < 0)
      {
        c.overcost = new_pathA.cost - path.cost;
        cardinal_conflictsA.push_back(c);
      }
      else if(new_pathA.cost > path.cost && new_pathB.cost > old_cost)
      {
        c.overcost = std::min(new_pathA.cost - path.cost, new_pathB.cost - old_cost);
        cardinal_conflictsA.push_back(c);
      }
      else if(new_pathA.cost > path.cost || new_pathB.cost > old_cost)
        semicard_conflictsA.push_back(c);
      else
        conflictsA.push_back(c);
      low_level_searches += 2;
      low_level_expanded += (new_pathA.expanded + new_pathB.expanded);
    }
  }

  node.conflicts = conflictsA;
  node.semicard_conflicts = semicard_conflictsA;
  node.cardinal_conflicts = cardinal_conflictsA;
  node.conflicts_num = conflictsA.size() + semicard_conflictsA.size() + cardinal_conflictsA.size();
  return;
}

std::list<Constraint> CBS::get_constraints(CBS_Node *node, int agent_id)
{
  CBS_Node* curNode = node;
  std::list<Constraint> constraints(0);
  while(curNode->parent != nullptr)
  {

    for (Constraint c:curNode->constraint){
    	if(agent_id < 0 || c.agent == agent_id)
    		constraints.push_back(c);
    }
    //if(agent_id < 0 || curNode->constraint.agent == agent_id)
    //  constraints.push_back(curNode->constraint);

    if(curNode->positive_constraint.agent == agent_id)
      constraints.push_back(curNode->positive_constraint);
    curNode = curNode->parent;
  }
  return constraints;
}

Conflict CBS::check_paths(const sPath &pathA, const sPath &pathB)
{
  unsigned int a(0), b(0);
  auto nodesA = pathA.nodes;
  auto nodesB = pathB.nodes;
  while(a < nodesA.size() - 1 || b < nodesB.size() - 1)
  {
    //cout<<"a: "<<a<<" b: "<<b<<endl<<flush;
    double temp=pow(map->get_i(nodesA[a].id) - map->get_i(nodesB[b].id), 2) + pow(map->get_j(nodesA[a].id) - map->get_j(nodesB[b].id), 2);
    if (abs(temp)<CN_EPSILON) temp=0;
    assert(temp>=0);
    double dist = sqrt(temp); 
    if(a < nodesA.size() - 1 && b < nodesB.size() - 1) // if both agents have not reached their goals yet
    {
      if(dist < (nodesA[a+1].g - nodesA[a].g) + (nodesB[b+1].g - nodesB[b].g) + config.agent_size*2 - CN_EPSILON)
        if(check_conflict(Move(nodesA[a], nodesA[a+1]), Move(nodesB[b], nodesB[b+1]))){
          return Conflict(pathA.agentID, pathB.agentID, Move(nodesA[a], nodesA[a+1]), Move(nodesB[b], nodesB[b+1]), std::min(nodesA[a].g, nodesB[b].g));
        }
    }
    else if(a == nodesA.size() - 1) // if agent A has already reached the goal
    {
      if(dist < (nodesB[b+1].g - nodesB[b].g) + config.agent_size*2-CN_EPSILON)
        if(check_conflict(Move(nodesA[a].g, CN_INFINITY, nodesA[a].id, nodesA[a].id), Move(nodesB[b], nodesB[b+1]))){
          return Conflict(pathA.agentID, pathB.agentID, Move(nodesA[a].g, CN_INFINITY, nodesA[a].id, nodesA[a].id), Move(nodesB[b], nodesB[b+1]), std::min(nodesA[a].g, nodesB[b].g));
        }
    }
    else if(b == nodesB.size() - 1) // if agent B has already reached the goal
    {
      if(dist < (nodesA[a+1].g - nodesA[a].g) + config.agent_size*2 -CN_EPSILON)
        if(check_conflict(Move(nodesA[a], nodesA[a+1]), Move(nodesB[b].g, CN_INFINITY, nodesB[b].id, nodesB[b].id))){
          return Conflict(pathA.agentID, pathB.agentID, Move(nodesA[a], nodesA[a+1]), Move(nodesB[b].g, CN_INFINITY, nodesB[b].id, nodesB[b].id), std::min(nodesA[a].g, nodesB[b].g));
        }
    }
    //if (a==3 && b==3)
    //{
      //cout<<"A cost: "<<nodesA[a+1].g<<" B cost: "<<nodesB[b+1].g<<endl;
      //cout<<"     |a-b|"<< fabs(nodesA[a+1].g - nodesB[b+1].g)<<endl;
      //cout<<"CN_EPSILON"<<CN_EPSILON<<endl;
      //cout<<"opt1: "<<(nodesA[a+1].g >= nodesB[b+1].g - CN_EPSILON)<<endl;
      //cout<<"opt2: "<<(nodesA[a+1].g <= nodesB[b+1].g + CN_EPSILON)<<endl;
      //cout<<"opt3a: "<<((nodesA[a+1].g  < nodesB[b+1].g + CN_EPSILON)) << endl;
      //cout<<"opt3c: "<<(nodesA[a+1].g == nodesB[b+1].g + CN_EPSILON)<<endl;
      //cout<<"opt3b: "<<(nodesA[a+1].g > nodesB[b+1].g-CN_EPSILON) << endl;
      //cout<<"opt4a: "<<(nodesA[a+1].g <= nodesB[b+1].g + CN_EPSILON && nodesA[a+1].g >= nodesB[b+1].g - CN_EPSILON) << endl;
      //cout<<"opt4: "<<(abs(nodesA[a+1].g - nodesB[b+1].g) <= CN_EPSILON)<<endl;
      //cout<<"opt5: "<<(nodesA[a+1].g - nodesB[b+1].g < - CN_EPSILON)<<endl;
      //cout<<"opt6: "<<(nodesB[b+1].g  - nodesA[a+1].g < - CN_EPSILON)<<endl<<flush;
      //assert(false);
    //}
    if(a == nodesA.size() - 1)
    {
      b++;
    }
    else if(b == nodesB.size() - 1)
    {
      a++;
    }
    else if(eq(nodesA[a+1].g, nodesB[b+1].g))
    {
      a++;
      b++;
    }
    else if(lt(nodesA[a+1].g,nodesB[b+1].g))
    {
      a++;
    }
    else if(gt(nodesA[a+1].g,nodesB[b+1].g))
    {
      b++;
    }
  }
  return Conflict();
}

Conflict CBS::check_paths_ori(const sPath &pathA, const sPath &pathB)
{
  unsigned int a(0), b(0);
  auto nodesA = pathA.nodes;
  auto nodesB = pathB.nodes;
  while(a < nodesA.size() - 1 || b < nodesB.size() - 1)
  {
    double temp=pow(map->get_i(nodesA[a].id) - map->get_i(nodesB[b].id), 2) + pow(map->get_j(nodesA[a].id) - map->get_j(nodesB[b].id), 2);
    if (abs(temp)<CN_EPSILON) temp=0;
    assert(temp>=0);
    double dist = sqrt(temp); 
    if(a < nodesA.size() - 1 && b < nodesB.size() - 1) // if both agents have not reached their goals yet
    {
      if(dist < (nodesA[a+1].g - nodesA[a].g) + (nodesB[b+1].g - nodesB[b].g) + config.agent_size*2 - CN_EPSILON)
        if(check_conflict_ori(Move(nodesA[a], nodesA[a+1]), Move(nodesB[b], nodesB[b+1]))){
          return Conflict(pathA.agentID, pathB.agentID, Move(nodesA[a], nodesA[a+1]), Move(nodesB[b], nodesB[b+1]), std::min(nodesA[a].g, nodesB[b].g));
        }
    }
    else if(a == nodesA.size() - 1) // if agent A has already reached the goal
    {
      if(dist < (nodesB[b+1].g - nodesB[b].g) + config.agent_size*2-CN_EPSILON)
        if(check_conflict_ori(Move(nodesA[a].g, CN_INFINITY, nodesA[a].id, nodesA[a].id), Move(nodesB[b], nodesB[b+1]))){
          return Conflict(pathA.agentID, pathB.agentID, Move(nodesA[a].g, CN_INFINITY, nodesA[a].id, nodesA[a].id), Move(nodesB[b], nodesB[b+1]), std::min(nodesA[a].g, nodesB[b].g));
        }
    }
    else if(b == nodesB.size() - 1) // if agent B has already reached the goal
    {
      if(dist < (nodesA[a+1].g - nodesA[a].g) + config.agent_size*2 -CN_EPSILON)
        if(check_conflict_ori(Move(nodesA[a], nodesA[a+1]), Move(nodesB[b].g, CN_INFINITY, nodesB[b].id, nodesB[b].id))){
          return Conflict(pathA.agentID, pathB.agentID, Move(nodesA[a], nodesA[a+1]), Move(nodesB[b].g, CN_INFINITY, nodesB[b].id, nodesB[b].id), std::min(nodesA[a].g, nodesB[b].g));
        }
    }
    if(a == nodesA.size() - 1)
      b++;
    else if(b == nodesB.size() - 1)
      a++;
    else if(fabs(nodesA[a+1].g - nodesB[b+1].g) < CN_EPSILON)
    {
      a++;
      b++;
    }
    else if(nodesA[a+1].g < nodesB[b+1].g - CN_EPSILON)
      a++;
    else if(nodesB[b+1].g  < nodesA[a+1].g + CN_EPSILON)
      b++;
  }
  return Conflict();
}
std::vector<Conflict> CBS::get_all_conflicts(const std::vector<sPath> &paths, int id)
{
  std::vector<Conflict> conflicts;
  //check all agents
  if(id < 0)
    for(unsigned int i = 0; i < paths.size(); i++)
      for(unsigned int j = i + 1; j < paths.size(); j++)
      {
        Conflict conflict = check_paths(paths[i], paths[j]);
        if(conflict.agent1 >= 0)
          conflicts.push_back(conflict);
      }
  else
  {
    for(unsigned int i = 0; i < paths.size(); i++)
    {
      if(int(i) == id)
        continue;
      Conflict conflict = check_paths(paths[i], paths[id]);
      if(conflict.agent1 >= 0)
        conflicts.push_back(conflict);
    }
  }
  return conflicts;
}

double CBS::get_cost(CBS_Node node, int agent_id)
{
  while(node.parent != nullptr)
  {
    if(node.paths.begin()->agentID == agent_id)
      return node.paths.begin()->cost;
    node = *node.parent;
  }
  return node.paths.at(agent_id).cost;
}

std::vector<sPath> CBS::get_paths(CBS_Node *node, unsigned int agents_size)
{
  CBS_Node* curNode = node;
  std::vector<sPath> paths(agents_size);
  while(curNode->parent != nullptr)
  {
    if(paths.at(curNode->paths.begin()->agentID).cost < 0)
      paths.at(curNode->paths.begin()->agentID) = *curNode->paths.begin();
    curNode = curNode->parent;
  }
  for(unsigned int i = 0; i < agents_size; i++)
    if(paths.at(i).cost < 0)
      paths.at(i) = curNode->paths.at(i);
  return paths;
}






Vector2D CBS::ind2Vec(int nodeId)
{
  return Vector2D(map->get_i(nodeId),map->get_j(nodeId));
}




