#include "sipp.h"
void SIPP::clear()
{
    open.clear();
    close.clear();
    collision_intervals.clear();
    landmarks.clear();
    constraints.clear();
    visited.clear();
    path.cost = -1;
}

double SIPP::dist(const Node& a, const Node& b)
{
    return std::sqrt(pow(a.i - b.i, 2) + pow(a.j - b.j, 2));
}

void SIPP::find_successors(Node curNode, const Map &map, std::list<Node> &succs, Heuristic &h_values, Node goal)
{
    Node newNode;
    std::vector<Node> valid_moves = map.get_valid_moves(curNode.id);
    if(config.debug>1 || p)
    {
        std::cout<<"expanding: "<<curNode.id<<" g:"<<curNode.g<<"interval:("<<curNode.interval.first<<", "<<curNode.interval.second<<") "<<std::endl;
        for (Node n:valid_moves){
            std::cout<<"->("<<n.id<<", dist: "<<dist(curNode,n)<<")";
        }
        std::cout<<std::endl;
    }
    for(int i=0;i<valid_moves.size();++i)
    {
        Node move=valid_moves[i];
        if (config.debug>1 || p)
        {
            cout<<"-------------------------"<<endl;
            cout<<"checking: "<<move.id<<endl;
        }
        newNode.i = move.i;
        newNode.j = move.j;
        newNode.id = move.id;
        double cost = dist(curNode, newNode);
        vector<Interval> curIntervals(0);
        //cout<<"cons idx: "<<curNode.id<<" ~ "<<i<<endl;
        auto moveColls_it = constraints.find({curNode.id,i});
        if(moveColls_it != constraints.end())
        {
            Interval interval = curNode.interval;
            Interval moveCons;
            for(unsigned int i = 0; i < moveColls_it->second.size(); i++)
            {
                moveCons=moveColls_it->second[i];
                //cout<<"moveCons: "<<moveCons.first<<","<<moveCons.second<<endl;
                if(gt_raw(moveCons.first,curNode.interval.second))
                {
                    curIntervals.push_back(interval);
                    break;
                }
                if (le_raw(moveCons.first,interval.first))
                {
                    if(ge_raw(moveCons.second,curNode.interval.second))
                    {
                        //interval.second=curNode.interval.second;
                        //if()
                        //curIntervals.push_back(interval);
                        break;
                    }
                    //if(lt_raw(moveCons.second,interval.second))
                    if(ge_raw(moveCons.second,interval.first))
                        interval.first=moveCons.second;
                    //curIntervals.push_back(interval);
                }
                else {
                    interval.second=moveCons.first;
                    curIntervals.push_back(interval);
                    //cout<<interval.first<<", "<<interval.second<<endl;
                    if(ge_raw(moveCons.second,curNode.interval.second))
                        break;
                    interval.first=moveCons.second;
                    interval.second=curNode.interval.second;
                }
            }
            if(lt_raw(moveCons.second,interval.second))
            {
                interval.second=curNode.interval.second;
                curIntervals.push_back(interval);
            }
        }
        else
            curIntervals.push_back(curNode.interval);
        assert(curNode.interval.first!=-1);
        if(config.debug>1 || p)
        {
            cout<<"cur interval"<<endl;
            for (auto i:curIntervals)
                cout<<"("<<i.first<<","<<i.second<<")"<<endl;
        }

        std::vector<Interval> nextIntervals(0);
        auto waitColls_it = collision_intervals.find(newNode.id);
        if(waitColls_it != collision_intervals.end())
        {
            std::pair<double, double> interval = {0, CN_INFINITY};
            for(unsigned int i = 0; i < waitColls_it->second.size(); i++)
            {
                interval.second = waitColls_it->second[i].first;
                nextIntervals.push_back(interval);
                interval.first = waitColls_it->second[i].second;
            }
            interval.second = CN_INFINITY;
            nextIntervals.push_back(interval);
        }
        else
            nextIntervals.push_back({0, CN_INFINITY});

        if(config.debug>1 || p)
        {
            cout<<"next interval"<<endl;
            for (auto i:nextIntervals)
                cout<<"("<<i.first<<","<<i.second<<")"<<endl;
        }

        int id(0);
        /*
        auto itcur intervalcur interval=nextIntervals.begin();
        while(lt_raw(it->second,curIntervals[0].first+cost))
        {
            it++;
            id++;
        }
        */
        for(Interval next:nextIntervals)
        {
            newNode.interval_id=id;
            newNode.interval=next;
            id++;
            //for (;it!=nextIntervals.end();++it)
            for (Interval cur:curIntervals)
            {
                //cout<<"cur: "<<cur.first<<" ~ "<<cur.second<<endl;
                //cout<<"next: "<<next.first<<" ~ "<<next.second<<endl;
                //cout<<"cur+csot: "<<cur.first+cost<<" ~ "<<cur.second+cost<<endl;
                //cout<<"cur+cost: "<<cur.first+cost<<"   next-cost: "<<next.first-cost<<endl;
            //while (lt_raw(i.first+cost,it->second)  && gt_raw(i.second+cost,it->first)) // iterate all reachable next interval
                if(!(lt_raw(cur.first+cost,next.second) && gt_raw(cur.second,next.first-cost)))
                    continue;
                //cout<<"over lapping: "<<endl;
                auto vis_it = visited.find(make_tuple(move.id,newNode.interval_id,false));
                /*
                if (vis_it!=visited.end())
                    if(vis_it->second.second)
                        continue;
                */
                if(ge_raw(cur.first+cost,next.first))
                {
                    newNode.g=cur.first+cost;
                    newNode.prev_g=cur.first;
                    //cout<<"updating 1: "<<newNode.prev_g<<"->"<<newNode.g<<endl;
                }
                else {
                    newNode.g=next.first;
                    newNode.prev_g=next.first-cost;
                    //cout<<"updating 2: "<<newNode.prev_g<<"->"<<newNode.g<<endl;
                }
                newNode.interval.first=newNode.g;

                if(vis_it != visited.end())
                {
                    if(le_raw(vis_it->second.first, newNode.g)) // <=
                        continue;
                    else
                    {
                        vis_it->second.first = newNode.g;
                        vis_it->second.second=false;
                    }
                }
                else
                    visited.insert({make_tuple(newNode.id,newNode.interval_id,false), {newNode.g, false}});
                if(goal.id == agent.goal_id) //perfect heuristic is known
                    newNode.f = newNode.g + h_values.get_value(newNode.id, agent.id);
                else
                {
                    double h = sqrt(pow(goal.i - newNode.i, 2) + pow(goal.j - newNode.j, 2));
                    for(unsigned int i = 0; i < h_values.get_size(); i++) //differential heuristic with pivots placed to agents goals
                        h = std::max(h, fabs(h_values.get_value(newNode.id, i) - h_values.get_value(goal.id, i)));
                    newNode.f = newNode.g + h;
                }
                succs.push_back(newNode);
                break;
                //it++;
                //if (it==nextIntervals.end())
                //    break;
            }
            //if (it!=nextIntervals.end())
            //    break;
        }
        if(config.debug>1 || p)
        {
            for (Node n:succs){
                std::cout<<"("<<n.id<<"@"<<n.g<<")";
            }
            std::cout<<"]"<<std::endl;
        }
    }
    if(config.debug>1 || p)
        cout<<"----------"<<endl;
}

Node SIPP::find_min()
{
    Node min = *open.begin();
    open.pop_front();
    return min;
}

void SIPP::add_open(Node newNode)
{
    if (open.empty() || le_raw(open.back().f, newNode.f))
    {
        open.push_back(newNode);
        return;
    }
    for(auto iter = open.begin(); iter != open.end(); ++iter)
    {
        if(gt_raw(iter->f, newNode.f)) // if newNode.f has lower f-value  >
        {
            open.insert(iter, newNode);
            return;
        }
        else if(eq_raw(iter->f,newNode.f) && ge_raw(newNode.g, iter->g)) // if f-values are equal, compare g-values  == && >=
        {
            open.insert(iter, newNode);
            return;
        }
    }
    open.push_back(newNode);
    return;
}

std::vector<Node> SIPP::reconstruct_path(Node curNode)
{
    path.nodes.clear();
    if(curNode.parent != nullptr)
        do
        {
            if (config.debug>1 || p)
                cout<<curNode.parent->id<<"@"<<curNode.parent->g<<"->"<<(curNode.prev_g)<<"->"<<curNode.id<<"@"<<curNode.g<<endl;
            assert(curNode.prev_g!=-1);
            assert(curNode.prev_g>=curNode.parent->g);
            path.nodes.insert(path.nodes.begin(), curNode);
            if(!eq_raw(curNode.prev_g,curNode.parent->g))
            {
                Node add = *curNode.parent;
                add.g = curNode.prev_g;
                path.nodes.insert(path.nodes.begin(), add);
                if (config.debug>1 || p)
                    cout<<"insert a wait: "<<path.nodes[0].g<<endl;
            }
            curNode = *curNode.parent;
        }
        while(curNode.parent != nullptr);
    path.nodes.insert(path.nodes.begin(), curNode);
    /*if(!eq_raw(curNode.g,0))
    {
        Node add = curNode;
        add.g = 0;
        path.nodes.insert(path.nodes.begin(), add);
    }*/
    if(config.debug>1 || p)
    {
        cout<<"recons nodes: "<<endl;
        prt_nodes(path.nodes);
    }
    /*
       for(unsigned int i = 0; i < path.nodes.size(); i++)
       {
       unsigned int j = i + 1;
       if(j == path.nodes.size())
       break;
       cout<<"path dif: "<<path.nodes[j].g - path.nodes[i].g <<endl;
       cout<<"dist:     "<<dist(path.nodes[j], path.nodes[i])<<endl;
       if(!eq_raw(path.nodes[j].g - path.nodes[i].g, dist(path.nodes[j], path.nodes[i]))) //==
       {
       Node add = path.nodes[i];
       add.g = path.nodes[j].g - dist(path.nodes[j], path.nodes[i]);
       path.nodes.emplace(path.nodes.begin() + j, add);
       }
       cout<<"stuck here"<<endl;
       }
       */
    return path.nodes;
}

void SIPP::add_collision_interval(int id, std::pair<double, double> interval)
{
    std::vector<std::pair<double, double>> nextIntervals(0);
    if(collision_intervals.count(id) == 0)
        collision_intervals.insert({id, {interval}});
    else
        collision_intervals[id].push_back(interval);
    std::sort(collision_intervals[id].begin(), collision_intervals[id].end());
    if(config.debug>1 || p)
    {
        cout<<"id: "<<id<<endl;
        cout<<"----------------------"<<endl;
        cout<<"before: "<<endl;
        for(unsigned int i = 0; i < collision_intervals[id].size(); i++)
        {
            cout<<collision_intervals[id][i].first<<"~"<<collision_intervals[id][i].second<<endl;
            prt_double(collision_intervals[id][i].second);
        }
    }
    for(unsigned int i = 0; i + 1 < collision_intervals[id].size(); i++)
    {
        if (ge_raw(collision_intervals[id][i].second, collision_intervals[id][i+1].second)) // >=
        {
            /*
               cout<<"comb A"<<endl;
               cout<<collision_intervals[id][i].second<<"  "<<collision_intervals[id][i+1].second<<endl;
               cout<<gt_raw(collision_intervals[id][i].second, collision_intervals[id][i+1].second)<<endl;
               cout<<eq_raw(collision_intervals[id][i].second, collision_intervals[id][i+1].second)<<endl;
               cout<<collision_intervals[id][i].second-collision_intervals[id][i+1].second<<endl;
               cout<<CN_EPSILON<<endl;
               */
            collision_intervals[id].erase(collision_intervals[id].begin() + i + 1);
            i--;
        }
        else if(ge_raw(collision_intervals[id][i].second, collision_intervals[id][i+1].first))  //>=
        {
            //cout<<"comb B"<<endl;
            collision_intervals[id][i].second = collision_intervals[id][i+1].second;
            collision_intervals[id].erase(collision_intervals[id].begin() + i + 1);
            i--;
        }
    }
    if(config.debug>1 || p)
    {
        cout<<"after: "<<endl;
        for(unsigned int i = 0; i < collision_intervals[id].size(); i++)
        {
            cout<<collision_intervals[id][i].first<<"~"<<collision_intervals[id][i].second<<endl;
            prt_double(collision_intervals[id][i].second);
        }
        cout<<"----------------------"<<endl;
    }
}

void SIPP::add_move_constraint(Move move)
{
    //cout<<"inserting move: "<<move.t1<<"~"<<move.t2<<endl;
    std::vector<Interval> m_cons(0);
    if(constraints.count({move.id1, move.id2}) == 0)
        constraints.insert({{move.id1, move.id2}, {make_pair(move.t1,move.t2)}});
    else
    {
        m_cons = constraints.at({move.id1, move.id2});
        bool inserted(false);
        for(unsigned int i = 0; i < m_cons.size(); i++)
        {
            //cout<<m_cons[i].t1<<"~"<<m_cons[i].t2<<endl;
            if(inserted)
                break;
            if(gt_raw(m_cons[i].first, move.t1))  //>
            {
                //cout<<"got here A"<<endl;
                if(le_raw(m_cons[i].second, move.t2))  //<=
                {
                    //cout<<"got here B"<<endl;
                    m_cons[i].first = move.t1;
                    if(ge_raw(move.t2, m_cons[i].second)) // >=
                        m_cons[i].second = move.t2;
                    inserted = true;
                    if(i != 0)
                        if(ge_raw(m_cons[i-1].second, move.t1) && le_raw(m_cons[i-1].second, move.t2)) //>= && <=
                        {
                            m_cons[i-1].second = move.t2;
                            if(ge_raw(m_cons[i-1].second, m_cons[i].first) && le_raw(m_cons[i-1].second, m_cons[i].second)) //>= && <=
                            {
                                m_cons[i-1].second = m_cons[i].second;
                                m_cons.erase(m_cons.begin() + i);
                            }
                            inserted = true;
                        }
                }
                else
                {
                    //cout<<"got here C"<<endl;
                    if(i != 0)
                        if(gt_raw(m_cons[i-1].second, move.t1) && le_raw(m_cons[i-1].second, move.t2)) //>= && <=
                        {
                            m_cons[i-1].second = move.t2;
                            inserted = true;
                            break;
                        }
                    m_cons.insert(m_cons.begin() + i, {move.t1,move.t2});
                    inserted = true;
                }
            }
        }
        //cout<<m_cons.back().t1<<"~"<<m_cons.back().t2<<endl;
        //cout<<ge_raw(m_cons.back().t2 ,move.t1)<<" "<<le_raw(m_cons.back().t2, move.t2)<<endl;
        if(ge_raw(m_cons.back().second ,move.t1) && le_raw(m_cons.back().second, move.t2))  //>= && <=
            m_cons.back().second = move.t2;
        else if(!inserted)
            m_cons.push_back({move.t1,move.t2});
        constraints.at({move.id1, move.id2}) = m_cons;
    }
}

void SIPP::make_constraints(std::list<Constraint> &cons,const Map &map)
{
    for(auto con : cons)
    {
        if(con.positive == false)
        {
            if(con.id2 == -1) // wait consatraint
                add_collision_interval(con.id1, std::make_pair(con.t1, con.t2));
            else
                add_move_constraint(Move(con));
        }
        else
        {
            bool inserted = false;
            int id2;
            if (con.id2==-1)
               id2=con.id1;
            else
                id2=map.ind2id(con.id1,con.id2);
            for(unsigned int i = 0; i < landmarks.size(); i++)
                if(gt_raw(landmarks[i].t1, con.t1))  // >
                {
                    landmarks.insert(landmarks.begin() + i, Move(con.t1, con.t2, con.id1, id2));
                    inserted = true;
                    break;
                }
            if(!inserted)
                landmarks.push_back(Move(con.t1, con.t2, con.id1, id2));
        }
    }
}

Path SIPP::add_part(Path result, Path part)
{
    part.nodes.erase(part.nodes.begin());
    for(auto n: part.nodes)
        result.nodes.push_back(n);
    return result;
}

std::vector<Path> SIPP::find_partial_path(std::vector<Node> starts, std::vector<Node> goals, const Map &map, Heuristic &h_values, double max_f)
{
    open.clear();
    close.clear();
    path.cost = -1;
    visited.clear();
    call_find_pp+=1;

    std::vector<Path> paths(goals.size());
    int pathFound(0);
    for(auto s:starts)
    {
        s.parent = nullptr;
        if(goals[0].id == agent.goal_id) //perfect heuristic is known
            s.f = s.g + h_values.get_value(s.id, agent.id);
        else
        {
            double h = map.get_dist(goals[0].id,s.id);//sqrt(pow(goal.i - newNode.i, 2) + pow(goal.j - newNode.j, 2));
            for(unsigned int i = 0; i < h_values.get_size(); i++) //differential heuristic with pivots placed to agents goals
                h = std::max(h, fabs(h_values.get_value(s.id, i) - h_values.get_value(goals[0].id, i)));
            s.f = s.g + h;
        }
        //prt_node(s);
        add_open(s);
        node_idx temp_idx(make_tuple(s.id,s.interval_id,s.from_landMark));
        visited.insert({temp_idx, {s.g, false}});
    }
    Node curNode;
    while(!open.empty())
    {
        curNode = find_min();
        node_exp+=1;
        if(config.debug>1 || p)
        {
            cout<<"###########"<<endl;
            cout<<"poped: llid: "<<node_exp<<endl;
            prt_node(curNode);
        }
        auto v = visited.find(make_tuple(curNode.id,curNode.interval_id,curNode.from_landMark));
        if(v->second.second){
            continue;
        }
        if (!curNode.from_landMark)
        {
            v->second.second = true;
        }

        if(p)
        {
            cout<<"curNode.g"<<curNode.g<<endl;
            prt_double(curNode.g);
            cout<<endl;
        }
        auto temp=make_tuple(curNode.id,curNode.interval_id,curNode.from_landMark);
        boost::unordered_map<node_idx, Node>::iterator close_it;
        Node* parent;
        close_it=close.find(temp);
        if (close_it==close.end()) {
            parent = &close.insert({temp, curNode}).first->second;
        }
        else {
            close_it->second=curNode;
            parent=&(close_it->second);
        }

        if(p)
        {
            cout<<"parent g"<<parent->g<<endl;
            prt_double(parent->g);
            cout<<endl;
            if(curNode.id==88)
            {
                curNode.id+=1;
                curNode.id-=1;
            }
        }
        if(curNode.id == goals[0].id)
        {
            for(unsigned int i = 0; i < goals.size(); i++)
            {

                bool isGoal;
                if (goals[i].interval.first==goals[i].interval.second)
                    isGoal=(curNode.g==goals[i].interval.first);
                else
                    isGoal=(lt_raw(curNode.g, goals[i].interval.second) && le_raw(goals[i].interval.first,curNode.interval.second));
                if(isGoal && (paths[i].cost==-1 || paths[i].cost>curNode.g))//< && <=
                {
                    paths[i].nodes = reconstruct_path(curNode);
                    if(p)
                    {
                        cout<<"solution: "<<endl;
                        prt_nodes(paths[i].nodes);
                    }
                    if(paths[i].nodes.back().g < goals[i].interval.first)
                    {
                        curNode.g = goals[i].interval.first;
                        paths[i].nodes.push_back(curNode);
                    }
                    paths[i].cost = curNode.g;
                    paths[i].expanded = int(close.size());
                    pathFound++;
                }
            }
            if(pathFound == int(goals.size()))
            {
                return paths;
            }
        }
        std::list<Node> succs;
        succs.clear();
        find_successors(curNode, map, succs, h_values, Node(goals[0].id, 0, 0, goals[0].i, goals[0].j));
        std::list<Node>::iterator it = succs.begin();
        while(it != succs.end())
        {
            if(gt(it->f, max_f)) //>
            {
                it++;
                continue;
            }
            it->parent = parent;
            add_open(*it);
            node_open+=1;
            it++;
        }

        succs.clear();
    }
    return paths;
}

std::vector<Node> SIPP::get_endpoints(int node_id, double node_i, double node_j, double t1, double t2)
{
    std::vector<Node> nodes;
    int id(0);
    nodes = {Node(node_id, 0, 0, node_i, node_j, nullptr, t1, t2,id++)};
    if(collision_intervals[node_id].empty())
        return nodes;
    else
        for(unsigned int k = 0; k < collision_intervals[node_id].size(); k++)
        {
            unsigned int i(0);
            while(i < nodes.size())
            {
                Node n = nodes[i];
                auto c = collision_intervals[node_id][k];
                //cout<<"n: "<<n.interval.first<<"~"<<n.interval.second<<endl;
                //cout<<"c: "<<c.first<<"~"<<c.second<<endl;
                bool changed = false;
                if(le_raw(c.first, n.interval.first) && ge_raw(c.second, n.interval.second)) //<= && >=
                {
                    //cout<<"A"<<endl;
                    nodes.erase(nodes.begin() + i);
                    changed = true;
                }
                else if(le_raw(c.first, n.interval.first) && (gt_raw(c.second, n.interval.first))) //<= && >
                {
                    //cout<<"B"<<endl;
                    nodes[i].interval.first = c.second;
                    changed = true;
                }
                else if(gt_raw(c.first, n.interval.first) && lt_raw(c.second, n.interval.second)) //> && <
                {
                    //cout<<"C"<<endl;
                    nodes[i].interval.second = c.first;
                    nodes.insert(nodes.begin() + i + 1, Node(node_id, 0, 0, node_i, node_j, nullptr, c.second, n.interval.second,id++,true));
                    changed = true;
                }
                else if((lt_raw(c.first, n.interval.second)) && ge_raw(c.second, n.interval.second)) //< && >=
                {
                    //cout<<"D"<<endl;
                    nodes[i].interval.second = c.first;
                    changed = true;
                }
                if(changed)
                {
                    i = -1;
                    k = 0;
                }
                i++;
            }
        }
    return nodes;
}

pair<double,double> SIPP::check_endpoint(Node start, Node goal,const Map &map)
{
    //double cost = sqrt(pow(start.i - goal.i, 2) + pow(start.j - goal.j, 2));
    double cost=map.get_dist(start.id,goal.id);
    double ind=map.id2ind(start.id,goal.id);
    if(lt_raw(start.g+cost, goal.interval.first)) //<=
        start.g = goal.interval.first - cost;

    if(constraints.count({start.id, ind}) != 0)
    {
        auto it = constraints.find({start.id, ind});
        for(unsigned int i = 0; i < it->second.size(); i++)
        {
            if (config.debug>1 || p)
            {
                cout<<start.interval.first<<"~"<<start.interval.second<<endl;
               cout<<"start g: "<<start.g<<" it t1: "<<it->second[i].first<<endl;
               cout<<"g:     0b";
               prt_double(start.g);

               std::cout << std::endl;
               cout<<"it t1: 0b";
               prt_double(it->second[i].first);
               std::cout << std::endl;
            }
            if(ge_raw(start.g, it->second[i].first)  && lt_raw(start.g, it->second[i].second)) //>= && <
            {
                //cout<<"got here"<<endl;
                start.g = it->second[i].second;
                //cout<<"new start g"<<start.g<<endl;
                //prt_double(start.g);
                //cout<<endl;
            }
        }
    }
    if(ge_raw(start.g, start.interval.second)  || ge_raw(start.g+cost, goal.interval.second))
        return {CN_INFINITY,CN_INFINITY};
    else
        return {start.g, start.g + cost};
}

Path SIPP::find_path(Agent agent, const Map &map, std::list<Constraint> cons, Heuristic &h_values, bool p)
{
    auto temp=find_path_aux(agent, map, cons, h_values, p);
    if(config.debug>1 || p)
        prt_info();
    return temp;
}

void SIPP::prt_intervals()
{
    //std::unordered_map<int, std::vector<std::pair<double, double>>>::iterator it;//stores sets of collision nextIntervals associated with cells
    for (auto it = collision_intervals.begin(); it != collision_intervals.end(); it++)
    {
        std::cout << it->first<< ":[ ";
        for (auto m : it->second){
            std::cout<<m.first<<"~"<<m.second<<" + ";
        }
        std::cout<<"]"<<std::endl;
    }
    cout<<"---------------"<<endl;
}

void SIPP::prt_landmarks()
{
    for (Move m: landmarks){
        cout<<"from "<<m.id1<<" -> "<<m.id2<<"[t:"<<m.t1<<"~"<<m.t2<<"]"<<endl;
    }
}

void SIPP::prt_cons()
{
    std::map<std::pair<int, int>, std::vector<Interval>>::iterator it;//stores sets of constraints associated with moves

    for (it = constraints.begin(); it != constraints.end(); it++)
    {
        std::cout << it->first.first<<"->"<<it->first.second
            << ":[ ";
        for (Interval m : it->second){
            std::cout<<m.first<<"~"<<m.second<<" + ";
            cout<<endl;
            prt_double(m.first);
            cout<<" ~"<<endl;
            prt_double(m.second);
            cout<<" +"<<endl;
        }

        std::cout<<"]"<<std::endl;
    }
}

void SIPP::prt_node(Node n)
{
    //cout<<"id:"<<n.id<<", f:"<<n.f<<" g:"<<n.g<<", i:"<<n.i<<", j:"<<n.j<<" interval:"<<n.interval_id<<"("<<n.interval.first<<" , "<<n.interval.second<<") "<<(n.from_landMark ? "++" : "--")<<endl;
    cout<<"id:"<<n.id<<", f:"<<n.f<<" g:"<<n.g<<", prev_g:"<<n.prev_g<<" interval:"<<n.interval_id<<"("<<n.interval.first<<" , "<<n.interval.second<<") "<<(n.from_landMark ? "++" : "--")<<endl;
}

void SIPP::prt_nodes(std::vector<Node> nodes)
{
    for (Node n:nodes){
        prt_node(n);
        prt_double(n.g);
        cout<<endl;
    }
}

void SIPP::prt_constraints(std::list<Constraint> constraints){
    for(Constraint c:constraints)
        std::cout<<"Constraint "<<c.positive<<" a:"<<c.agent<<" from:"<<c.id1<<"to:"<<c.id2<<"[t:"<<c.t1<<"~"<<c.t2<<"]"<<std::endl;
}

Path SIPP::find_path_aux(Agent agent, const Map &map, std::list<Constraint> cons, Heuristic &h_values, bool p)
{
    //cout<<"start planning"<<endl;
    this->clear();
    this->agent = agent;
    this->p=p;
    this->call_find_pp=0;
    this->node_open=0;
    this->node_exp=0;
    if(config.debug>1 || p)
    {
        cout<<"inital constraints:"<<endl;
        prt_constraints(cons);
    }
    make_constraints(cons,map);
    if(config.debug>1 || p)
    {
        cout<<"converted: "<<endl;
        prt_cons();
        cout<<"unsafe Intervals: "<<endl;
        prt_intervals();
        cout<<"landmarks: "<<endl;
        prt_landmarks();
    }
    //cout<<"nextIntervals: "<<endl;
    //prt_intervals();
    //cout<<"finish make constraint"<<endl;

    std::vector<Node> starts, goals;
    std::vector<Path> parts, results, new_results;
    Path part, result;
    int expanded(0);
    if(!landmarks.empty())
    {
        for(unsigned int i = 0; i <= landmarks.size(); i++)
        {
            if(i == 0)
            {
                starts = {get_endpoints(agent.start_id, agent.start_i, agent.start_j, 0, CN_INFINITY).at(0)};
                if (starts.at(0).interval.first!=0)
                    return Path();
                goals = get_endpoints(landmarks[i].id1, map.get_i(landmarks[i].id1), map.get_j(landmarks[i].id1), landmarks[i].t1, landmarks[i].t2);
            }
            else
            {
                starts.clear();
                for(auto p:results)
                    starts.push_back(p.nodes.back());
                if(i == landmarks.size())
                    goals = {get_endpoints(agent.goal_id, agent.goal_i, agent.goal_j, 0, CN_INFINITY).back()};
                else
                {
                    goals = get_endpoints(landmarks[i].id1, map.get_i(landmarks[i].id1), map.get_j(landmarks[i].id1), landmarks[i].t1, landmarks[i].t2);
                }
            }
            if(goals.empty())
            {
                return Path();
            }
            if(config.debug>1 || p)
            {
                cout<<"start planning"<<endl;
                cout<<"start: "<<endl;
                prt_nodes(starts);
                cout<<"goal: "<<endl;
                prt_nodes(goals);
            }
            parts = find_partial_path(starts, goals, map, h_values, goals.back().interval.second);
            if(config.debug>1 || p)
            {
                cout<<"finish planning"<<endl;
                for (Path p:parts)
                {
                    sPath tempPath;
                    tempPath=p;
                    cout<<tempPath<<endl;
                }
            }
            expanded += int(close.size());
            new_results.clear();
            if(i == 0)
                for(unsigned int k = 0; k < parts.size(); k++)
                {
                    if(parts[k].nodes.empty())
                        continue;
                    new_results.push_back(parts[k]);
                }
            for(unsigned int k = 0; k < parts.size(); k++)
                for(unsigned int j = 0; j < results.size(); j++)
                {
                    if(parts[k].nodes.empty())
                        continue;
                    if(eq_raw(parts[k].nodes[0].interval.first, results[j].nodes.back().interval.first) && eq_raw(parts[k].nodes[0].interval.second, results[j].nodes.back().interval.second))
                    {
                        new_results.push_back(results[j]);
                        new_results.back() = add_part(new_results.back(), parts[k]);
                        /*
                           cout<<"old part"<<endl;
                           prt_nodes(results[j].nodes);
                           cout<<"add part"<<endl;
                           prt_nodes(new_results.back().nodes);
                           */
                    }
                }
            results = new_results;
            if(results.empty())
                return Path();
            if(i < landmarks.size())
            {
                starts.clear();
                for(auto p:results)
                    starts.push_back(p.nodes.back());
                /*
                   cout<<"starts"<<endl;
                   prt_nodes(starts);
                   */
                //double offset = sqrt(pow(map.get_i(landmarks[i].id1) - map.get_i(landmarks[i].id2), 2) + pow(map.get_j(landmarks[i].id1) - map.get_j(landmarks[i].id2), 2));
                double offset=map.get_dist(landmarks[i].id1,landmarks[i].id2);
                goals = get_endpoints(landmarks[i].id2, map.get_i(landmarks[i].id2), map.get_j(landmarks[i].id2), landmarks[i].t1 + offset, landmarks[i].t2 + offset);
                /*
                   cout<<"goals"<<endl;
                   prt_nodes(goals);
                   */
                if(goals.empty())
                    return Path();
                new_results.clear();
                for(unsigned int k = 0; k < goals.size(); k++)
                {
                    double best_g(CN_INFINITY),best_prev_g(CN_INFINITY);
                    int best_start_id = -1;
                    for(unsigned int j = 0; j < starts.size(); j++)
                    {
                        pair<double,double> g_val = check_endpoint(starts[j], goals[k],map);
                        /*
                           cout<<"picking g"<<endl;
                           cout<<g<<endl;
                           prt_double(g);
                           cout<<endl;
                           */
                        if(g_val.second < best_g)
                        {
                            best_start_id = j;
                            best_g = g_val.second;
                            //cout<<"best_g: "<<g_val.first<<endl;
                            //prt_double(g_val.first);
                            //cout<<endl;
                            //cout<<"best_preb_g: "<<g_val.second<<endl;
                            //prt_double(g_val.second);
                            //cout<<endl;
                            best_prev_g=g_val.first;
                        }
                    }
                    if(best_start_id >= 0)
                    {
                        goals[k].g = best_g;
                        goals[k].interval.first=best_g;
                        //cout<<"best_g: " <<best_g<<endl;
                        goals[k].prev_g=best_prev_g;
                        /*
                           cout<<"prev_g: " <<goals[k].prev_g<<endl;
                           prt_double(goals[k].prev_g);
                           */
                        if(collision_intervals[goals[k].id].empty())
                            goals[k].interval.second = CN_INFINITY;
                        else
                        {
                            goals[k].interval.second = CN_INFINITY;
                            for(auto c:collision_intervals[goals[k].id])
                                if(goals[k].g < c.first)
                                {
                                    goals[k].interval.second = c.first;
                                    break;
                                }
                        }
                        new_results.push_back(results[best_start_id]);
                        /*
                           cout<<"old back"<<endl;
                           prt_node(new_results.back().nodes.back());
                           prt_double(new_results.back().nodes.back().g);
                           */
                        if(gt_raw(goals[k].g,starts[best_start_id].g+offset) && !eq_raw(goals[k].g,new_results.back().nodes.back().g+offset))
                        {

                            //cout<<"add extra wait: "<<endl;
                            new_results.back().nodes.push_back(new_results.back().nodes.back());
                            new_results.back().nodes.back().g = goals[k].prev_g;
                            //new_results.back().nodes.back().prev_g = goals[k].prev_g;
                            new_results.back().nodes.back().interval.first = goals[k].prev_g;
                            //cout<<"new back1"<<endl;
                            //prt_node(new_results.back().nodes.back());
                        }
                        new_results.back().nodes.push_back(goals[k]);
                        if(config.debug>1 || p)
                        {
                            cout<<"new ressult "<<endl;
                            for (Path p:parts)
                            {
                                sPath tempPath;
                                tempPath=p;
                                cout<<tempPath<<endl;
                            }
                        }
                        //cout<<"new back2"<<endl;
                        //prt_node(new_results.back().nodes.back());
                    }
                }

                results = new_results;
                if(results.empty())
                    return Path();
            }
        }
        result = results[0];
    }
    else
    {

        starts = {get_endpoints(agent.start_id, agent.start_i, agent.start_j, 0, CN_INFINITY).at(0)};
        goals = {get_endpoints(agent.goal_id, agent.goal_i, agent.goal_j, 0, CN_INFINITY).back()};
        //cout<<"goals: "<<endl;
        //prt_nodes(goals);
        parts = find_partial_path(starts, goals, map, h_values);
        expanded = int(close.size());
        if(parts[0].cost < 0)
            return Path();
        result = parts[0];
    }
    result.nodes.shrink_to_fit();
    result.cost = result.nodes.back().g;
    result.agentID = agent.id;
    result.expanded = expanded;
    if(config.debug>1 || p)
    {
        cout<<"returned path"<<endl;
        prt_nodes(result.nodes);

    }
    return result;
}
