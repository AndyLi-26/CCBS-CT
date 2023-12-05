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
    if (config.debug>1)
    {
        std::cout<<"expanding: "<<curNode.id<<" g:"<<curNode.g<<"interval:("<<curNode.interval.first<<", "<<curNode.interval.second<<") "<<std::endl;
        for (Node n:valid_moves){
            std::cout<<"->("<<n.id<<","<<dist(curNode,n)<<")";
        }
        std::cout<<std::endl;
    }
    for(int i=0;i<valid_moves.size();++i)
    {
        Node move = valid_moves[i];
        if (config.debug>1)
            std::cout<<"considering: "<<move.id<<"  ->";
        //if (move.agent.find(-1)==move.agent.end() && move.agent.find(agent.id)==move.agent.end())
        //  continue;
        newNode.i = move.i;
        newNode.j = move.j;
        newNode.id = move.id;
        double cost = dist(curNode, newNode);
        newNode.g = curNode.g + cost;
        std::vector<std::pair<double, double>> intervals(0);
        auto colls_it = collision_intervals.find(newNode.id);
        if(colls_it != collision_intervals.end())
        {
            std::pair<double, double> interval = {0, CN_INFINITY};
            for(unsigned int i = 0; i < colls_it->second.size(); i++)
            {
                interval.second = colls_it->second[i].first;
                intervals.push_back(interval);
                interval.first = colls_it->second[i].second;
            }
            interval.second = CN_INFINITY;
            intervals.push_back(interval);
        }
        else
            intervals.push_back({0, CN_INFINITY});
        if (config.debug>1)
        {
            for (auto intv:intervals){
                std::cout<<"("<<intv.first<<"~"<<intv.second<<")";
            }
            std::cout<<"->";
        }
        //int to_index=map.valid_moves
        auto cons_it = constraints.find({curNode.id,i});
        int id(0);
        for(auto interval: intervals)
        {
            newNode.interval_id = id;
            id++;
            auto it = visited.find(make_tuple(newNode.id,newNode.interval_id,false));
            if(it != visited.end())
                if(it->second.second)
                    continue;
            if(lt(interval.second, newNode.g)) // <
                continue;
            if(gt(interval.first, newNode.g)) // >
                newNode.g = interval.first;
            if(cons_it != constraints.end())
                for(unsigned int i = 0; i < cons_it->second.size(); i++)
                {
                    //cout<<"g: "<<newNode.g<<" cost: "<<cost<<endl;
                    //cout<<"g-cost: "<<newNode.g-cost<<" cons_it->second[i].t1: "<<cons_it->second[i].t1<<" cons_it->second[i].t2: "<<cons_it->second[i].t2<<endl;
                    //cout<<": ="<<(newNode.g - cost + CN_EPSILON > cons_it->second[i].t1)<<" & "<<(newNode.g - cost < cons_it->second[i].t2 + CN_EPSILON)<<endl;
                    if(ge(newNode.g - cost, cons_it->second[i].t1) && lt(newNode.g - cost, cons_it->second[i].t2)) // >= && <
                        newNode.g = cons_it->second[i].t2 + cost;
                    //cout<<"new g: "<<newNode.g<<endl;
                }
            newNode.interval = interval;
            if(gt(newNode.g - cost ,curNode.interval.second) || gt(newNode.g, newNode.interval.second)) // > || >
                continue;
            if(it != visited.end())
            {
                if(le(it->second.first, newNode.g)) // <=
                    continue;
                else
                    it->second.first = newNode.g;
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
        }
        if (config.debug>1)
        {
            for (Node n:succs){
                std::cout<<"("<<n.id<<"@"<<n.g<<")";
            }
            std::cout<<"]"<<std::endl;
        }
    }
}

Node SIPP::find_min()
{
    Node min = *open.begin();
    open.pop_front();
    return min;
}

void SIPP::add_open(Node newNode)
{
    if (open.empty() || le(open.back().f, newNode.f))
    {
        open.push_back(newNode);
        return;
    }
    for(auto iter = open.begin(); iter != open.end(); ++iter)
    {
        if(gt(iter->f, newNode.f)) // if newNode.f has lower f-value  >
        {
            open.insert(iter, newNode);
            return;
        }
        else if(eq(iter->f,newNode.f) && ge(newNode.g, iter->g)) // if f-values are equal, compare g-values  == && >=
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
            path.nodes.insert(path.nodes.begin(), curNode);
            curNode = *curNode.parent;
        }
        while(curNode.parent != nullptr);
    path.nodes.insert(path.nodes.begin(), curNode);
    //store sipp path
    path.SIPP_nodes.clear();
    for (Node n:path.nodes)
        path.SIPP_nodes.push_back(n);
    for(unsigned int i = 0; i < path.nodes.size(); i++)
    {
        //cout<<"new i: "<<i<<endl<<flush;
        //assert(i<=10);
        unsigned int j = i + 1;
        if(j == path.nodes.size())
            break;
        if(!eq(path.nodes[j].g - path.nodes[i].g, dist(path.nodes[j], path.nodes[i]))) //==
        {
            //cout<<path.nodes[i].id<<"->"<<path.nodes[j].id<<" : "<<path.nodes[i].g<<"~"<<path.nodes[j].g<<" dis:"<<dist(path.nodes[j], path.nodes[i])<<endl<<flush;
            Node add = path.nodes[i];
            add.g = path.nodes[j].g - dist(path.nodes[j], path.nodes[i]);
            path.nodes.emplace(path.nodes.begin() + j, add);
        }
    }
    return path.nodes;
}

void SIPP::add_collision_interval(int id, std::pair<double, double> interval)
{
    std::vector<std::pair<double, double>> intervals(0);
    if(collision_intervals.count(id) == 0){
        collision_intervals.insert({id, {interval}});
    }
    else
    {
        collision_intervals[id].push_back(interval);
    }
    std::sort(collision_intervals[id].begin(), collision_intervals[id].end());
    for(unsigned int i = 0; i + 1 < collision_intervals[id].size(); i++)
    {
        if (ge(collision_intervals[id][i].second, collision_intervals[id][i+1].second)) // >=
        {
            collision_intervals[id].erase(collision_intervals[id].begin() + i + 1);
            i--;
        }
        else if(ge(collision_intervals[id][i].second, collision_intervals[id][i+1].first))  //>=
        {
            collision_intervals[id][i].second = collision_intervals[id][i+1].second;
            collision_intervals[id].erase(collision_intervals[id].begin() + i + 1);
            i--;
        }
    }
}

void SIPP::add_move_constraint(Move move)
{
    std::vector<Move> m_cons(0);
    if(constraints.count({move.id1, move.id2}) == 0)
        constraints.insert({{move.id1, move.id2}, {move}});
    else
    {
        m_cons = constraints.at({move.id1, move.id2});
        bool inserted(false);
        for(unsigned int i = 0; i < m_cons.size(); i++)
        {
            if(inserted)
                break;
            if(gt(m_cons[i].t1, move.t1))  //>
            {
                if(le(m_cons[i].t1, move.t2))  //<=
                {
                    m_cons[i].t1 = move.t1;
                    if(ge(move.t2, m_cons[i].t2)) // >=
                        m_cons[i].t2 = move.t2;
                    inserted = true;
                    if(i != 0)
                        if(ge(m_cons[i-1].t2, move.t1) && le(m_cons[i-1].t2, move.t2)) //>= && <=
                        {
                            m_cons[i-1].t2 = move.t2;
                            if(ge(m_cons[i-1].t2, m_cons[i].t1) && le(m_cons[i-1].t2, m_cons[i].t2)) //>= && <=
                            {
                                m_cons[i-1].t2 = m_cons[i].t2;
                                m_cons.erase(m_cons.begin() + i);
                            }
                            inserted = true;
                        }
                }
                else
                {
                    if(i != 0)
                        if(gt(m_cons[i-1].t2, move.t1) && le(m_cons[i-1].t2, move.t2)) //>= && <=
                        {
                            m_cons[i-1].t2 = move.t2;
                            inserted = true;
                            break;
                        }
                    m_cons.insert(m_cons.begin() + i, move);
                    inserted = true;
                }
            }
        }
        if(ge(m_cons.back().t2 ,move.t1) && le(m_cons.back().t2, move.t2))  //>= && <=
            m_cons.back().t2 = move.t2;
        else if(!inserted)
            m_cons.push_back(move);
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
            //auto temp=map.get_valid_moves(con.id1,agent.id);
            //int id2=temp[con.id2].id;
            for(unsigned int i = 0; i < landmarks.size(); i++)
                if(gt(landmarks[i].t1, con.t1))  // >
                {
                    landmarks.insert(landmarks.begin() + i, Move(con.t1, con.t2, con.id1, con.id2));
                    inserted = true;
                    break;
                }
            if(!inserted)
                landmarks.push_back(Move(con.t1, con.t2, con.id1, con.id2));
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
    //cout<<"start planning partial: "<<endl;
    open.clear();
    close.clear();
    path.cost = -1;
    visited.clear();
    std::vector<Path> paths(goals.size());
    int pathFound(0);
    //cout<<"start node in partial planner:"<<endl;
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
        //prt_open();
        curNode = find_min();
        auto v = visited.find(make_tuple(curNode.id,curNode.interval_id,curNode.from_landMark));
        if(v->second.second){
            continue;
        }
        if (config.debug>1)
        {
            cout<<endl<<endl<<"***"<<endl<<"poped: ";
            prt_node(curNode);
        }
        if (!curNode.from_landMark)
        {
            v->second.second = true;
        }
        auto parent = &close.insert({make_tuple(curNode.id,curNode.interval_id,curNode.from_landMark), curNode}).first->second;
        if(curNode.id == goals[0].id)
        {
            for(unsigned int i = 0; i < goals.size(); i++)
            {
                if(lt(curNode.g, goals[i].interval.second) && le(goals[i].interval.first,curNode.interval.second))//< && <=
                {
                    paths[i].nodes = reconstruct_path(curNode);
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
                return paths;
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
            it++;
        }

        succs.clear();
    }
    return paths;
}

std::vector<Node> SIPP::get_endpoints(int node_id, double node_i, double node_j, double t1, double t2)
{
    //cout<<"received t1: "<<t1<<" t2:"<<t2<<endl;
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
                bool changed = false;
                if(le(c.first, n.interval.first) && ge(c.second, n.interval.second)) //<= && >=
                {
                    //cout<<"A"<<endl;
                    nodes.erase(nodes.begin() + i);
                    changed = true;
                }
                else if(le(c.first, n.interval.first) && gt(c.second, n.interval.first)) //<= && >
                {
                    //cout<<"B"<<endl;
                    nodes[i].interval.first = c.second;
                    changed = true;
                }
                else if(gt(c.first, n.interval.first) && lt(c.second, n.interval.second)) //> && <
                {
                    //cout<<"C"<<endl;
                    //cout<<"c.f: "<<c.first<<" n.int.f:"<< n.interval.first<<"   c.s:"<<c.second<<" n.int.s: "<<n.interval.second<<endl;
                    nodes[i].interval.second = c.first;
                    nodes.insert(nodes.begin() + i + 1, Node(node_id, 0, 0, node_i, node_j, nullptr, c.second, n.interval.second,id++,true));
                    changed = true;
                }
                else if(lt(c.first, n.interval.second) && ge(c.second, n.interval.second)) //< && >=
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

double SIPP::check_endpoint(Node start, Node goal, int exit_id)
{
    if(config.debug>1)
        cout<<start.id<<"@"<<start.g<<"->"<<exit_id<<"("<<goal.id<<")@"<<goal.g<<"( "<<goal.interval.first<<","<<goal.interval.second<<") "<<endl;
    double cost = sqrt(pow(start.i - goal.i, 2) + pow(start.j - goal.j, 2));
    if(le(start.g + cost, goal.interval.first)) //<=
        start.g = goal.interval.first - cost;
    if(config.debug>1)
        cout<<start.id<<"@"<<start.g<<"->"<<exit_id<<"("<<goal.id<<")@"<<goal.g<<endl;

    if(constraints.count({start.id, exit_id}) != 0)
    {
        auto it = constraints.find({start.id,exit_id});
        for(unsigned int i = 0; i < it->second.size(); i++)
        {
            if(config.debug>1)
            {
                cout<<"start g:"<<start.g<<" it->second[i].t1:"<<it->second[i].t1<<" it->second[i].t2:"<<it->second[i].t2<<endl;
                //cout<<": ="<<(start.g +CN_EPSILON > it->second[i].t1) <<" and "<< (start.g < it->second[i].t2) <<endl;
            }
            if(ge(start.g, it->second[i].t1)  && le(start.g, it->second[i].t2)) //>= && <
            {
                start.g = it->second[i].t2;
                if(config.debug>1)
                {
                    cout<<"changed g:"<<endl;
                    cout<<start.id<<"@"<<start.g<<"->"<<exit_id<<"("<<goal.id<<")@"<<goal.g<<endl;
                }
            }
        }
    }
    //cout<<"start.g: "<<start.g<<"start.interval.second: "<<start.interval.second<<" goal.interval.second"<<goal.interval.second<<endl;
    if(ge(start.g, start.interval.second)  || ge(start.g + cost, goal.interval.second)) // >= && >=
        return CN_INFINITY;
    else
        return start.g + cost;
}


Path SIPP::find_path(Agent agent, const Map &map, std::list<Constraint> cons, Heuristic &h_values, bool p){
    this->clear();
    this->agent = agent;
    if (config.debug>1){
        std::cout<<"planning:"<<std::endl;
        std::cout<<agent.id<<std::endl;
        cout<<"map:"<<endl;
        map.prt_validmoves();
        cout<<"constraints:"<<endl;
        prt_constraints(cons);
        cout<<"h:"<<endl;
        h_values.prt();
    }
    this->p =p;
    make_constraints(cons,map);
    if (config.debug>1){
        cout<<"constraints: "<<endl;
        prt_cons();
        cout<<"intervals: "<<endl;
        prt_intervals();
        cout<<"landmarks: "<<endl;
        prt_landmarks();
    }
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
                return Path();
            if (config.debug>1){
                cout<<"starts before plan"<<endl;
                prt_nodes(starts);

                cout<<"goals"<<endl;
                for (Node n : goals)
                {
                    cout<<n.id<<": ("<<n.interval.first<<","<<n.interval.second<<")"<<endl;
                }
            }
            parts = find_partial_path(starts, goals, map, h_values, goals.back().interval.second);
            if (config.debug>1){
                cout<<"parts"<<endl;
                prt_paths(parts);
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
                    if(eq(parts[k].nodes[0].interval.first, results[j].nodes.back().interval.first) && eq(parts[k].nodes[0].interval.second, results[j].nodes.back().interval.second))
                    {
                        new_results.push_back(results[j]);
                        new_results.back() = add_part(new_results.back(), parts[k]);
                    }
                }
            results = new_results;

            if (config.debug>1){
                cout<<"results:"<<endl;
                prt_paths(results);
            }
            if(results.empty())
                return Path();
            if(i < landmarks.size())
            {
                starts.clear();
                for(auto p:results)
                    starts.push_back(p.nodes.back());

                auto temp=map.get_valid_moves(landmarks[i].id1);
                int id2=temp[landmarks[i].id2].id;
                double offset = sqrt(pow(map.get_i(landmarks[i].id1) - map.get_i(id2), 2) + pow(map.get_j(landmarks[i].id1) - map.get_j(id2), 2));
                //cout<<"id2: "<<id2<<"offset: "<<offset<<endl;
                goals = get_endpoints(id2, map.get_i(id2), map.get_j(id2), landmarks[i].t1 + offset, landmarks[i].t2 + offset);

                if (config.debug>1){
                    cout<<"starts:"<<endl;
                    prt_nodes(starts);
                    cout<<"goals:"<<endl;
                    prt_nodes(goals);
                }
                if(goals.empty())
                    return Path();
                new_results.clear();
                for(unsigned int k = 0; k < goals.size(); k++)
                {
                    double best_g(CN_INFINITY);
                    int best_start_id = -1;
                    for(unsigned int j = 0; j < starts.size(); j++)
                    {
                        double g = check_endpoint(starts[j], goals[k],landmarks[i].id2);
                        if(config.debug>1)
                            cout<<"check_endpoint g:"<<g<<endl;
                        if(g < best_g)
                        {
                            best_start_id = j;
                            best_g = g;
                        }
                    }
                    if (config.debug>1){
                        cout<<"best_g: "<<best_g<<endl;
                    }
                    if(best_start_id >= 0)
                    {
                        goals[k].g = best_g;
                        if(collision_intervals[goals[k].id].empty())
                            goals[k].interval.second = CN_INFINITY;
                        else
                        {
                            for(auto c:collision_intervals[goals[k].id])
                            {
                                //cout<<"g:"<<goals[k].g<<" c:"<<c.first<<endl;
                                if(goals[k].g < c.first)
                                {
                                    goals[k].interval.second = c.first;
                                    break;
                                }
                            }
                        }
                        new_results.push_back(results[best_start_id]);
                        if(config.debug>1)
                        {
                            cout<<"offset: "<<offset<<endl;
                            cout<<"goals[k].g: "<<goals[k].g <<" offset:"<<offset<<" new_results.back().nodes.back().g"<<new_results.back().nodes.back().g<<endl;
                            cout<<": ="<<(goals[k].g-offset)-new_results.back().nodes.back().g<<endl;
                        }
                        if(gt(goals[k].g - starts[best_start_id].g, offset) && !eq((goals[k].g-offset),new_results.back().nodes.back().g))
                        {
                            if(config.debug>1)
                            {
                                cout<<"old path:"<<endl;
                                prt_path(new_results.back());
                            }
                            new_results.back().nodes.push_back(new_results.back().nodes.back());
                            new_results.back().nodes.back().g = goals[k].g - offset;
                            if(config.debug>1)
                            {
                                cout<<"old path:"<<endl;
                                prt_path(new_results.back());
                            }
                        }
                        new_results.back().nodes.push_back(goals[k]);
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
    //Path p=result;
    return result;
}

void SIPP::prt_constraint(Constraint c){
    std::cout<<"Constraint "<<c.positive<<" a:"<<c.agent<<" from:"<<c.id1<<"to:"<<c.id2<<"[t:"<<c.t1<<"~"<<c.t2<<"]"<<std::endl;
}

void SIPP::prt_constraints(std::list<Constraint> constraints){
    std::cout<<std::endl<<"constraints:"<<std::endl;
    for(Constraint c:constraints)
        prt_constraint(c);
}

void SIPP::prt_cons()
{
    std::map<std::pair<int, int>, std::vector<Move>>::iterator it;//stores sets of constraints associated with moves

    for (it = constraints.begin(); it != constraints.end(); it++)
    {
        std::cout << it->first.first<<"->"<<it->first.second
            << ":[ ";
        for (Move m : it->second){
            std::cout<<m.t1<<"~"<<m.t2<<" + ";
        }

        std::cout<<"]"<<std::endl;
    }
}

void SIPP::prt_intervals()
{

    std::unordered_map<int, std::vector<std::pair<double, double>>>::iterator it;//stores sets of collision intervals associated with cells

    for (it = collision_intervals.begin(); it != collision_intervals.end(); it++)
    {
        std::cout << it->first<< ":[ ";
        for (auto m : it->second){
            std::cout<<m.first<<"~"<<m.second<<" + ";
        }


        std::cout<<"]"<<std::endl;
    }
    cout<<"---------------"<<endl;

}

void SIPP::prt_path(Path p)
{
    for (sNode n:p.nodes){
        cout<<"("<<n.id<<","<<n.g<<")->";
    }
    cout<<endl;
}

void SIPP::prt_paths(std::vector<Path> paths)
{
    for(Path p:paths){
        prt_path(p);
        cout<<endl;
    }
}
void SIPP::prt_nodes(std::vector<Node> nodes)
{
    for (Node n:nodes){
        prt_node(n);
        cout<<endl;
    }
}
void SIPP::prt_node(Node n)
{
    cout<<"id:"<<n.id<<", f:"<<n.f<<" g:"<<n.g<<", i:"<<n.i<<", j:"<<n.j<<" interval:"<<n.interval_id<<"("<<n.interval.first<<" , "<<n.interval.second<<") "<<(n.from_landMark ? "++" : "--")<<endl;
}

void SIPP::prt_open()
{
    cout<<"open list:"<<endl;
    cout<<"vvvvvvvvvvvvvvvvvvvvvvv"<<endl;
    for (Node n:open){
        prt_node(n);
    }
    cout<<"^^^^^^^^^^^^^^^^^^^^^^"<<endl;
}
void SIPP::prt_landmarks()
{
    for (Move m: landmarks){
        cout<<"from "<<m.id1<<" -> "<<m.id2<<"[t:"<<m.t1<<"~"<<m.t2<<"]"<<endl;
    }
}
