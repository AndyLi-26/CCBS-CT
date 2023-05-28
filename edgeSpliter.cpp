#include "edgeSpliter.h"
void edgeSpliter::find_deltas(const Conflict &conf, Map_deltas &deltasR, Map_deltas &deltasL, Map &map, Heuristic &h_values)
{
  int node11=conf.move1.id1;
  int node12=conf.move1.id2;
  int node21=conf.move2.id1;
  int node22=conf.move2.id2;
  int a1=conf.agent1,a2=conf.agent2;

  if ((node11==node22) && (node12==node21)) //cross each other
    return;

  Move m1(conf.move1),m2(conf.move2);
  if (node11==node12 && node21==node22) return;
  if (writeNode)
    output.open(config.F_debug_info, std::ios::app);
  if (node11==node12)
  {
    waiting(m1,m2,deltasR, map, h_values,a1);
    moving(m2,m1,deltasL,map, h_values,a2);
  }
  else if (node21==node22)
  {
    waiting(m2,m1,deltasL,map, h_values,a2);
    moving(m1,m2,deltasR,map, h_values,a1);
  }
  else
  {
    moving(m1,m2,deltasR,map, h_values,a1);
    moving(m2,m1,deltasL,map, h_values,a2);
  }
  if (writeNode)
    output.close();
}

void edgeSpliter::waiting(Move m1, Move m2, Map_deltas &deltas, Map &map, Heuristic &h_values,int a)
{
  int waitingID=m1.id1;
  int node21(m2.id1),node22(m2.id2);
  Vector2D P0(map.get_coord(waitingID));
  Vector2D P1(map.get_coord(node21)),P2(map.get_coord(node22));

  std::vector<Node> succ=map.get_valid_moves(waitingID);
  for (Node n: succ){
    int fromID=n.id;
    if (m1.id1==m2.id2 && fromID==m2.id1) continue;
    assert(n.i!=-1 && n.j!=-1);
    Vector2D nextNode(n.i,n.j);
    Move tempM1(-1,-1,fromID,waitingID);
    
    moving(tempM1,m2,deltas,map, h_values,a);
  }
}

void edgeSpliter::moving(Move m1, Move m2, Map_deltas &deltas, Map &map, Heuristic &h_values,int a)
{
  int node11=m1.id1;
  int node12=m1.id2;
  int node21=m2.id1;
  int node22=m2.id2;
  Vector2D P_new;
  Vector2D P11(map.get_coord(node11)),P12(map.get_coord(node12)),P21(map.get_coord(node21)),P22(map.get_coord(node22));
  double totalT=map.get_dist(node12,node11);
  Vector2D v((P12-P11)/totalT);
  double vx(v.i),vy(v.j);

  if (m2.id1 == m2.id2) //other agent is waiting
  {
    if (m2.id1==m1.id2) //other agent is waiting on the moving agent's destination
    {
      double new_t=round_down(totalT-2*r) - CN_PRECISION;
      P_new=Vector2D(P11+v*(new_t));
    }
    else
    {
      P_new=case2(P11,v, P21);
    }
  }
  else
  {
    //case1
    double x0(P11.i),y0(P11.j),x1(P21.i),y1(P21.j),x2(P22.i),y2(P22.j);
    double dx(x2-x1),dy(y2-y1);
    double D2(dx*dx+dy*dy);

    double a(pow(dx*vy-dy*vx,2));

    if (abs(a)<=CN_EPSILON)
    {
      if(P11.i<=P21.i && P11.i>=P22.i || P11.j <= P21.j && P11.j >= P22.j ||
          P11.i<=P22.i && P11.i>=P21.i || P11.j <= P22.j && P11.j >= P21.j)
        P_new=Vector2D(-1,-1);
      else
        P_new=case2(P11,v,P22);
    }
    else
    {
      double C(4*r*r*D2 - pow(dx*y1-dy*x1,2));
      double b(2*(dx*dx*y0*vy-dx*dx*y1*vy+dy*dy*x0*vx-dy*dy*x1*vx+dx*dy*y1*vx+dx*dy*x1*vy-dx*dy*y0*vx-dx*dy*x0*vy));
      double c( -(C-pow((dx*y0 - dy*x0),2) + 2*dx*dx*y1*y0  + 2*dy*dy*x1*x0 - 2*dx*dy*y1*x0 - 2*dx*dy*x1*y0 ));

      double t2=solveQuad(a,b,c) - CN_PRECISION;

      if (t2<CN_EPSILON)
        P_new=Vector2D(-1,-1);
      else
      {
        P_new=P11+v*t2;

        Vector2D v2(P22-P21);
        double v2x(v2.i),v2y(v2.j);
        if (eq(v2x,0))
        {
          if(ge(P_new.j,P21.j) && ge(P_new.j,P22.j))
            if(gt(P21.j,P22.j))
              if (gt(P21.dis(P_new),2*r))
                P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
          else if(le(P_new.j,P21.j) && le(P_new.j,P22.j))
            if(lt(P21.j,P22.j))
              if (gt(P21.dis(P_new),2*r))
                P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);

        }
        else if(eq(v2y,0))
        {
          if(ge(P_new.i,P21.i) && ge(P_new.i,P22.i))
            if(gt(P21.i,P22.i))
              if (gt(P21.dis(P_new),2*r))
                P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
          else if(le(P_new.i,P21.i) && le(P_new.i,P22.i))
            if(lt(P21.j,P22.i))
              if (gt(P21.dis(P_new),2*r))
                P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
        }
        else if(v2x*v2y>0) // dir:/ checking if two sign are the same
        {
          if(ge(P_new.i,P21.i) && ge(P_new.i,P22.i) && ge(P_new.j,P21.j) && ge(P_new.j,P22.j))
            if(gt(P21.i,P22.i))
              if (gt(P21.dis(P_new),2*r))
                P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
          else if(le(P_new.i,P21.i) && le(P_new.i,P22.i) && le(P_new.j,P21.j) && le(P_new.j,P22.j))
            if(lt(P21.j,P22.i))
              if (gt(P21.dis(P_new),2*r))
                 P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
        }
        else // dir:\  
        {
          if(ge(P_new.i,P21.i) && ge(P_new.i,P22.i) && le(P_new.j,P21.j) && le(P_new.j,P22.j))
            if(gt(P21.i,P22.i))
              if (gt(P21.dis(P_new),2*r))
                P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
          else if(le(P_new.i,P21.i) && le(P_new.i,P22.i) && ge(P_new.j,P21.j) && ge(P_new.j,P22.j))
            if(lt(P21.j,P22.i))
              if (gt(P21.dis(P_new),2*r))
                 P_new=case2(P11,v,P21);
            else
              if (gt(P22.dis(P_new),2*r))
                P_new=case2(P11,v,P22);
        }
        /*
        if (eq(P21.i,P22.i))
        //if (P_new.i==P21.i && P_new==P22.i)
        {
          cout<<"into case2a"<<endl;
          if (gt(P21.j,P22.j))
          {
            if (ge(P_new.j - 2*r,P21.j))
              P_new=case2(P11,v,P21);
            else if(le(P_new.j + 2*r,P22.j))
              P_new=case2(P11,v,P22);
          }
          else 
          {
            if (ge(P_new.j - 2*r,P22.j))
              P_new=case2(P11,v,P22);
            else if(le(P_new.j + 2*r,P21.j))
              P_new=case2(P11,v,P21);
          }

        }
        else if(gt(P21.i,P22.i))
        {
          cout<<"into case2b"<<endl;
          cout<<"P_new:"<<P_new<<" P22"<<P22<<" P21"<<P21<<endl;
          if (gt(P_new.i,P21.i) && gt(P21.dis(P_new),2*r))
            P_new=case2(P11,v,P21);
          else if (lt(P_new.i,P22.i) && gt(P22.dis(P_new),2*r))
            P_new=case2(P11,v,P22);
        }
        else if(lt(P21.i,P22.i))
        {
          cout<<"into case2c"<<endl;
          if (gt(P_new.i,P22.i) && gt(P22.dis(P_new),2*r))
            P_new=case2(P11,v,P22);
          else if (lt(P_new.i,P21.i) && gt(P21.dis(P_new),2*r))
            P_new=case2(P11,v,P21);
        }
      }
      */
    }
  }

  //final adding
  if (P_new.i!=-1 && validNewNode(P11,P12,P_new))
  {
    int new_id=map.add_node(P_new.i,P_new.j,node11,node12);
    if (new_id!=-1)
    {
      if(writeNode)
        output <<new_id <<","<< P_new.i<<","<<P_new.j<<","<<endl;
      Map_delta new_delta(new_id,{node11,node12});
      deltas.push_back(new_delta);
      h_values.add_node(new_id,a,node11);
    }
  }
}

Vector2D edgeSpliter::case2(Vector2D P0,Vector2D v, Vector2D P2)
{
  double x0(P0.i),y0(P0.j);
  double vx(v.i),vy(v.j);
  double x2(P2.i),y2(P2.j);

  double C(4*r*r-x2*x2-y2*y2);
  double a(vx*vx+vy*vy);
  double b(2*(x0*vx - x2*vx + y0*vy - y2*vy));
  double c(-(C-x0*x0 +2*x2*x0-y0*y0+2*y2*y0));

  double t2=solveQuad(a,b,c) - CN_PRECISION;
  Vector2D P(P0+v*t2);
  return P;
}

double edgeSpliter::round_down(double f)
{
  return std::floor(f/CN_EPSILON)*CN_EPSILON;
}


bool edgeSpliter::validNewNode(Vector2D node1,Vector2D node2,Vector2D New)
{
  if ((node1==New) || (node2 == New)) return false;
  if (New.i>node1.i && New.i>node2.i) return false;
  if (New.i<node1.i && New.i<node2.i) return false;
  if (New.j>node1.j && New.j>node2.j) return false;
  if (New.j<node1.j && New.j<node2.j) return false;
  return true;
}

double edgeSpliter::solveQuad(double a, double b, double c)
{
  double delta(b*b-4*a*c);
  assert(delta>-CN_EPSILON);
  if (abs(delta)<CN_EPSILON) delta=0;
  double sqrtDelta(sqrt(delta));
  double t1( round_down( (-b+sqrtDelta)/(2*a) ));
  double t2( round_down( (-b-sqrtDelta)/(2*a) ));
  
  return t2;
}
