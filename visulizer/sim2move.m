clc;clear all; close all;
global resolution
pos11=[96.9881,165.795];
pos12=[153.478, 157.941];
pos21=[94,175.1];
pos22=[100.3,173.3];
t1=99.9188;t2=99.9307;
agent_size=4.5;
resolution =0.01;

%rang=abs([pos11(1)-pos12(1),pos11(2)-pos12(2),pos21(1)-pos22(1),pos21(2)-pos22(2)]);
%temprang=max(rang);
%rang=1.2*temprang;
%minx=min([pos11(1),pos12(1),pos21(1),pos22(1)]);
%miny=min([pos11(2),pos12(2),pos21(2),pos22(2)]);
figure();
xlim([85,120]);
ylim([145,180]);
%xlim([minx-0.5*temprang,minx+rang]);ylim([miny-0.5*temprang,miny+rang]);

base=min(t1,t2);
t1=t1-base;t2=t2-base;
d1=norm(pos11-pos12);d2=norm(pos21-pos22);
info1=[0,0,0,0,t1+d1;pos11,pos11,t1;
    pos11,pos12,d1];

info2=[0,0,0,0,t2+d2;pos21,pos21,t2;
    pos21,pos22,d2];
M1=genPathMatrix(info1);
M2=genPathMatrix(info2);
colors=[1,0,0;0,0,1];
plots =plot([],[]);agentNum =plot([],[]);
for t=1:max(t1+d1,t2+d2)/resolution
    title(sprintf("t=%.1f",t*resolution))
    delete(plots);delete(agentNum);
    plots=[];agentNum=[];
    
    p=M1;
    t;
    if size(p,1)>t
        p1=[p(t,1),p(t,2)];
        plots=[plots;circle(p(t,1),p(t,2),agent_size,colors(1,:))];
        agentNum=[agentNum,text(p(t,1),p(t,2),"0",'Color',[1,1,1]-colors(1,:),'FontSize',22,'FontWeight','bold')];
    else
        p1=[p(end,1),p(end,2)];
        plots=[plots;circle(p(end,1),p(end,2),agent_size,colors(1,:))];
        agentNum=[agentNum,text(p(end,1),p(end,2),"0",'Color',[1,1,1]-colors(1,:),'FontSize',22,'FontWeight','bold')];
    end


    p=M2;
    t;
    if size(p,1)>t
        p2=[p(t,1),p(t,2)];
        plots=[plots;circle(p(t,1),p(t,2),agent_size,colors(2,:))];
        agentNum=[agentNum,text(p(t,1),p(t,2),"1",'Color',[1,1,1]-colors(2,:),'FontSize',22,'FontWeight','bold')];
    else
        p2=[p(end,1),p(end,2)];
        plots=[plots;circle(p(end,1),p(end,2),agent_size,colors(2,:))];
        agentNum=[agentNum,text(p(end,1),p(end,2),"1",'Color',[1,1,1]-colors(2,:),'FontSize',22,'FontWeight','bold')];
    end
    
    title(sprintf("t=%.1f,dis=%.3f",t*resolution,norm(p1-p2)-2*agent_size))
    if (norm(p1-p2)-2*agent_size)<0
        t*resolution+base
    end
    drawnow;
    grid off;
    %exportgraphics(gcf,'ccbs_unsolve.gif','Append',true);
end








function M=genPathMatrix(info)
    global resolution
    M=zeros(ceil(info(1,5)/resolution),2);
    Y=[0;ceil(cumsum(info(2:end,5))/resolution)];
    for i=2:size(info,1)
        x=linspace(info(i,1),info(i,3),Y(i)-Y(i-1));
        y=linspace(info(i,2),info(i,4),Y(i)-Y(i-1));
        M(Y(i-1)+1:Y(i),1)=x;
        M(Y(i-1)+1:Y(i),2)=y;
    end
end

function drawline(s,t,spec)
    global nodes
    xs=[nodes(s,1),nodes(t,1)];ys=[nodes(s,2),nodes(t,2)];
    plot(xs,ys,spec,LineWidth=2);
end

function out = circle(x0,y0,radius,color)
pos = [[x0,y0]-radius 2*radius 2*radius];
out = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', color, 'Edgecolor','none');
end