clc; clear all; close all;
recording=true;
if recording    
obj = VideoWriter("zoomin.mp4", 'MPEG-4');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj);
end
agent_size=0.5;
%draw mapi
global nodes resolution
resolution=0.01;
nodes=readmatrix('nodes.csv');

mapSize=[min(nodes,[],'all'),max(nodes,[],'all')];range=mapSize(2)-mapSize(1);
mapSize=[mapSize(1)-0.01*range,mapSize(2)+0.01*range];
mapSize=[5,300];
edges=readmatrix('edges.csv');
%summary=readmatrix(' summary.txt');%run-time, SoC, makespan

figure('Renderer', 'painters', 'Position', [10 10 900 900]);
hold on; grid on; xlim(mapSize);ylim(mapSize);
%xlim([150,250]);ylim([130,230]);
%xlabel(sprintf("run-time=%f,    SoC=%f,    makespan=%f",summary));
plot(nodes(:,1),nodes(:,2),'ro',LineWidth=2,MarkerSize=1)
for i=1:length(nodes)
    i
    text(nodes(i,1),nodes(i,2),num2str(i-1),'FontSize',8,'FontWeight','bold');
end

% new=readmatrix('newNode.csv');
% plot(new(:,2),new(:,3),'bo',LineWidth=2)
% for i=1:size(new,1)
%     i
%     text(new(i,2),new(i,3),num2str(new(i,1)),'FontSize',20,'FontWeight','bold','Color','blue');
% end


dcm = datacursormode;          
dcm.Enable = 'on';
dcm.UpdateFcn = @displayind;
for i=1:size(edges,1)
    drawline(edges(i,1)+1,edges(i,2)+1,'r-')
end
%xlim([39,311]);ylim([39,611])
%I = imread('warehouseMap.png');
%h = image(xlim,ylim,I);
%uistack(h,'bottom')
%xlim([30,320]);ylim([30,620])
%% plot tasks
tasks=readmatrix('tasks.csv');
%agents=size(tasks,1);
a=dir(['paths/*.csv']);
agents=length(extractfield(a,"name"));
colors=jet(agents+1);
%colors=double([0xb6,0xd7,0xa8;0xa4,0xc2,0xf4])/255;
for i=1:agents
    plot(nodes(tasks(i,1)+1,1),nodes(tasks(i,1)+1,2),'o',MarkerFaceColor=colors(i,:),MarkerSize=10);
    text(nodes(tasks(i,1)+1,1),nodes(tasks(i,1)+1,2),num2str(i),'Color','black','FontSize',9)
end
for i=1:agents
    plot(nodes(tasks(i,2)+1,1),nodes(tasks(i,2)+1,2),'p',MarkerEdgeColor=colors(i,:),MarkerSize=8,LineWidth=3);
end
if recording
    f = getframe(gcf);
    writeVideo(obj,f);
end
%gen animation
paths=cell(1,agents);
tmax=-1;
for i=1:agents
    info=readmatrix(['paths/',num2str(i-1),'.csv']);
    if info(1,5)>tmax
        tmax=info(1,5);
    end
    paths(i)={genPathMatrix(info)};
end

plots =plot([],[]);agentNum =plot([],[]);
for t=173/resolution:180/resolution %tmax
    title(sprintf("small ds t=%.1f",t*resolution))
    t*resolution
    if (t*resolution>155)
        xlim([185,195]);ylim([168,178]);
    end
    delete(plots);delete(agentNum);
    plots=[];agentNum=[];
    checka1=0;checka2=3;
    for a=1:agents
        p=paths{a};
        if a==checka1+1
            tempp1=p(t,:);
        end
        if a==checka2+1
            tempp2=p(t,:);
        end
        if size(p,1)>t
            plots=[plots;circle(p(t,1),p(t,2),agent_size,colors(a,:))];
            agentNum=[agentNum,text(p(t,1),p(t,2),num2str(a-1),'Color',[1,1,1]-colors(a,:),'FontSize',22,'FontWeight','bold')];
        else
            plots=[plots;circle(p(end,1),p(end,2),agent_size,colors(a,:))];
            agentNum=[agentNum,text(p(end,1),p(end,2),num2str(a-1),'Color',[1,1,1]-colors(a,:),'FontSize',22,'FontWeight','bold')];
        end
    end
    drawnow;
    norm(tempp1-tempp2)
    grid off;
    if recording
        f = getframe(gcf);
        writeVideo(obj,f);
    end
    %exportgraphics(gcf,'ccbs_unsolve.gif','Append',true);
end
if recording
    f = getframe(gcf);
    %exportgraphics(gcf,'paper.gif','Append',true);
    obj.close();
end
function out = circle(x0,y0,radius,color)
pos = [[x0,y0]-radius 2*radius 2*radius];
out = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', color, 'Edgecolor','none');
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
    plot(xs,ys,spec,LineWidth=1);
end

function txt=displayind(~,info)
    global nodes
    x = info.Position(1);
    y = info.Position(2);
    txt=num2str(find(ismember(nodes, [x y],'rows'))-1);
end

