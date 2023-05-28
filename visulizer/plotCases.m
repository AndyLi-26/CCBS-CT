clc; clear all; close all;
r=4.5;range=[-5,60];
%a1=[5,5;40,40]; a2=[40,5;5,40]; %case1
%a1=[40,5;40,40]; a2=[5,5;40,40]; %case2
%a1=[40,5;40,40]; a2=[30,10;37,30]; %case3
%a1=[40,5;40,40]; a2=[40,40;40,40]; %case4
%a1=[40,5;40,30]; a2=[20,50;37,36]; %case5
%a1=[10,5;40,40]; a2=[30,18;30,18]; %case6
%a1=[10,5;40,40]; a2=[25,20;30,30]; %case7
%a1=[10,5;40,40]; a2=[30,30;25,20]; %case8
%a1=[30,5;30,30]; a2=[5,40;30,30]; %case9
%a1=[30,5;30,30]; a2=[30,40;30,20]; %case10
%a1=[5,5;20,20]; a2=[40,35;15,10]; %case11
%a1=[5,20;20,20]; a2=[40,5;20,20]; %case12
%a1=[30,5;30,30]; a2=[40,40;30,30]; %case13
%a1=[10,15;30,35]; a2=[40,40;10,10]; %case14
newNode=calcNewP(a1,a2,r)

figure; hold on; xlim(range);ylim(range);
%draw a1
circle(a1(1,:),r,[0,0,1],'-');circle(a1(2,:),r,[0,0,1],'-'); 
anArrow = annotation('arrow','Color',[0,0,1]) ;
anArrow.Parent = gca;
anArrow.X = a1(:,1);
anArrow.Y = a1(:,2);

%draw a2
circle(a2(1,:),r,[0,1,0],'-');
if ~isequal( a2(1,:),a2(2,:))
    circle(a2(2,:),r,[0,1,0],'-');
    anArrow = annotation('arrow','Color',[0,1,0]) ;
    anArrow.Parent = gca;
    anArrow.X = a2(:,1);
    anArrow.Y = a2(:,2);
    %draw boarder
    v2=a2(2,:)-a2(1,:);
    v2=v2/norm(v2);
    vp=[-v2(2),v2(1)];
    vp=r*[vp;vp];
    b1=a2-vp;
    plot(b1(:,1),b1(:,2),'g--');
    b2=a2+vp;
    plot(b2(:,1),b2(:,2),'g--');
end
%draw new point
circle(newNode,r,[1,0,0],'--');
%exportgraphics(gca,'sparse.png','Resolution',300)




function out = circle(P,radius,color,Lstyle)
pos = [P-radius 2*radius 2*radius];
out = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'none', 'Edgecolor',color,'LineStyle',Lstyle);
end