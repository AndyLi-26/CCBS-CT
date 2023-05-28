clc;clear all;close all;
a1=[5,5;25,25]; a2=[5,40;20,5];
r=4.5;

plot(a1(:,1),a1(:,2));
hold on;plot(a2(:,1),a2(:,2));
xlim([0,35]);ylim([0,35]);

x_0=a1(1,1);y_0=a1(1,2);
x_1=a2(1,1);y_1=a2(1,2);
x_2=a2(2,1);y_2=a2(2,2);

v=(a1(2,:)-a1(1,:));
v=v/norm(v);
v_x=v(1);v_y=v(2);

d_x=x_2-x_1;d_y=y_2-y_1;
D2=d_x^2+d_y^2;
C=4*r^2*D2-(d_x*y_1-d_y*x_1)^2;

a=(d_x*v_y-d_y*v_x)^2;
b=2*(d_x^2*y_0*v_y-d_x^2*y_1*v_y+d_y^2*x_0*v_x-d_y^2*x_1*v_x+d_x*d_y*y_1*v_x+d_x*d_y*x_1*v_y-d_x*d_y*y_0*v_x-d_x*d_y*x_0*v_y);
c=-(C - (d_x*y_0 - d_y*x_0)^2 + 2*d_x^2*y_1*y_0  + 2*d_y^2*x_1*x_0 - 2*d_x*d_y*y_1*x_0 - 2*d_x*d_y*x_1*y_0 );
delta=b^2-4*a*c
t1=(-b+sqrt(delta))/(2*a)
t2=(-b-sqrt(delta))/(2*a)

t=0:0.1:20;

xa=x_0+v_x*t2;
ya=y_0+v_y*t2;
circle(xa,ya,r,[0,0,1])
t=0:0.1:norm(a2(2,:)-a2(1,:));
a2x=linspace(a2(1,1),a2(2,1),length(t));
a2y=linspace(a2(1,2),a2(2,2),length(t));
plots=[];
for i=1:length(t)
    delete(plots);
    plots=[];
    plots=[circle(a2x(i),a2y(i),r,[0,1,0])];
    drawnow;
end
%d_der=d_x^2*(y_0^2+v_y^2*t.^2+2*y_0*v_y*t) ...
%    -2*d_x^2*y_1*(y_0+v_y*t) ...
%    + d_y^2*(x_0^2+v_x^2*t.^2+2*x_0*v_x*t) ...
%    -2*d_y^2*x_1*(x_0+v_x*t) ...
%     + 2 * d_x*d_y*y_1*(x_0+v_x*t) ...
%     + 2 * d_x*d_y*x_1*(y_0+v_y*t) ...
%     - 2 * d_x*d_y*(y_0+v_y*t).*(x_0+v_x*t) ...
%    -C;

d_der=(d_x*v_y - d_y*v_x)^2*t.^2 ...
     + 2*(d_x^2*y_0*v_y - d_x^2*y_1*v_y + d_y^2*x_0*v_x - d_y^2*x_1*v_x ...
     + d_x*d_y*y_1*v_x + d_x*d_y*x_1*v_y - d_x*d_y*y_0*v_x - d_x*d_y*x_0*v_y)*t ...
    -(C - (d_x*y_0 - d_y*x_0)^2 ...
    + 2*d_x^2*y_1*y_0  + 2*d_y^2*x_1*x_0 - 2*d_x*d_y*y_1*x_0 ...
    - 2*d_x*d_y*x_1*y_0 );

d_ori=(d_x*(y_1-ya)-(x_1-xa)*d_y).^2-4*r^2*D2;
figure
hold on;
plot(t,d_ori)
plot(t,d_der)
legend("ori","der")


function out = circle(x0,y0,radius,color)
pos = [[x0,y0]-radius 2*radius 2*radius];
out = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', color, 'Edgecolor','none');
end