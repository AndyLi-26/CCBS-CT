clc;clear all;close all;
for i=0:4
    figure('Name',sprintf("r=%i.5",i));
    f=strcat('r_',int2str(i));
    new=readmatrix(f);
    y=new(:,6:27);
    x=5:26;
    subplot(3,1,1);hold on;
    xlabel("agent number");ylabel("Succ rate");
    plot(x,[y(1,:);y(4,:);y(7,:);y(10,:)])
    legend("Vanilla","DS","CR","DS+CR");
    subplot(3,1,2);hold on;
    semilogy(x,[y(2,:);y(5,:);y(8,:);y(11,:)])
    legend("Vanilla","DS","CR","DS+CR");
    xlabel("agent number");ylabel("Node expanded");
    subplot(3,1,3);hold on;
    semilogy(x,[y(3,:);y(6,:);y(9,:);y(12,:)])
    legend("Vanilla","DS","CR","DS+CR");
    xlabel("agent number");ylabel("lower bound");
end