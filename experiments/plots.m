clc; clear all; close all;
%read all data
algo=["ORI","my"];
map=["sparse","dense","super-dense"];
all_data=cell(2,3);
a_num=2:34;

%plot succ rate
figure; hold on;
title("succ rate")
for i=1:3
    subplot(3,1,i); hold on; title(map(i));
    T=readtable(strcat(map(i),"_ORI",".csv"));
    T=T{:,:};
    plot(a_num,T(1,:),'b-','LineWidth',2);
    plot(a_num,T(2,:),'rx','LineWidth',3);

    T=readtable(strcat(map(i),"_my",".csv"));
    T=T{:,:};
    plot(a_num,T(1,:),'g-','LineWidth',2);
    legend("their version","their version stuck in inf loop","our version");
end

figure;hold on;
T=readtable("lower Bound.csv");
T=T{:,:};
x=T(1,:);y=T(2,:);
plot(x,y,"o");
set(gca,'xscale','log');set(gca,'yscale','log')