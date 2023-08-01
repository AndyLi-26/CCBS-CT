clc; clear all; close all;
%read all data
algo=["sparse-0-0-0-0","sparse-0-0-ds-0","sparse-0-0-0-icp","sparse-0-0-ds-icp"];
algo_name=["vanillia","ds","icp","ds+icp"];
all_data=cell(2,4);
a_size=[0.5,4.5];
for alg=1:4
    T = readtable(strcat(algo(alg),".csv"));
    T=T{:,:};
    for a=1:2
        temp=T(T(:,2)==a_size(a),:);
        all_data{a,alg}=temp; 
    end
end
%plot succ rate
figure; hold on;
title("succ rate")
for a=1:2
    subplot(2,1,a); hold on; title(sprintf("agent size=%f",a_size(a)));
    for algo=1:4
        data=all_data{a,algo};
        x=unique(data(:,1));
        tsum = accumarray(data(:,1),data(:,5));
        y = tsum(x(:))/25*100;
        plot(x,y);
    end
    legend(algo_name);
end

%plot node expansion
comparsion=[1,3;1,4;2,3;2,4];
for a=1:2
    figure('Name','Node expansion')
    for i=1:size(comparsion,1)
        subplot(2,2,i); hold on; title(sprintf("agent size=%f",a_size(a)));
        datax=all_data{a,comparsion(i,1)};
        datay=all_data{a,comparsion(i,2)};
        maxNodex=max(datax(:,9));maxNodey=max(datay(:,9));
        datax(datax(:,5)==0,:)=maxNodex;
        datay(datay(:,5)==0,:)=maxNodey;
        check1=datax(:,9)<datay(:,9);
        check2=datax(:,9)>datay(:,9);
        c=zeros(size(datax,1),3);
        c(:,3)=1;
        c(check1(:,1),:)=repmat([1,0,0],sum(check1),1);
        c(check2(:,1),:)=repmat([0,1,0],sum(check2),1);

        scatter(datax(:,9),datay(:,9),30,c,'filled');
        
        set(gca,'xscale','log');set(gca,'yscale','log')
        m=max([maxNodex,maxNodey]);
        tempx=linspace(1,m,1000);
        plot(tempx,tempx,"r-");hold on;
        xlabel(algo_name(comparsion(i,1)));
        ylabel(algo_name(comparsion(i,2))); 
    end
end