clc; clear all; close all;
%read all data
algo=["super-dense-0-0-0-0","super-dense-0-0-ds-0","super-dense-0-ct-0-0","super-dense-0-ct-ds-0","super-dense-0-ct_abs-0-0","super-dense-0-ct_abs-ds-0"];
algo_name=["vanillia","ds","ct","ds+ct","ct_abs","ds+ct_abs"];
all_data=cell(2,6);
a_size=[0.5,4.5];
for alg=1:6
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
    for algo=1:6
        data=all_data{a,algo};
        x=unique(data(:,1));
        tsum = accumarray(data(:,1),data(:,5));
        y = tsum(x(:))/25*100;
        plot(x,y);
    end
    legend(algo_name);
end

%plot node expansion
comparsion=[1,5;1,6;2,5;2,6];
for a=1:2
    figure('Name','Node expansion')
    for i=1:size(comparsion,1)
        subplot(2,2,i); hold on; 
        datax=all_data{a,comparsion(i,1)};
        datay=all_data{a,comparsion(i,2)};
        maxNodex=max(datax(:,9));maxNodey=max(datay(:,9));
        datax(datax(:,5)==0,9)=maxNodex;
        datay(datay(:,5)==0,9)=maxNodey;
        check1=datax(:,9)<datay(:,9);
        check2=datax(:,9)>datay(:,9);
        check3=datax(:,5)==0 & datay(:,5)==0;
        c=zeros(size(datax,1),3);
        c(:,3)=1;
        c(check1(:,1),:)=repmat([1,0,0],sum(check1),1);
        c(check2(:,1),:)=repmat([0,1,0],sum(check2),1);
        c(check3(:,1),:)=repmat([1,1,1],sum(check3),1);


        scatter(datax(:,9),datay(:,9),30,c,'filled');
        title(sprintf("agent size=%f, better with%d",a_size(a),sum(check2 & ~check3)));
        set(gca,'xscale','log');set(gca,'yscale','log')
        m=max([maxNodex,maxNodey]);
        tempx=linspace(1,m,1000);
        plot(tempx,tempx,"r-");hold on;
        xlabel(algo_name(comparsion(i,1)));
        ylabel(algo_name(comparsion(i,2))); 
    end
end