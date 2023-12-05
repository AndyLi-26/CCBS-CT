clc; clear all; close all;
%read all data
map=["sparse","dense","super-dense"];
a_num=[37,38,32];
algo=["-0-0-0-0","-0-0-ds-0","-0-ct-ds-0","-0-ct_abs-ds-0","-0-0-ds-icp","-0-ct_abs-ds-icp"];
algo_name=["vanillia","ds","ds+ct","ds+ct-abs","ds+icp","ds+ct-abs+icp"];
%a_size=[0.353553385,0.5,1,4.5];
a_size=[0.5];
%y_tot=zeros(38,3);
%scatterx_tot=zeros(950,3);
%scattery_tot=zeros(950,3);
%c_tot=cell(1,3);
%tot_stat=[0,0,0];
figure(1);figure(2);
mark=["-o","-x","-+","-hexagram"]
for m=1:3
    all_data=cell(2,1);
    for a=1:length(algo)
        T = readtable(strcat(map(m),algo(i),".csv"));
    T=T{:,:};
        temp=T(T(:,2)==a_size(a),:);
        all_data{a,1}=temp; 
    end
    %plot succ rate
%     figure(1); 
%     subplot(3,1,m);
%     hold on;
% 
%     title(strcat(map(m),' succ rate'));
%     
%     for a=1:length(a_size)
%         data=all_data{a,1};
%         x=unique(data(:,1));
%         tsum = accumarray(data(:,1),data(:,5));
%         y = tsum(x(:));
%         %y_tot(:,1)=y_tot(:,1)+tsum(x(:));
%         x=x(1:a_num(m)-2);y=y(1:a_num(m)-2);
%         plot(x,y,mark(a),'MarkerSize',12,'LineWidth', 3);
%     end
%     legend(["r=0.353553385","r=0.5","r=1","r=4.5"]);
%     if m==3
%         xlabel("agents");
%     end
%     ylabel("solved instance");
%     xticks(x);
    %saveas(gcf,strcat(map(m),' succ rate r=',num2str(a_size(a)),".png"))
    figure(2);
    %plot node expansion
    comparsion=[2,3;2,4;2,6];
    for a=1:length(a_size)
        for i=1:size(comparsion,1)
            subplot(3,3,i+(m-1)*3); hold on; 
            datax=all_data{a,comparsion(i,1)};
            datay=all_data{a,comparsion(i,2)};
            maxNodex=max(datax(:,9));maxNodey=max(datay(:,9));
            datax(datax(:,5)==0,9)=maxNodex;
            datay(datay(:,5)==0,9)=maxNodey;
            check1=datax(:,9)<datay(:,9);
            check2=datax(:,9)>datay(:,9);
            eqcheck=datax(:,9)==datay(:,9);
            check3=datax(:,5)==0 & datay(:,5)==0;
            c=zeros(size(datax,1),3);
            c(:,3)=1;
            c(check1(:,1),:)=repmat([1,0,0],sum(check1),1);
            c(check2(:,1),:)=repmat([0,1,0],sum(check2),1);
            c(check3(:,1),:)=repmat([0,0,0],sum(check3),1);
    
    
            scatter(datax(:,9),datay(:,9),30,c,'filled');
            if i==4
                scatterx_tot(:,m)=datax(:,9);
                scattery_tot(:,m)=datay(:,9);
                c_tot{m}=c
                tot_stat(1)=tot_stat(1)+sum(check1 & ~check3);
                tot_stat(2)=tot_stat(2)+sum(eqcheck & ~check3);
                tot_stat(3)=tot_stat(3)+sum(check2 & ~check3);
            end
            title(sprintf("↑ %d, = %d, ↓ %d",sum(check2 & ~check3), sum(eqcheck & ~check3), sum(check1 & ~check3)));
            set(gca,'xscale','log');set(gca,'yscale','log');
            tempm=max([maxNodex,maxNodey]);
            tempx=linspace(1,tempm,1000);
            plot(tempx,tempx,"r-");hold on;
            xlabel(algo_name(comparsion(i,1)));
            ylabel(algo_name(comparsion(i,2))); 
        end
        %saveas(gcf,strcat(map(m),' Node expension r=',num2str(a_size(a)),'.png'))
    end

end
figure(1);set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'success rate.png','Resolution',300)
fontsize(gcf, 18, "points")
figure(2);set(gcf, 'Position', get(0, 'Screensize'));
fontsize(gcf, 24, "points")
exportgraphics(gcf,'node expansion.png','Resolution',300)
%%
figure; hold on; title("total succ rate")
plot(x,y_tot(:,1),'b-o');
plot(x,y_tot(:,2),'r-s');
plot(x,y_tot(:,3),"g-diamond");
legend(["vanillia","DS","DS+CT"]);
figure();hold on; title("")
title(sprintf("total node expanded:\n %d is worse, %d is same, %d is better", tot_stat));
for i=1:3
    scatter(scatterx_tot(:,i),scattery_tot(:,i),30,c_tot{i},'filled');
end
set(gca,'xscale','log');set(gca,'yscale','log');
tempm=max([max(max(scatterx_tot)),max(max(scattery_tot))]);
plot(linspace(1,tempm,1000),linspace(1,tempm,1000),'r-');