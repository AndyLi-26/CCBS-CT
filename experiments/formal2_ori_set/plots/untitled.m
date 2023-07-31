clc; clear all;
%read all data
map=["sparse","dense","super-dense"];
algo=["-0-0-0-0","-0-0-ds-0","-0-ct-ds-0","-0-ct_abs-ds-0","-0-0-ds-icp","-0-ct_abs-ds-icp"];
algo_name=["vanillia","ds","ds+ct","ds+ct_abs","ds+icp","ds+ct_abs+icp"];
mark=["","-o","-x","-+","-*","-hexagram"]
a_size=[0.353553];

y_tot=zeros(38,3);
scatterx_tot=zeros(950,3);
scattery_tot=zeros(950,3);
c_tot=cell(1,3);
tot_stat=[0,0,0];

for m=1:3
    all_data=cell(2,length(algo));
    for alg=1:length(algo)
        T = readtable(strcat(map(m),algo(alg),".csv"));
        T=T{:,:};
        for a=1:1
            temp=T(T(:,2)==a_size(a),:);
            all_data{a,alg}=temp; 
        end
    end
    %plot succ rate
    figure('Name',strcat(map(m),' succ rate')); hold on;
    title("succ rate")
    for a=1:1
        subplot(2,1,a); hold on; title(sprintf("agent size=%f",a_size(a)));
        for alg=1:length(algo)
            data=all_data{a,alg};
            x=unique(data(:,1));
            tsum = accumarray(data(:,1),data(:,5));
            y = tsum(x(:));
            if a==1 && (alg==1 || alg==2 || alg==6)
                switch alg
                    case 1
                        y_tot(:,1)=y_tot(:,1)+y;
                    case 2
                        y_tot(:,2)=y_tot(:,2)+y;
                    case 6   
                        y_tot(:,3)=y_tot(:,3)+y;
                end
            end
            plot(x,y,mark(alg));
        end
        legend(algo_name);
    end
    
    %plot node expansion
    comparsion=[2,3;2,4 ;2,5;2,6];
    for a=1:1
        figure('Name',strcat(map(m),' Node expansion, r=',num2str(a_size(a))));
        for i=1:size(comparsion,1)
            subplot(2,2,i); hold on; 
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
            c(check3(:,1),:)=repmat([1,1,1],sum(check3),1);
    
    
            scatter(datax(:,9),datay(:,9),30,c,'filled');
            if i==4
                scatterx_tot(:,m)=datax(:,9);
                scattery_tot(:,m)=datay(:,9);
                c_tot{m}=c;
                tot_stat(1)=tot_stat(1)+sum(check1 & ~check3);
                tot_stat(2)=tot_stat(2)+sum(eqcheck & ~check3);
                tot_stat(3)=tot_stat(3)+sum(check2 & ~check3);
            end
            title(sprintf("%d is better",sum(check2 & ~check3)));
            set(gca,'xscale','log');set(gca,'yscale','log')
            tempm=max([maxNodex,maxNodey]);
            tempx=linspace(1,tempm,1000);
            plot(tempx,tempx,"r-");hold on;
            xlabel(algo_name(comparsion(i,1)));
            ylabel(algo_name(comparsion(i,2))); 
        end
    end

end


figure(); hold on; title("total succ rate")
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