clc; clear all; close all;
%read all data
map=["sparse","dense","super-dense"];
a_num=[37,38,22];
algo=["-0-0-0-0","-0-0-ds-0","-0-ct_abs-ds-0","-0-icp-ds-0","-0-ct_abs-ds-2","-0-icp-ds-2"];
algo_name=["vanillia","ds","ds+ct","ds+icp","ds+ct-2","ds+icp+2"];
a_size=[0.5];
figure(1);figure(2);
for m=1:3
    all_data=cell(2,length(algo));
    for alg=1:length(algo)
        T = readtable(strcat(map(m),algo(alg),".csv"));
        T=T{:,:};
        for a=1:length(a_size)
            temp=T(T(:,2)==a_size(a),:);
            all_data{a,alg}=temp; 
        end
    end
  
    figure(1);
    %plot node expansion
    comparsion=[2,3;2,4;];
    for a=1:length(a_size)
        for i=1:size(comparsion,1)
            subplot(3,length(comparsion),i+(m-1)*length(comparsion)); hold on; 
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
fontsize(gcf, 24, "points")
exportgraphics(gcf,'node expansion.png','Resolution',300)