clc; clear all; close all;
%read all data
map=["sparse","dense","super-dense"];
a_num=[37,38,22];
algo=["-0-0-0","-0-0-ds","-0-ct_abs-ds","-0-icp-ds","-0-ct_abs-0","-0-icp-0"];
algo_name=["vanillia","ds","ds+ct","ds+icp","ct","icp"];
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
%     plot succ rate
    figure(1); 
    subplot(3,1,m);
    x=0:0.01:30;
    hold on;
    toPlot=[1 2 3 4];
    title(strcat(map(m),' solved instance'));
    
    for p=1:length(toPlot)
        y=zeros(size(x));
        data=all_data{a,p};
        for i=1:length(x)
            y(i)=sum(1);
        
        end
    
    end
    data=all_data{a,1};
    x=unique(data(:,1));
    tsum = accumarray(data(:,1),data(:,5));
    y = tsum(x(:));

    x=x(1:a_num(m)-2);y=y(1:a_num(m)-2);
    plot(x,y,'-o','MarkerSize',12,'LineWidth', 3);

    data=all_data{a,2};
    x=unique(data(:,1));
    tsum = accumarray(data(:,1),data(:,5));
    y = tsum(x(:));

    x=x(1:a_num(m)-2);y=y(1:a_num(m)-2);
    plot(x,y,'-s','MarkerSize',12,'LineWidth', 3);

    data=all_data{a,3};
    x=unique(data(:,1));
    tsum = accumarray(data(:,1),data(:,5));
    y = tsum(x(:));
    x=x(1:a_num(m)-2);y=y(1:a_num(m)-2);
    plot(x,y,"-x",'MarkerSize',12,'LineWidth', 3);

    data=all_data{a,4};
    x=unique(data(:,1));
    tsum = accumarray(data(:,1),data(:,5));
    y = tsum(x(:));
    x=x(1:a_num(m)-2);y=y(1:a_num(m)-2);
    plot(x,y,"-diamond",'MarkerSize',12,'LineWidth', 3);

%     data=all_data{a,6};
%     x=unique(data(:,1));
%     tsum = accumarray(data(:,1),data(:,5));
%     y = tsum(x(:));
% 
%     x=x(1:a_num(m)-2);y=y(1:a_num(m)-2);
%     plot(x,y,"-*",'MarkerSize',12,'LineWidth', 3);

    legend(["vanillia","DS","DS+CT-abs","DS+CT-abs+ICP"]);
    if m==3
        xlabel("agents");
    end
    ylabel("solved instance");
    xticks(x);
    %saveas(gcf,strcat(map(m),' succ rate r=',num2str(a_size(a)),".png"))
    figure(2);
    %plot node expansion
    comparsion=[2,3;2,4;1,4];
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