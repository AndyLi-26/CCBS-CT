clc; clear all; close all;
maps=["sparse","dense","super-dense"];
for m=1:3
    figure; hold on; title(maps(m));
    T = readtable(strcat(maps(m),".csv"));
    T=T{:,:};
    for ds=[1,0]
        for es=[1,0]
            for cr=[1,0]
                temp=T(T(:,4)==ds,:);
                temp=temp(temp(:,5)==es,:);
                temp=temp(temp(:,6)==cr,:);
                temp=temp(temp(:,2)==4.5,:);
                succ=zeros(1,21);
                for a=[5,10,15,25]
                    aaa=temp(temp(:,1)==a,:);
                    if (size(aaa,1)==0)
                        succ(a-4)=0;
                    else
                        succ(a-4)=sum(aaa(:,8));
                    end
                end
                plot(5:25,succ);
            end
        end
    end
    legend("DS+ES+CR","DS+ES","DS+CR","DS","ES+CR","ES","CR");
%     run_t=zeros(1,17);suc=zeros(1,17);fail=zeros(1,17);
%     for i=1:length(T.Var1)
%         if T.Var2(i)<300
%             run_t(T.Var1(i)/2)=run_t(T.Var1(i)/2)+T.Var2(i);
%             suc(T.Var1(i)/2)=suc(T.Var1(i)/2)+1;
%         else
%             fail(T.Var1(i)/2)=fail(T.Var1(i)/2)+1;
%         end
%     end
%     figure;
%     run_t=run_t./suc;
%     plot(2:2:34,suc/25);
%     xticks(2:2:34)
%     title(maps(m))
%     xlabel("number of agents");
%     ylabel("Success rate");
end
