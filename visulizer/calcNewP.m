function P=calcNewP(a1,a2,r)
    x_0=a1(1,1);y_0=a1(1,2);
    x_1=a2(1,1);y_1=a2(1,2);
    x_2=a2(2,1);y_2=a2(2,2);
    
    v=(a1(2,:)-a1(1,:));
    totalT=norm(v);
    v=v/totalT;
    v_x=v(1);v_y=v(2);
    
    if a2(1,:)==a2(2,:)
        if a2(1,:)==a1(2,:) %case3
            fprintf("case3")
            P=a1(1,:)+v*(totalT-2*r);
        else %case2
            fprintf("case2")
            P=case2(r,a2(1,:),v,a1(1,:));
        end
        return
    end
    
    
    d_x=x_2-x_1;d_y=y_2-y_1;
    D2=d_x^2+d_y^2;
    C=4*r^2*D2-(d_x*y_1-d_y*x_1)^2;
    
    a=(d_x*v_y-d_y*v_x)^2;
    if a==0
        if ( a1(1,1)<=a2(1,1) && a1(1,1) >= a2(2,1) ) || (a1(1,2)<=a2(1,2) && a1(1,2) >= a2(2,2)) || ( a1(1,1)<=a2(2,1) && a1(1,1) >= a2(1,1) ) || (a1(1,2)<=a2(2,2) && a1(1,2) >= a2(1,2))
            P=[-1,-1];
        else
            P=case2(r,a2(2,:),v,a1(1,:));
        end
        return
    end
    
    b=2*(d_x^2*y_0*v_y-d_x^2*y_1*v_y+d_y^2*x_0*v_x-d_y^2*x_1*v_x+d_x*d_y*y_1*v_x+d_x*d_y*x_1*v_y-d_x*d_y*y_0*v_x-d_x*d_y*x_0*v_y);
    c=-(C - (d_x*y_0 - d_y*x_0)^2 + 2*d_x^2*y_1*y_0  + 2*d_y^2*x_1*x_0 - 2*d_x*d_y*y_1*x_0 - 2*d_x*d_y*x_1*y_0 );
    delta=b^2-4*a*c;
    t1=(-b+sqrt(delta))/(2*a)
    t2=(-b-sqrt(delta))/(2*a)
    P=a1(1,:)+v*t2;
    if P(1)==a2(1,1) && P(1)==a2(2,1)
        fprintf("case2")
        if P(2)>=a2(1,2) && P(2)>=a2(2,2)
            if a2(2,2)>a2(1,2)
                P=case2(r,a2(2,:),v,a1(1,:));
            else
                P=case2(r,a2(1,:),v,a1(1,:));
            end
        elseif P(2)<=a2(1,2) && P(2)<=a2(2,2)
            if a2(2,2)<a2(1,2)
                P=case2(r,a2(2,:),v,a1(1,:));
            else
                P=case2(r,a2(1,:),v,a1(1,:));
            end
        end
    
    elseif P(1)>=a2(1,1) && P(1)>=a2(2,1)
        fprintf("case2")
        if a2(1,1) > a2(2,1)
            P=case2(r,a2(1,:),v,a1(1,:));
        else
            P=case2(r,a2(2,:),v,a1(1,:));
        end
    elseif P(1)<=a2(1,1) && P(1)<=a2(2,1)
        fprintf("case2")
        if a2(1,1) < a2(2,1)
            P=case2(r,a2(1,:),v,a1(1,:));
        else
            P=case2(r,a2(2,:),v,a1(1,:));
        end
    
    end
end


function P=case2(r,P2,v,P0)
    x_0=P0(1);y_0=P0(2);
    v_x=v(1);v_y=v(2);
    x_2=P2(1);y_2=P2(2);

    C=4*r^2-x_2^2-y_2^2;
    a=v_x^2+v_y^2;
    b=2*(x_0*v_x - x_2*v_x + y_0*v_y - y_2*v_y);
    c=-(C-x_0^2 +2*x_2*x_0-y_0^2+2*y_2*y_0);
    delta=b^2-4*a*c;
    t1=(-b+sqrt(delta))/(2*a)
    t2=(-b-sqrt(delta))/(2*a)
    P=P0(1,:)+v*t2;
end