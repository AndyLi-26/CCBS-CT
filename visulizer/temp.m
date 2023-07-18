clc;
a([206.714000000,35.774400000]);


function a(C)
A=[203.687000000,36.623700000];B=[182.616000000,32.115700000];
a=norm(C-B);b=norm(A-C);c=norm(A-B);
fprintf("a: %.9f, b: %.9f, c %.9f\n",a,b,c);
cos_theta=(b^2+c^2-a^2)/(2*b*c);
fprintf("cos: %.9f\n",cos_theta);
t=  sqrt((8*4.5*4.5)/(1-cos_theta));
fprintf("t: %.9f",t);
end

