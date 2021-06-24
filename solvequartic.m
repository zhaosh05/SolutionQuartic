function [xr,yr]=solvequartic(A,B)
%input: coefficients of the two equations as follows
%A=[a1,b1,c1,d1,e1,f1]; B=[a2,b2,c2,d2,e2,f2];
%the two equations should be like this:
% a1*x^2 + b1*x*y + c1*y^2 + d1*x + e1*y + f1=0
% a2*x^2 + b2*x*y + c2*y^2 + d2*x + e2*y + f2=0

zerothr=1e-20;%threshold to judge if a value is zero. adjustable according to specific problems
xr=NaN(4,1);%four root storage for x
yr=NaN(4,1);%four root storage for y
%a*x^2 + b*x*y + c*y^2 + d*x + e*y + f=0
a1=A(1);b1=A(2);c1=A(3);d1=A(4);e1=A(5);f1=A(6);
a2=B(1);b2=B(2);c2=B(3);d2=B(4);e2=B(5);f2=B(6);
% eliminate y^2 term to form (t1*x+t2)y=t3*x^2+t4*x+t5
t1=b1*c2-b2*c1;t2=e1*c2-e2*c1;t3=-(a1*c2-a2*c1);t4=-(d1*c2-d2*c1);t5=-(f1*c2-f2*c1);
%handle t1*x+t2==0
%two cases: t1~=0 and t1==0. 分两种情况，t1不是0 和 t1是0
if abs(t1)>zerothr %if t1 is non-zero
    %determine if -t2/t1 is a solution. 先判断-t2/t1是不是解
    xtmp=t2/t1;
    if abs(t3*xtmp^2+t4*xtmp+t5)<zerothr
        % ay*y^2+by*y+cy=0
        ay=c1; by=b1*xtmp+e1;cy=a1*xtmp^2+d1*xtmp+f1;
        if abs(ay)<zerothr
            if abs(by)<zerothr
                % no solution
                xr=NaN;
                yr=NaN;
            else%one solution
                xr=xtmp;
                yr=-cy/by;
            end
        else%two solutions
            xr(1)=xtmp; xr(2)=xtmp;
            yr(1)=(-by+sqrt(by^2-4*ay*cy))/2/ay;
            yr(2)=(-by-sqrt(by^2-4*ay*cy))/2/ay;
        end
        return
    else
        % no such solution
    end
else
    %t1==0 and t2==0, then degenerate to t3*x^2+t4*x+t5=0
    %t1是0，则t2也是0，此时方程退化为t3*x^2+t4*x+t5=0
    if abs(t3)<zerothr
        if abs(t4)<zerothr
            %no solution
            return
        else
            xr(1)=-t5/t4;
            ay=c2;by=e2+b2*xr(1);cy=f2+a2*xr(1)^2+d2*xr(1);
            if abs(ay)<zerothr
                if abs(by)<zerothr
                    %no solution
                    xr(1)=NaN;
                else
                    yr(1)=-cy/by;
                end
                return
            else
                xr(2)=xr(1);
                yr(1)=(-by+sqrt(by^2-4*ay*cy))/2/ay;
                yr(2)=(-by-sqrt(by^2-4*ay*cy))/2/ay;
                return
            end
        end
    else
        xr(1)=(-t4+sqrt(t4^2-4*t3*t5))/2/t3;
        xr(2)=(-t4-sqrt(t4^2-4*t3*t5))/2/t3;
        for j=1:2
            ay=c2;by=e2+b2*xr(j);cy=f2+a2*xr(j)^2+d2*xr(j);%use 2nd equation to solve y
            if abs(ay)<zerothr
                if abs(by)<zerothr
                    %no solution
                    xr=NaN;yr=NaN;
                else
                    yr(j)=-cy/by;
                end
                return
            else
                yr(2*(j-1)+1)=(-by+sqrt(by^2-4*ay*cy))/2/ay;
                yr(2*(j-1)+2)=(-by-sqrt(by^2-4*ay*cy))/2/ay;
            end
        end
        xr(3)=xr(2);xr(4)=xr(3);xr(2)=xr(1);
        return
    end
end

% case t1*x+t2 != 0, then form the quartic equation in 1 variable: 
% a*x^4+b*x^3+c*x^2+d*x + e = 0
%general solution for four roots. 一般形式的四个根解法
e=f1*t2^2 + e1*t2*t5 + c1*t5^2;
d=d1*t2^2 + b1*t2*t5 + 2*c1*t4*t5 + e1*t1*t5 + e1*t2*t4 + 2*f1*t1*t2;
c=c1*(t4^2 + 2*t3*t5) + a1*t2^2 + f1*t1^2 + b1*t1*t5 + b1*t2*t4 + 2*d1*t1*t2 + e1*t1*t4 + e1*t2*t3;
b=d1*t1^2 + 2*a1*t1*t2 + b1*t1*t4 + b1*t2*t3 + 2*c1*t3*t4 + e1*t1*t3;
a=a1*t1^2 + b1*t1*t3 + c1*t3^2;

delta0 = c^2 - 3*b*d + 12*a*e;
delta1 = 2*c^3 - 9*b*c*d + 27*b^2*e + 27*a*d^2 - 72*a*c*e;
disc=-(delta1^2-4*delta0^3)/27;

Q=((delta1+sqrt(-disc*27))/2)^(1/3);
p=(8*a*c-3*b^2)/(8*a*a);
q=(b^3 - 4*a*b*c + 8*a*a*d)/(8*a^3);
tmp=(-2/3)*p+1/(3*a)*(Q+delta0/Q);
S=0.5*sqrt(tmp);

tmp1=-4*S*S-2*p+q/S;
tmp2=-4*S*S-2*p-q/S;
stmp1=sqrt(tmp1);
stmp2=sqrt(tmp2);
ba=b/4/a;
xr(1)=-ba-S+0.5*stmp1;
xr(2)=-ba-S-0.5*stmp1;
xr(3)=-ba+S+0.5*stmp2;
xr(4)=-ba+S-0.5*stmp2;

yr=(t3*xr.^2+t4*xr+t5)./(t1*xr+t2);
