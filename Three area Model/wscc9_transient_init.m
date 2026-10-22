%  This program obtains th Bus Admittance Matrix for power flow solution
%  Copyright (c) 1998-2010 by H. Saadat
%linedata for extended ybus

%               |  From |  To   |   R     |   X     |   B/2     |  1 for line code 
%               |  Bus  | Bus   |         |         |           | or tap setting value |
     linedataa = [ 1       4        0        0.1184       0           1;
                   4       6        0.017    0.092      0.079         1;
                   4       5        0.010    0.085       0.088        1;
                   6       9        0.039    0.170      0.179         1;
                   5       7        0.032    0.161      0.153         1;                                          ; 
                   7       2        0        0.1823     0             1;
                   7       8        0.0085   0.072      0.0745        1;
                   8       9        0.0119   0.1008     0.1045        1;
                   9       3        0       0.2399      0             1] ;

j=sqrt(-1); i = sqrt(-1);
nl = linedataa(:,1); nr = linedataa(:,2); R = linedataa(:,3);
X = linedataa(:,4); Bc = j*linedataa(:,5); a = linedataa(:, 6);
nbr=length(linedataa(:,1)); nbus = max(max(nl), max(nr));
Z = R + j*X; y= ones(nbr,1)./Z;        %branch admittance
for n = 1:nbr
if a(n) <= 0  
    a(n) = 1; else end
Ybusext=zeros(nbus,nbus);     % initialize Ybus to zero
               % formation of the off diagonal elements
for k=1:nbr;
       Ybusext(nl(k),nr(k))=Ybusext(nl(k),nr(k))-y(k)/a(k);
       Ybusext(nr(k),nl(k))=Ybusext(nl(k),nr(k));
    end
end
              % formation of the diagonal elements
for  n=1:nbus
     for k=1:nbr
         if nl(k)==n
         Ybusext(n,n) = Ybusext(n,n)+y(k)/(a(k)^2) + Bc(k);
         elseif nr(k)==n
         Ybusext(n,n) = Ybusext(n,n)+y(k) +Bc(k);
         else, end
     end
end
clear Pgg
Ybusext;
bmva=100;
Ybusext(5,5)=Ybusext(5,5)+((125.0-50.0i)/bmva)/(0.996^2);
Ybusext(6,6)=Ybusext(6,6)+((90.0-30.0i)/bmva)/(1.013^2);
Ybusext(8,8)=Ybusext(8,8)+((100.0-35.0i)/bmva)/(1.016^2);
Ybusext
n=3;
[r,c]=size(Ybusext);
K=Ybusext(1:n,1:n);
L=Ybusext(1:n,n+1:c);
M=Ybusext(n+1:r,n+1:c);
LT=transpose(L);
Y1=K-(L*inv(M)*LT)
%Faults Occurs at bus 7
y2=Ybusext;
y2(7,:)=[];
y2(:,7)=[];
[r,c]=size(y2);
K=y2(1:n,1:n);
L=y2(1:n,n+1:c);
M=y2(n+1:r,n+1:c);
LT=transpose(L);
Y2=K-(L*inv(M)*LT)
%Post fault
y3=Ybusext;
% Remove line 5-7: cancel its series admittance and its shunt charging (B/2 at each end)
assert(isequal(linedataa(5,1:2),[5 7]), 'Row 5 of linedataa must be line 5-7');
Bc57=1i*linedataa(5,5);
y3(5,5)=y3(5,5)+y3(5,7)-Bc57;
y3(7,7)=y3(7,7)+y3(7,5)-Bc57;
y3(5,7)=0;
y3(7,5)=0;
[r,c]=size(y3);
K=y3(1:n,1:n);
L=y3(1:n,n+1:c);
M=y3(n+1:r,n+1:c);
LT=transpose(L);
Y3=K-(L*inv(M)*LT)


%          bus  v    theta   p            q
loadflow =[1  1.040  0     71.6/100     27/100;
           2  1.025  9.3   163/100      6.7/100;
           3  1.026  4.7   85/100       -10.9/100];

V=loadflow(:,2);
theta=loadflow(:,3);
Pg=loadflow(:,4);
Qg=loadflow(:,5);

sb=100;%base mva
Pm1=(Pg(1,1));
Pm2=(Pg(2,1));
Pm3=(Pg(3,1));
    qq1=(Qg(1,1));
    qq2=(Qg(2,1));
    qq3=(Qg(3,1));
    Pm=[Pm1;Pm2;Pm3];
    qg=[qq1;qq2;qq3];
    xd1=0.0608;
    xd2=0.1198;
    xd3=0.1813;
    X=[xd1;xd2;xd3];
    v=V;
    alpha=theta;
    e=zeros(n,1);
    n=3
    for i=1:n
        e(i)=v(i)+(qg(i)*X(i))/(v(i))+(1i*(Pm(i))*X(i)/(v(i)));
    E=abs(e)
    del(i)=(angle(e(i))*180/pi+alpha(i));
    del(i)=del(i)*pi/180
end
E1=E(1);
E2=E(2);
E3=E(3);
D10=del(1);
D20=del(2);
D30=del(3);
H=[23.645 ; 6.4 ; 3.01];
f=60;


