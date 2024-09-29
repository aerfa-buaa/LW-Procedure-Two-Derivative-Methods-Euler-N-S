function coneq = Order_MSTDRK(A,Ahat,v,vhat,d,b,L,chars)
%%
if chars==230
    n=2;
    nd=2;
    order=3;
end
if chars==240
    n=2;
    nd=2;
    order=4;
end
if chars==340
    n=3;
    nd=3;
    order=4;
end
% n=step+stage-1;
es=ones(n,1);
es1=ones(nd,1);
% L=[1-step:1:0]';%  L=[ -1; 0] ;
c=A*es+ d*L;   C1=diag(c);
jj=1;
for i = nd
    coneq(jj)  =d(i,:)*es1-1;
    jj=jj+1;
end
stage=nd+1;
coneq(stage)  =b*es1-1;

tau_2 = (c.^2-d*(L).^2)/2 - A*c- Ahat*es;
tau_3 = (c.^3-d*(L).^3)/6 - A*c.^2/2- Ahat*c;
tau_4 = (c.^4-d*(L).^4)/24 -  A*c.^3/6- Ahat*c.^2/2;
tau_5 = (c.^5-d*(L).^5)/120 -  A*c.^4/24- Ahat*c.^3/6;

%     coneq(1)=sum( d*es1)-3;
%     coneq(2)=1-b*es1;
%
%     coneq(3)=1-b*L-v*es;
%     coneq(4)=1/2-b* L.^2/2 -v*c-vhat*es;
%     coneq(5)=1/6-b* L.^3/6 -v*c.^2/2-vhat*c;
%     tau_2 = (c.^2-d*(L).^2)/2 - A*c- Ahat*es;
%     coneq(6)=v*tau_2;



if order==3
    coneq(stage+1)=1-b*L-v*es;
    coneq(stage+2)=1/2-b* L.^2/2 -v*c-vhat*es;
    coneq(stage+3)=1/6-b* L.^3/6 -v*c.^2/2-vhat*c;
    
    coneq(stage+4)=v*tau_2;
end



if order==4
    
    coneq(stage+1)=1-b*L-v*es;
    coneq(stage+2)=1/2-b* L.^2/2 -v*c-vhat*es;
    coneq(stage+3)=1/6-b* L.^3/6 -v*c.^2/2-vhat*c;
    coneq(stage+4)=v*tau_2;
    
    coneq(stage+5)=(1 -b* L.^4)/24 -v*c.^3/6-vhat*c.^2/2;
    coneq(stage+6)=v*C1*tau_2;
    coneq(stage+7)=v*A*tau_2;
    coneq(stage+8)=v* tau_3;
    coneq(stage+9)=vhat* tau_2;
    
end
 
if order==5
    
    coneq(stage+1)=1-b*L-v*es;
    coneq(stage+2)=1/2-b* L.^2/2 -v*c-vhat*es;
    coneq(stage+3)=1/6-b* L.^3/6 -v*c.^2/2-vhat*c;
%     coneq(stage+4)=v*tau_2;
    
    coneq(stage+4)=(1 -b* L.^4)/24 -v*c.^3/6-vhat*c.^2/2;
%     coneq(stage+6)=v*C1*tau_2;
%     coneq(stage+7)=v*A*tau_2;
    coneq(stage+5)=v* tau_3;
%     coneq(stage+9)=vhat* tau_2;
    
% coneq(stage+6)=sum( abs(tau_2));
jj=1;
for kk=step+1 : step+stage-1
coneq(stage+5+jj)=tau_2(kk);
jj=jj+1;
end
 coneq(stage+jj+5)=(1 -b* L.^5)/120 -v*c.^4/24-vhat*c.^3/6;
coneq(stage+jj+6)=v*C1*tau_3;
 coneq(stage+jj+7)=v*A*tau_3;
 coneq(stage+jj+8)=vhat* tau_3;
  coneq(stage+jj+9)=v* tau_4;

end
 
if order==6
    
    coneq(stage+1)=1-b*L-v*es;
    coneq(stage+2)=1/2-b* L.^2/2 -v*c-vhat*es;
    coneq(stage+3)=1/6-b* L.^3/6 -v*c.^2/2-vhat*c;
%     coneq(stage+4)=v*tau_2;
    
    coneq(stage+4)=(1 -b* L.^4)/24 -v*c.^3/6-vhat*c.^2/2;
%     coneq(stage+6)=v*C1*tau_2;
%     coneq(stage+7)=v*A*tau_2;
    coneq(stage+5)=v* tau_3;
%     coneq(stage+9)=vhat* tau_2;
    
% coneq(stage+6)=sum( abs(tau_2));
jj=1;
for kk=step+1 : step+stage-1
coneq(stage+5+jj)=tau_2(kk);
jj=jj+1;
end
 coneq(stage+jj+5)=(1 -b* L.^5)/120 -v*c.^4/24-vhat*c.^3/6;
coneq(stage+jj+6)=v*C1*tau_3;
 coneq(stage+jj+7)=v*A*tau_3;
 coneq(stage+jj+8)=vhat* tau_3;
  coneq(stage+jj+9)=v* tau_4;

 coneq(stage+jj+10)=(1 -b* L.^6)/720 -v*c.^5/120-vhat*c.^4/24;
 coneq(stage+jj+11)=v*A*A*tau_3;
 coneq(stage+jj+12)=v*C1*A*tau_3;
 coneq(stage+jj+13)=v*A*C1*tau_3;
 coneq(stage+jj+14)=v*C1*C1*tau_3;
 coneq(stage+jj+15)=v*Ahat*tau_3;
 coneq(stage+jj+16)=vhat*A*tau_3;
coneq(stage+jj+17)=vhat*C1*tau_3;
coneq(stage+jj+18)=v*C1*tau_4;
 coneq(stage+jj+19)=v*A*tau_4;
 coneq(stage+jj+20)=vhat* tau_4;
  coneq(stage+jj+21)=v* tau_5;

end





end

