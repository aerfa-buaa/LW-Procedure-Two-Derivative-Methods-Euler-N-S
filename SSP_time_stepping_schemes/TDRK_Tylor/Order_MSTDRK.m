function coneq = Order_MSTDRK(A,Ahat,v,vhat, stage,order)
%%

 
es=ones(stage,1);
 
c=A*es ;   C1=diag(c);

tau_2 = (c.^2 )/2 - A*c- Ahat*es;
tau_3 = (c.^3 )/6 - A*c.^2/2- Ahat*c;
tau_4 = (c.^4 )/24 -  A*c.^3/6- Ahat*c.^2/2;
tau_5 = (c.^5 )/120 -  A*c.^4/24- Ahat*c.^3/6;

 

int=0;

if order==3
    coneq(int+1)=1 -v*es;
    coneq(int+2)=1/2  -v*c -vhat*es;
    coneq(int+3)=1/6  -v*c.^2/2 -vhat*c;
    
    coneq(int+4)=v*tau_2;
end



if order==4
    
    coneq(int+1)=1 -v*es;
    coneq(int+2)=1/2  -v*c -vhat*es;
    coneq(int+3)=1/6  -v*c.^2/2 -vhat*c;
    coneq(int+4)=v*tau_2;
    
    coneq(int+5)=1/24 -v*c.^3/6-vhat*c.^2/2;
    coneq(int+6)=v*C1*tau_2;
    coneq(int+7)=v*A*tau_2;
    coneq(int+8)=v* tau_3;
    coneq(int+9)=vhat* tau_2;
    
end
    
if order>4
disp('order <5');
end
end

