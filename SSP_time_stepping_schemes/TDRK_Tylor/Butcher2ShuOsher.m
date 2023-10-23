function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat, r,K)
%function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,K);
%Converting Butcher form to Modified Shu Osher
% w=S*x+dt*T*F(Un)+ dt^2*That*Fdot(Un)
%Y= Un+dt*A*F(Un)+dt^2*Ahat*Fdot(Un)
%Y= R*e*Un+ P*(Un+(dt/r)*F(Un))+Q*(Un+(dt^2/r2)*Fdot(Un))     
%  second derivatives 
%     r2=r^2/K^2; s=length(A); z=zeros(s+1,1); I=eye(s+1);e=ones(s+1,1);
% 	S=[[A;b],z];Shat=[[Ahat;bhat],z];          
%     %Converting butcher to Modified Shu Osher
%     R=inv((I+r*S+r2*Shat));  
%     Re=R*e;
%     P=R*r*S;
%     Q=R*r2*Shat;
% Taylor
    rhat=r /K;
    n=length(A);   z=zeros(n+1,1);   I=eye(n+1);    e=ones(n+1,1);
	   S=[[A;v],z];   Shat=[[Ahat;vhat],z];       
    R=inv( (I+r*S+ 2*rhat*(rhat -r) *Shat) );  
    Re=R*e;
    P=r*R*(S - 2*rhat*Shat   );
    Q=2*rhat*rhat*R *Shat;
end
