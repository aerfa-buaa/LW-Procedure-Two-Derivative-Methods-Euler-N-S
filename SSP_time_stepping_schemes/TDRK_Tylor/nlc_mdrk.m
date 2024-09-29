function [con,coneq]=nlc_mdrk(x, stage,order,K)
% The Nonlinear Constraints for our optimization routine are a combination
% of the equality constraints, coneq, which come from the order conditions,
% and Inequality Constraints, con, which are from our SSP conditions. 
% IE Modified Shu Osher Decomposition having positive Coeficients

%=====================================================
% Extract arrays A,d,b,theta from x
%Unpacking decision variables from x
% KK=0.6888921;
 r=-x(end);
%r=-x(end-2:end);
 [A,Ahat,v,vhat ] =  unpackMSMDRK_all(x, stage,order);
% [A,Ahat,v,vhat,d,b] =  unpackMSMDRK(x,s);
[Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat, r,K);
% xx=[Re P Q];  
xx=[Re P Q ];  
 con=-xx(:); %require entries in xx and r to be positive
%con=[ ];
% coneq = Order_MSTDRK(A,Ahat,v,vhat, stage,order);
  coneq=[ ];

% [A,Ahat,b,bhat] =  unpackMSMDRK(x,s);
% [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,K);
% xx=[Re P Q];  
% coneq= oc_mdrk(p,x,s);


%This follows our SSP conditions defined in our paper. 
%Fmincon Syntax requires -xx<=0 to guarantee positivity of our coeficients

end % of function

%%% Converting Butcher form to Modified Shu Osher
%Y= Un+dt*A*F(Un)+dt^2*Ahat*Fdot(Un)
%Y= R*e*Un+ P*(Un+(dt/r)*F(Un))+Q*(Un+(dt^2/r2)*Fdot(Un))

%function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,K);
%   
%     r2=r^2/K^2; s=length(A); z=zeros(s+1,1); I=eye(s+1);e=ones(s+1,1);
%     S=[[A;b],z];Shat=[[Ahat;bhat],z];          
%     %Converting butcher to Modified Shu Osher
%     R=inv((I+r*S+r2*Shat));  
%     Re=R*e;
%     P=R*r*S;
%     Q=R*r2*Shat;
%end
