%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimization Driver File for Finding optimal SSP Two Derivative Runge
%Kutta Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   rr =opt_mdrk( stage,order,K)
restart=0;   % Start with randomly generated coeficients

% clc
% clear
%
%   restart=2;   % Start with randomly generated coeficients
% X=all_coff(18,:);

% stage=5;          %Number of Stages
% step=3;            %Number of step
% order=5;           %Number of order
%  K=1/sqrt(2);   %Second Derivative Coefficient (dtVV/dtFE)

minreff =1.e-4; %Keep looking until method with at least this value is found
if restart==0
    %Because its generated from random starting point
    %we are being less restrictive on satisfying constraints.
    %     clear A Ahat b bhat x X r   %Be sure all variables are reset
    X1=[0];
    opts=optimset('MaxFunEvals',10000000,'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
        'GradObj','on','MaxIter',10000000,'Diagnostics','on','Display','on',...
        'UseParallel','never','Algorithm','sqp'); % 'UseParallel','never','Algorithm','active-set');  %
    
else
    %     'interior-point'（默认�?�）
    %     'trust-region-reflective'
    %     'sqp'
    %     'sqp-legacy'（仅限于 optimoptions�?
    %     'active-set'
    %Restarting from previously found method
    %In our finetuning of methods use more restrictive tolerences to be sure
    %our conditions are sharply met
    X1 =  intX(step,stage,order);
    %     X1=X;% stores original coefficients from loaded method
    opts=optimset('MaxIter',1000000,'MaxFunEvals',1000000,...
        'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
        'GradObj','on','Diagnostics','off','Display','off',...
        'Algorithm','sqp','UseParallel','never');
end
%%
% n=stage*step+(2*step+stage-2)*(stage-1)+2*(step+stage-1)+1;
if(order==3)
    n=3 ;
end
if(order==4)
    n=6 ;
end
lb=0+zeros(1,n);    lb(end)=-2.8;
ub=1+zeros(1,n);     ub(end)=-0.101;        %requires r>=0


%==============================================
count=0;                                     %Count tracks the number of times optimizer has failed to find a method
info=-2;

%This While loop requires optimizer to keep running until r>minreff is found while satisfying all constraints
while (info==-2 || (r)<minreff || info==0)
    if count==1 %If fails to find a method after 100 times, stop routine
        ('exceed count')
        x=X1;
        r=-x(end);
        r=101;
        break
    end
    %defining initial starting point for Fmincon
    if restart==1
        x0=X1;
    elseif restart==2
        x0=X1+.4*rand(1,length(X1));
    else
        x0=[(2*rand(1,n-1)),-.01];
    end
    
    %==============================================
    %The optimization call:
    %  [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],Aeq,beq,lb,ub,@(x) nlc_mdrk(x,s,p,CC),opts);
    [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],[],[],lb,ub,@(x) nlc_mdrk(x, stage,order,K),opts);
    r=-FVAL;
    count=count+1;
end %while loop
%==============================================
% [con,coneq]=nlc_mdrk(X,s,K,KK);
% MaxViolation=max([con;coneq']);

[A,Ahat,v,vhat] =  unpackMSMDRK_all(X, stage,order);
%  [A,Ahat,v,vhat,d,b] =  unpackMSMDRK(X,s);
coneq = Order_MSTDRK(A,Ahat,v,vhat,stage,order);
r0=-X(end);
[Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat, r0,K);

%  rall=-X(end-2:end);
% [Re,P,Q,T] =ThD_Taylor(A,Ahat,Az,b,bhat,bz,rall,CC,KK);
% [Re,P,Q,T] =ThD_Taylor(A,Ahat,Az,b,bhat,bz,r,K,KK);

rr=X;
% [A ,Ahat, b, bhat] =  unpackMSMDRK(X,s);
% [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,CC);

