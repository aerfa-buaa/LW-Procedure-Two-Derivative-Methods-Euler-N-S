%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimization Driver File for Finding optimal SSP Two Derivative Runge
%Kutta Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function   rr =opt_mdrk(chars, w ,K)
clc
clear 
 restart=0;   % Start with randomly generated coeficients
chars=340; 
w = [1 1] ;
K=1;

%   restart=1;   % Start with randomly generated coeficients
% X=all_coff(18,:);

%  K=1/sqrt(2);   %Second Derivative Coefficient (dtVV/dtFE)

minreff =1.e-8; %Keep looking until method with at least this value is found
if restart==0
    %Because its generated from random starting point
    %we are being less restrictive on satisfying constraints.
%     clear A Ahat b bhat x X r   %Be sure all variables are reset
    X1=[0];
    opts=optimset('MaxFunEvals',10000000,'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
        'GradObj','on','MaxIter',10000000,'Diagnostics','on','Display','on',...
      'UseParallel','never','Algorithm','sqp'); % 'UseParallel','never','Algorithm','active-set');  %    
  
else
%     'interior-point'锛堥粯璁わ拷?锟斤級
%     'trust-region-reflective'
%     'sqp'
%     'sqp-legacy'锛堜粎闄愪簬 optimoptions锟??
%     'active-set'
    %Restarting from previously found method
    %In our finetuning of methods use more restrictive tolerences to be sure
    %our conditions are sharply met
    X1=X;% stores original coefficients from loaded method
    opts=optimset('MaxIter',100000,'MaxFunEvals',100000,...
        'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
        'GradObj','on','Diagnostics','off','Display','off',...
        'Algorithm','sqp','UseParallel','never');
end
%% 
% n=stage*step+(2*step+stage-2)*(stage-1)+2*(step+stage-1)+1;
%     lb=0+zeros(1,n);    lb(end)=-0.97929615965; 
%     ub=1+zeros(1,n);     ub(end)=-0.199969693973;        %requires r>=0
if chars==23|| chars==241
     n=2;
    lb=-0.2+zeros(1,n);    lb(end)=-0.95; 
    ub=1.5+zeros(1,n);     ub(end)=-0.01 ;        %requires r>=0
end
if chars==231 || chars==239
     n=3;
    lb=-0.2+zeros(1,n);    lb(end)=-0.8; 
    ub=1.2+zeros(1,n);     ub(end)=-0.1 ;        %requires r>=0
end

if chars==230
     n=6;
    lb=-0.2+zeros(1,n);    lb(end)=-0.95; 
    ub=1.5+zeros(1,n);     ub(end)=-0.1 ;        %requires r>=0
end
if chars==240
     n=6;
    lb=-0.2+zeros(1,n);    lb(end)=-0.95; 
    ub=3.+zeros(1,n);     ub(end)=-0.01 ;        %requires r>=0
end
if chars==34
     n=4;
    lb=-0.0001+zeros(1,n);    lb(end)=-0.85; 
    ub=2.0+zeros(1,n);     ub(end)=-0.01 ;        %requires r>=0
end
if chars==341||chars==349
     n=2;
    lb=-0.2+zeros(1,n);    lb(end)=-0.7; 
    ub=1.5+zeros(1,n);     ub(end)=-0.1 ;        %requires r>=0
end

if chars==340
     n=9;
    lb=-0.2+zeros(1,n);    lb(end)=-0.9; 
    ub=2.5+zeros(1,n);     ub(end)=-0.1 ;        %requires r>=0
end


%==============================================
count=0;                                     %Count tracks the number of times optimizer has failed to find a method
info=-2;

%This While loop requires optimizer to keep running until r>minreff is found while satisfying all constraints
while (info==-2 || (r)<minreff || info==0)
    if count==50 %If fails to find a method after 100 times, stop routine
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
        x0=X1+.4*rand(1,length(X));
    else
%         x0=[(2*rand(1,n-1)),-.01];
           x0=[(1*rand(1,n-1)),-.01];
    end

    %==============================================
    %The optimization call:
    %  [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],Aeq,beq,lb,ub,@(x) nlc_mdrk(x,s,p,CC),opts);
    [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],[],[],lb,ub,@(x) nlc_mdrk(x,chars, w,K),opts);
    r=-FVAL;
    count=count+1;
end %while loop
%==============================================
r=-X(end);
[A,Ahat,v,vhat,d,b,L] =  unpackTDMSRK_all(X,chars, w);
[Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat,d,b,r,K);
xx=[Re P Q ];
con=-xx(:); %require entries in xx and r to be positive
rr=X;




 