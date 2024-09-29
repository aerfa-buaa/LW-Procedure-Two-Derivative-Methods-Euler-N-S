%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimization Driver File for Finding optimal SSP Two Derivative Runge
%Kutta Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function   rr =opt_mdrk(step,stage,order,K)
%  restart=0;   % Start with randomly generated coeficients

% clc
% clear
% 
  restart=1;   % Start with randomly generated coeficients
 
X=all_coff(4,1:end);
% restart=0;
% stage=4;          %Number of Stages
% step=2;            %Number of step
% order=6;           %Number of order
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
%     'interior-point'锛堥粯璁�??硷級
%     'trust-region-reflective'
%     'sqp'
%     'sqp-legacy'锛堜粎闄愪簬 optimoptions�??
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
n=6 ;
    lb=0+zeros(1,n);    lb(end)=-2.95; 
    ub=1+zeros(1,n);     ub(end)=-0.11;        %requires r>=0


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
        x0=[(2*rand(1,n-1)),-.01];
    end

    %==============================================
    %The optimization call:
    %  [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],Aeq,beq,lb,ub,@(x) nlc_mdrk(x,s,p,CC),opts);
    [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],[],[],lb,ub,@(x) nlc_mdrk(x, stage,order,K),opts);
    r=-FVAL;
    count=count+1;
end %while loop
%%
 

%  xx=[ 0.0578228747503574	0.0565972858288410	0.146905532577135	0.942177125249643	0.943402714171159	0.853094467422865	0.0828598199039427	0.917140180096057	0.0407464568134064	0.415458092169921	0.0193495826801274	0.557557422416298	0.0633604062678195	0.0802727079208418	0.483066613222576	0.154564026679913	0.380955590363909	0.0391297353822429	0.504773804553181	0.197658457652033	0.274086523156778	0.0672112991597077	0.000198358920408367	0.0909906397001514	0.00752748797227276	0.0397638165626754	0.0888534182229604	0.00875585529516799	0.0436172954725402	0.0670832348077086	0.0565433263071460	0.00536789300350426	0.0523173463834105	0.0579243154859367	0.0684151932895744	0.0109734294730797	-1.40000000000000];
%      X=xx;
% X=all_coff(19,1:end-1);

  [A,Ahat,v,vhat ] =  unpackMSMDRK_all(X, stage,order);
%  [A,Ahat,v,vhat,d,b] =  unpackMSMDRK(X,s);
   coneq = Order_MSTDRK(A,Ahat,v,vhat, stage,order)';
    r0=-X(end);
  [Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat, r0,K);

%  rall=-X(end-2:end);
% [Re,P,Q,T] =ThD_Taylor(A,Ahat,Az,b,bhat,bz,rall,CC,KK);
% [Re,P,Q,T] =ThD_Taylor(A,Ahat,Az,b,bhat,bz,r,K,KK);

rr=X;
% [A ,Ahat, b, bhat] =  unpackMSMDRK(X,s);
% [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,CC);

 