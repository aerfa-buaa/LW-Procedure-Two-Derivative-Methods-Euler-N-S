clear all
% close all
warning off
clc

stage=2;          %Number of Stages
order=3;           %Number of order
K=1/sqrt(2); 
K=1.0;

% for ii =1:10
%     ii
%    all_coff(ii,:)=opt_mdrk( stage,order,K);
% end


 ii=0;
for iqq0=1:20
 
       K=0.1*iqq0
       opt_r =opt_mdrk( stage,order,K);
        ii=ii+1;
        all_coff(ii,1)=K;   
        [hh21,r32]=size(opt_r);
        for  ci=1:r32
        all_coff(ii,1+ci)= opt_r(ci) ;  
        end
    end
