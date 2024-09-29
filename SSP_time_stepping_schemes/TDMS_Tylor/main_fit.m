clear all
% close all
warning off
clc

K=1/sqrt(2);
% K=1
chars=2239;  %TDTS34, TDTS23,TDTSRK223
%%
 if chars==223||chars==2239
    wq=0.75;
    for ii =1:17
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.05
    end
    y_b1 = all_coff(:,2);   y_v1 = all_coff(:,3); y_v3 = all_coff(:,4);   y_a22 = all_coff(:,5);
 end
 
 if chars==224  || chars==2249
    wq=1;
    for ii =1: 9
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.0
    end
    y_b1 = all_coff(:,2);   y_v3 = all_coff(:,3); y_vv3 = all_coff(:,4); y_a21 = all_coff(:,5);  y_d31 = all_coff(:,6);   
 end
 

if chars==341 || chars==349||chars==34||chars==340
    a_jj=9;
    wq1=0.65;   %wq1=0.65;
    for ii =1:17
        wq2=0.75;   %wq1=0.65;
        for jj =1:a_jj
            ii
            jj
            w=[wq1, wq2];
            all_coff((ii-1)*a_jj+jj,:)=[w,opt_mdrk(chars, w ,K)];
            wq2=wq2+0.1
        end
        wq1=wq1+0.1
    end
    x_b1 = all_coff(:,3);   x_b2 = all_coff(:,4);  
    x_wq2=all_coff(:,2);
end
 

 if chars==231||chars==23||chars==230 || chars==239
    wq=0.7;
    for ii =1:19
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.05
    end
    y_b1=all_coff(:,2);
end
 
x_wq1=all_coff(:,1);


% k1=[0.1, 0.2, 0.3,  0.4,  0.5, 0.6, 1/sqrt(2), 0.8  1,  1.5,  2,  3,  4]; ii=0;
% [hh,nb0]=size(k1);
% nb1=40;
%
% method=25;
% for iqq0=1:nb1
% %           K=k1(iqq0);
%           K=0.1*iqq0
%        opt_r =opt_mdrk(method,K);
%         ii=ii+1;
%         all_coff(ii,1)=K;
%         [hh21,r32]=size(opt_r);
%         for  ci=1:r32
%         all_coff(ii,1+ci)= opt_r(ci) ;
%         end
%     end
