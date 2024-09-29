%%%Unpacking Explicit 2 Derivative Runge Kutta Method
%&& We are taking our coeficient vector X which is being fed through
%%% the optimizer and unpacking the pieces back into the original form

%A      SxS matrix and lower triangular because its explicit;
%Ahat   SxS matrix and lower triangular because its explicit;
%b      Vector of Length S
%bhat   Vector of Length S
% n=stage*step+(2*step+stage-2)*(stage-1)+2*(step+stage-1)+1;
%Count defines the number of unknows in each matrix we are building
function [A,Ahat,v,vhat ] =  unpackMSMDRK_all(x, stage,order)


if(order==3)
a21 = x(1);   v2 = x(2); 
 aa21 = a21^2/2; 
 v1 = 1 - v2;
vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2); vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
 
    A=[0   0;
        a21  0;];
    Ahat =[0  0;
        aa21  0;];
    v=[ v1, v2] ;
    vhat=[ vv1, vv2] ;   
end



if(order==4)
a21 = x(1);   a31 = x(2);  a32 = x(3);
aa31 = x(4);    v3 =x(5) ; 
 
 aa21 = a21^2/2; 
 aa32 = (1/2).*(a31.^2+(-2).*a21.*a32+2.*a31.*a32+a32.^2+(-2).*aa31);
 uu = a31 + a32;
 v2=(-1).*a21.^(-3).*(3.*a21.^2.*a32+6.*a21.*aa31+(-3).*a21.*uu.^2+  uu.^3).*v3;
 v1 = 1 - v2 - v3;
 vv1=(-1/12).*a21.^(-1).*uu.^(-1).*((-1)+2.*uu+(-2).*a21.*((-1)+3.*uu+ ...
  a21.*(a21+(-3).*uu).*v2)+2.*(3.*a21+(-1).*uu).*uu.^2.*v3);

vv2=(-1/12).*a21.^(-1).*(a21+(-1).*uu).^(-1).*((-1)+2.*uu+4.*a21.^3.* v2+(-6).*a21.^2.*uu.*v2+(-2).*uu.^3.*v3);
vv3=(-1/12).*uu.^(-1).*((-1).*a21+uu).^(-1).*((-1)+2.*a21+(-2).* ...
  a21.^3.*v2+(-6).*a21.*uu.^2.*v3+4.*uu.^3.*v3);

 v1 = 1 - v2 - v3;
    A=[0  0 0;
        a21  0 0;
        a31  a32 0;];
    Ahat =[0  0 0;
        aa21  0 0;
        aa31  aa32 0;];
    v=[ v1, v2, v3] ;
    vhat=[ vv1, vv2, vv3] ;   
end

end

