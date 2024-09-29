function [r,g]=mdrk_am_obj(coeffs)
% Used by opt_mdrk for Objective Function
% Same as Ketchesons original code

 r=coeffs(end);
g=zeros(size(coeffs));
g(end)=1;

%  rall=coeffs(end-2:end);
%  r=max(rall);
% g=zeros(size(coeffs));
% g(end-2:end)=1;