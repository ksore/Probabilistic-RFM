function [f,g] = gamma_gamma_ll( param )
%GAMMA_GAMMA_LL Summary of this function goes here
%   Detailed explanation goes here

global p1x zbar

p = param(1);
q = param(2);
gam = param(3);

func = gammaln(p*p1x+q) - gammaln(p*p1x) - gammaln(q) + ...
    q.*log(gam)+(p*p1x-1).*log(zbar) + (p*p1x) .* log(p1x) - ...
    (p*p1x+q) .* log(gam+p1x.*zbar);
    
accum = 0;

for i = 1:length(func)
    if p1x(i)>0
        accum = accum + func(i);
    end
end

f = - accum;

[f/1000 param];
% g=[];

end

