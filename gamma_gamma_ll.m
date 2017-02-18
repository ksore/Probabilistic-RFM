function [f,g] = gamma_gamma_ll( param )
%GAMMA_GAMMA_LL Summary of this function goes here
%   Detailed explanation goes here

global p1x zbar

p = param(1);
q = param(2);
gam = param(3);

func = real(gamma(p*p1x+q)./(gamma(p*p1x)*gamma(q)) .* ...
        zbar.^(p*p1x-1) .* gam.^q .* p1x.^(p*p1x) ./ ...
        (gam+p1x.*zbar).^(p*p1x+q));
    
accum = 0;

for i = 1:length(func)
    if func(i) > 0 && zbar(i)>0
        if func(i) == Inf
            a = 1;
        end
        accum = accum + real(log(func(i)));
    end
end

f = - accum;

[f/1000 param];
% g=[];

end

