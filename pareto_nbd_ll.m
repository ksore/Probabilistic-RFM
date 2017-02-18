function [f,g]=pareto_nbd_ll(param)
% pareto_nbd_ll -- Pareto/NBD model log-likelihood
%
% Computes the log likelihood function for the Pareto/NBD model
%
% Syntax: pareto_nbd_ll(param) where the elements of param are r, alpha, s,
% and beta, respectively.
%
% Peter S. Fader (http://petefader.com)
% Bruce G.S. Hardie (http://brucehardie.com)
% Ka Lok Lee (http://kaloklee.com)
%
% Last modified 2005-03-16

global p1x tx T

r     = param(1);
alpha = param(2);
s     = param(3);
beta  = param(4);

maxab = max(alpha,beta);
absab = abs(alpha-beta);
param2 = s+1;
if alpha < beta
    param2 = r+p1x;
end    

part1 = (alpha^r*beta^s/gamma(r))*gamma(r+p1x);
part2 = 1./((alpha+T).^(r+p1x).*(beta+T).^s);
if absab == 0 
   F1=1./((maxab+tx).^(r+s+p1x));
   F2=1./((maxab+T).^(r+s+p1x));
else
   F1=h2f1(r+s+p1x,param2,r+s+p1x+1,absab./(maxab+tx))./...
       ((maxab+tx).^(r+s+p1x));
   F2=h2f1(r+s+p1x,param2,r+s+p1x+1,absab./(maxab+T))./...
       ((maxab+T).^(r+s+p1x));
end

f = -sum(log(part1.*(part2+(s./(r+s+p1x)).*(F1-F2))));
[f/1000 param];
g=[];

end
