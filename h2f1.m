function y = h2f1(a,b,c,z)
% h2f1 -- Gaussian hypergeometric function
%
% Computes the Gaussian hypergeometric function by series expansion,
% iterating to machine epsilon 
%
% Syntax: h2f1(a,b,c,z) where a,b,c,z are scalars or column vectors.
%
% WARNING: this is *very* crude code
%   -- it doesn't perform basic checks such as |z| < 1 or c-a-b > 0 
%      for |z| = 1
%   -- it doesn't recognize special cases such as a = c and b = c 
%   -- it doesn't apply the relevant transformations when |z| is close
%      to 1 (so as to facilitate reliable convergence)
%   etc.
%
% Peter S. Fader (http://petefader.com)
% Bruce G.S. Hardie (http://brucehardie.com)
% Ka Lok Lee (http://kaloklee.com)
%
% Last modified 2005-03-16

lenz = length(z);
j = 0;
uj = ones(lenz,1);
y = uj;
lteps = 0;

while (lteps<lenz/10)
   lasty = y;
   j = j+1;
   uj = uj .*(a+j-1) .*(b+j-1) ./(c+j-1) .*z ./j;
   y = y + uj;
   lteps = sum(y==lasty); 
end


% b= b * ones(lenz,1);
% y = hypergeom([a,b],c,z);