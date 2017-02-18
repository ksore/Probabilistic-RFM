function [ r alpha s beta ] = Pareto_NBD_params( p1x,  tx  ,T , p2x )
%PARETO_NBD_PARAMS Summary of this function goes here
%   Detailed explanation goes here

lb = .001 * ones(1,4);
ub = 100 * ones(1,4);

initial = ones(1,4);

[params ll] = fmincon('pareto_nbd_ll',initial,[],[],[],[],lb,ub)

r = params(1); alpha = params(2);
s = params(3); beta  = params(4);

disp(strcat('Poisson transaction rate gamma-heterogeneity shape: ',num2str(r)));
disp(strcat('Poisson transaction rate gamma-heterogeneity scale: ',num2str(alpha)));
disp(strcat('Average number of weeks between visits: ',num2str(alpha/r)));

disp(strcat('Exponential dropout rate Gamma-heterogeneity shape: ',num2str(s)));
disp(strcat('Exponential dropout rate Gamma-heterogeneity scale: ',num2str(beta)));
disp(strcat('Average number of weeks until dropout: ',num2str(beta/s)));


end

