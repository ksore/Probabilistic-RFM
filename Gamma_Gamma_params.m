function [ p q gamma ] = Gamma_Gamma_params( p1x, zbar )
%GAMMA_GAMMA_PARAMS Summary of this function goes here
%   Detailed explanation goes here

lb = [0.001,0.001,0.001];
ub = [50,50,500];

initial = [1,1,1];

lloptions = optimoptions(@fmincon,'Display','final','Algorithm','interior-point', ...
    'OptimalityTolerance', 0.0001, 'Display','off');

problem = createOptimProblem('fmincon', 'objective','gamma_gamma_ll','x0',initial, ...
    'lb',lb,'ub',ub,'options',lloptions);

ms = MultiStart;
stpoints  = RandomStartPointSet('NumStartPoints',50);
list(stpoints, problem);
[params,ll] = run(ms,problem,50);

p = params(1);
q = params(2);
gamma = params(3);

disp(strcat('Gamma spend shape (common): ',num2str(p)));
disp(strcat('Gamma spend scale gamma-heterogeneity shape: ',num2str(q)));
disp(strcat('Gamma spend scale gamma-heterogeneity scale: ',num2str(gamma)));

disp(strcat('Mean spend: ',num2str(p*gamma/(q-1))));
if q>2
    qstdev = num2str((p^2*gamma^2)/((q-1)^2*(q-2)));
else
    qstdev = 'Undefined';
end
disp(strcat('Standard Deviation of spend: ',qstdev));

disp(strcat('Mode spend: ',num2str(p*gamma/(q+1))));

end

