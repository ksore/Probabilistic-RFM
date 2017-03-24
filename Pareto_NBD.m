%% Building Study Set and Cumulative Purchases

% Basic assumptions for calibration duration and cohort size

weeks_in_cohort = 10;
weeks_in_study = 40;

cumrptslsfig = figure('Name','Cumulative repeat sales');
rptslsfig = figure('Name','Repeat sales');
spendmodel = figure('Name','Spend Model');
hold

first_date = min(patients.first_visit);
patients.cohort = ceil((patients.first_visit-first_date+1)/(weeks_in_cohort*7));

for i = 1:size(patients)
       visits_of_patient = sortrows(visits(visits.patient_number==patients.patient_number(i),'date_completed_serial'));
       
       study_begin = first_date + ((patients.cohort(i)-1) * weeks_in_cohort)*7;
       study_end = first_date + ((patients.cohort(i)-1) * weeks_in_cohort + weeks_in_study) *7;
       holdout_end= first_date + ((patients.cohort(i)-1) * weeks_in_cohort + weeks_in_study*2) *7;
       
       studyCondition = visits_of_patient.date_completed_serial >= study_begin & ...
           visits_of_patient.date_completed_serial < study_end;
       visits_of_patient_in_study = visits_of_patient( studyCondition, 'date_completed_serial');
       
       holdoutCondition = visits_of_patient.date_completed_serial >= study_end & ...
           visits_of_patient.date_completed_serial < holdout_end;
       visits_of_patient_in_holdout = visits_of_patient(holdoutCondition, 'date_completed_serial');

       patients(i,'study_visits') = num2cell(height(visits_of_patient_in_study)-1);
       patients(i,'holdout_visits') = num2cell(height(visits_of_patient_in_holdout));
       patients(i,'first_visit') = num2cell(visits_of_patient_in_study.date_completed_serial(1));
       patients(i,'observation_time') = num2cell((study_end-patients.first_visit(i) )/7);
       if patients.study_visits(i) ~= 0
           patients(i,'last_repeat') = num2cell((visits_of_patient_in_study.date_completed_serial( ...
               height(visits_of_patient_in_study)) - patients.first_visit(i)  )/7);
       end

end

visits = join(visits,patients,'RightVariables',{'cohort','first_visit'});
visits.Properties.VariableNames([7 8]) = {'cohort' 'first_visit'};
visits.week = ceil((visits.date_completed_serial-first_date+1)/7);
visits.is_first = visits.date_completed_serial== visits.first_visit;

clear i visits_of_patient study_begin study_end holdout_end studyCondition visits_of_patient_in_study ...
    holdoutCondition visits_of_patient_in_holdout;

for i = 1:max(patients.cohort)
    cohortCondition = visits.cohort == i & visits.is_first ==0;
    cohort_visits = accumarray(table2array(visits(cohortCondition,'week')), 1);
    weekly_visits(:,i) = cohort_visits(:,1);
    weekly_cumul_visits(:,i) = cumsum(cohort_visits(:,1));
end

clear i cohortCondition cohort_visits

%% Narrow Down Cohort

global p1x p2x tx T zbar
cohort_results = table;

cohort_n = 5;

for i = 1:cohort_n
    
    p1x = table2array(patients(patients.cohort == i,'study_visits'));
    tx = table2array(patients(patients.cohort == i,'last_repeat'));
    T = table2array(patients(patients.cohort == i,'observation_time'));
    p2x = table2array(patients(patients.cohort == i,'holdout_visits'));
    zbar = table2array(patients(patients.cohort == i,'mean_rev'));
    
    % Calculate Pareto-NBD Frequency and Churn Model Coefficients
    
    [r2, alpha2, s2, beta2] = Pareto_NBD_params(p1x,tx,T,p2x);
    
    % Calculate Gamma-Gamma Spend Model Coefficients
    
    [p2, q2, gam2] = Gamma_Gamma_params(p1x, zbar);
    
    T2 = table(i,r2,alpha2,s2,beta2, p2,q2,gam2, ...
        alpha2/r2, beta2/s2, p2*gam2/(q2-1),p2*gam2/(q2+1),...
        'VariableNames', ...
        {'cohort','r','alpha','s','beta', 'p','q','gamma', ...
        'weeks_between_visits', 'life_expectancy', 'average_spend', 'mode_spend'});
    
    cohort_results = [cohort_results;T2];
    
end

clear i T2 r2 s2 alpha2 beta2 p2 q2 gam2;

%% Plots

figure(spendmodel);
x = 1:1000;
clf(spendmodel); 

for i = 1:cohort_n
    p1x = table2array(patients(patients.cohort == i,'study_visits'));
    tx = table2array(patients(patients.cohort == i,'last_repeat'));
    T = table2array(patients(patients.cohort == i,'observation_time'));
    p2x = table2array(patients(patients.cohort == i,'holdout_visits'));
    zbar = table2array(patients(patients.cohort == i,'mean_rev'));
        
    subplot(floor(sqrt(cohort_n)),floor(sqrt(cohort_n))+1,i)
    histogram(zbar,'Normalization','pdf', 'BinLimits',[0,1000],'BinWidth',10);
    hold on;
    plot(x,(cohort_results.p(i).*cohort_results.gamma(i)).^cohort_results.q(i) ...
        .*x.^(-cohort_results.q(i)-1).*exp(-cohort_results.p(i).*cohort_results.gamma(i)./x) ...
        ./gamma(cohort_results.q(i)));
    title(['Cohort ',num2str(i), ...
        ' - Mode: $', num2str(round(cohort_results.mode_spend(i),0)), ...
        ' - Average: $', num2str(round(cohort_results.average_spend(i),0))]);
end

hold off;
clear i x;

%% Script to create the cumulative and incremental tracking plots

days_in_cohort = weeks_in_cohort * 7;
days_in_study = weeks_in_study * 7;

x_axis_max = 80;
y_axis_max_rpt=700;
y_axis_max_cum=20000;

for i = 1:cohort_n
    p1x = table2array(patients(patients.cohort == i,'study_visits'));
    tx = table2array(patients(patients.cohort == i,'last_repeat'));
    T = table2array(patients(patients.cohort == i,'observation_time'));
    p2x = table2array(patients(patients.cohort == i,'holdout_visits'));
    zbar = table2array(patients(patients.cohort == i,'mean_rev'));
    
    % determine cohort size by day of trial
    ns = [];
    for j = 1:days_in_cohort
        ns(j) = sum( round(T,2) == round((days_in_study+1-j)/7,2));
    end

    % generate sales cumulative forecast
    endwk = weeks_in_study * 2 ;
    endday = endwk*7; 

    tmp1 = cohort_results.r(i)*cohort_results.beta(i)/ ...
        (cohort_results.alpha(i)*(cohort_results.s(i)-1));
    tmpcumsls1 = [];
    for j = 1:endday
        tmp2 = (cohort_results.beta(i)/(cohort_results.beta(i)+j/7))^(cohort_results.s(i)-1);
        tmpcumsls1(j) = tmp1*(1-tmp2);
    end

    tmpcumsls2 = zeros(days_in_cohort,endday);
    for j = 1:days_in_cohort
        tmpcumsls2(j,:) = [ zeros(1,j) tmpcumsls1(1:endday-j) ];
    end

    cumrptsls = [];
    dailysls = ns*tmpcumsls2;
    for j = 1:endwk
        cumrptsls(j) = dailysls(j*7);
    end     

    chart_offset = (i-1)*weeks_in_cohort;

    % create tracking plot of cumulative repeat sales (pred. vs actual)

    figure(cumrptslsfig);
    subplot(floor(sqrt(cohort_n)),floor(sqrt(cohort_n))+1,i)

    plot(1:endwk,weekly_cumul_visits(1:endwk,i),'k', ...
        1+chart_offset:endwk+chart_offset,cumrptsls,'b--');
    axis([0 x_axis_max 0 y_axis_max_cum]);
    xlabel('Week'); ylabel('Cumulative Repeat Transactions');
    legend('Actual','Model','Location','SouthEast');
    % print -depsc 'cumrptsls.eps'

    % create tracking plot of weekly repeat sales (pred. vs actual)
    figure(rptslsfig);
    subplot(floor(sqrt(cohort_n)),floor(sqrt(cohort_n))+1,i)

    incrptsls = [ cumrptsls(1) diff(cumrptsls) ];
    plot(1:endwk,weekly_visits(1:endwk,i),'k', ...
        1+chart_offset:endwk+chart_offset,incrptsls,'b--');
    axis([0 x_axis_max 0 y_axis_max_rpt]);
    xlabel('Week'); ylabel('Weekly Repeat Transactions');
    legend('Actual','Model','Location','SouthEast');
    % print -depsc 'incrptsls.eps'

end

clear ns tmp1 tmp2 tmpcumsls1 tmpcumsls2 chart_offset i j x_axis_max ...
    y_axis_max_rpt y_axis_max_cum cumrptsls incrptsls dailysls;

% %% 
% % Script to compute P(active|p1x,tx,T) for Pareto/NBD model.
% % Also creates a plot comparing average P(active|p1x,tx,T) with the
% % observed proportion of customers active in the second period by p1x.
% %
% 
% % compute P(active|p1x,tx,T)
% 
% maxab = max(alpha,beta);
% absab = abs(alpha-beta);
% param2 = s+1;
% if alpha < beta
%     param2 = r+p1x;
% end    
% 
% F0 = (alpha+T).^(r+p1x).*(beta+T).^s;
% F1=h2f1(r+s+p1x,param2,r+s+p1x+1,absab./(maxab+tx))./...
%     ((maxab+tx).^(r+s+p1x));
% F2=h2f1(r+s+p1x,param2,r+s+p1x+1,absab./(maxab+T))./...
%     ((maxab+T).^(r+s+p1x));
% pactive = 1./(1+(s./(r+s+p1x)).*F0 .*(F1-F2));
% 
% % compute average P(active|p1x,tx,T) and determine the proportion of
% % customers buying in the second 39 weeks for each level of p1x
% pa_actual = zeros(max(p1x)+1,1);
% pa_est = zeros(max(p1x)+1,1);
% np1x = zeros(max(p1x)+1,1);
% for y = unique(p1x)'
%     isx = find(p1x==y);
%     np1x(y+1) = length(isx);
%     pa_actual(y+1) = sum(p2x(isx)>0)/np1x(y+1);
%     pa_est(y+1) = sum(pactive(isx))/np1x(y+1);
% end
% 
% % create right-censored version for plot
% censor = 12;  % right-censor at 12+
% denom = sum(np1x(censor+1:length(np1x)));
% 
% pa_act_cen = pa_actual(1:censor);
% pa_act_cen(censor+1) = (np1x(censor+1:length(np1x))'*...
%     pa_actual(censor+1:length(np1x)))/denom;
% 
% pa_est_cen = pa_est(1:censor);
% pa_est_cen(censor+1) = (np1x(censor+1:length(np1x))'*...
%     pa_est(censor+1:length(np1x)))/denom;
% 
% subplot(2,2,3)
% 
% plot(0:censor,pa_act_cen,'k',0:censor,pa_est_cen,'kp--');
% legend('Empirical','Pareto/NBD','Location','SouthEast');
% xlabel(strcat('Transactions in Study Period (',string(weeks_in_study),' weeks)')); 
% ylabel('Probability of still being a customer');
% axis([-.3 12.3 0 1]);
% xticks([0 1 2 3 4 5 6 7 8 9 10 11 12])
% xticklabels({ ' 0', ' 1', ' 2', ' 3', ' 4', ' 5', ' 6', '7', '8', '9', '10', '11', '12+' });
% grid on
% grid minor
% % print -depsc 'pactive_grouped.eps'
% 
% clear denom y isx np1x F0 F1 F2 absab

%% 
% Script to compute conditional expectations
%

t = weeks_in_study;  % period for which conditional expectations are to be computed

for i = 1:size(patients)
    cohort = patients.cohort(i);
    study_visits = patients.study_visits(i);
    obs_time = patients.observation_time(i);
    lst_rpt = patients.last_repeat(i);
    avg_purchase = patients.mean_rev(i);
    
    if cohort <=cohort_n
        
        r = cohort_results.r(cohort);
        alpha= cohort_results.alpha(cohort);
        s= cohort_results.s(cohort);
        beta= cohort_results.beta(cohort);
        p= cohort_results.p(cohort);
        q= cohort_results.q(cohort);
        gam= cohort_results.gamma(cohort);

        % Pareto/NBD Expectations
        
        maxab = max(alpha,beta);
        absab = abs(alpha-beta);
        param2 = s+1;
        if alpha < beta
            param2 = r+study_visits;
        end  

        F0 = (alpha+obs_time)^(r+study_visits)*(beta+obs_time)^s;
        F1 = h2f1(r+s+study_visits,param2,r+s+study_visits+1,absab/(maxab+lst_rpt))/...
            ((maxab+lst_rpt)^(r+s+study_visits));
        F2 = h2f1(r+s+study_visits,param2,r+s+study_visits+1,absab/(maxab+obs_time))/...
            ((maxab+obs_time)^(r+s+study_visits));
        A0 = F1-F2;
        pactive = (1+(s/(r+s+study_visits))*F0*A0)^(-1);

        tmp1 = (r+study_visits)*(beta+obs_time)/(alpha+obs_time)/(s-1);
        tmp2 = ((beta+obs_time)/(beta+obs_time + t))^(s-1);
        ce = tmp1*(1-tmp2)*pactive;

        patients.exp_tr(i) = ce;
        
        % Gamma-gamma Expectations

        patients.exp_sp(i) = p*(gam+(study_visits+1)*avg_purchase) / ...
            (p*(study_visits+1)+q-1);
        
        % Discounted Expected Transactions
        
        delta = log(1+0.1)/52;
        
        L = gamma(r+study_visits)*alpha^r*beta^s/gamma(r) * ...
            (1/F0 + (s/(r+s+study_visits)) * A0);
        patients.DET(i) = alpha^r * beta^s * delta^(s-1) * gamma(r+study_visits+1) * ...
            kummerU(s,s,delta*(beta+obs_time)) / ...
            (gamma(r)*(alpha+obs_time)^(r+study_visits+1) * L);
        
        % CLV !!!
        
        patients.CLV(i) = patients.DET(i) * patients.exp_sp(i);
        
    end
    
end

clear cohort study_visits  obs_time lst_rpt avg_purchase r alpha s beta ...
    p q gam F0 F1 F2 A0 pactive delta



% % compute average E[Y(t)|p1x,tx,T] and average actual number of 
% % transactions in the second 39 weeks for each level of p1x
% ce_act = zeros(max(p1x)+1,1);
% ce_est = zeros(max(p1x)+1,1);
% np1x = zeros(max(p1x)+1,1);
% for y = unique(p1x)'
%     isx = find(p1x==y);
%     np1x(y+1) = length(isx);
%     ce_act(y+1) = sum(p2x(isx))/np1x(y+1);
%     ce_est(y+1) = sum(ce(isx))/np1x(y+1);
% end
% 
% % create right-censored version for plot
% censor = 12;  % right-censor at 12+
% denom = sum(np1x(censor+1:length(np1x)));
% 
% ce_act_cen = ce_act(1:censor);
% ce_act_cen(censor+1) = (np1x(censor+1:length(np1x))'...
%     *ce_act(censor+1:length(np1x)))/denom;
% 
% ce_est_cen = ce_est(1:censor);
% ce_est_cen(censor+1) = (np1x(censor+1:length(np1x))'...
%     *ce_est(censor+1:length(np1x)))/denom;
% 
% subplot(2,2,4)
% 
% plot([0:censor],ce_act_cen,'k',[0:censor],ce_est_cen,'kp--');
% legend('Actual','Pareto/NBD','Location','SouthEast');
% xlabel(strcat('Transactions in Study Period (',string(weeks_in_study),' weeks)')); 
% ylabel(strcat('Average Transactions in Holdout Period (',string(weeks_in_study),' weeks)'));
% xticks([0 1 2 3 4 5 6 7 8 9 10 11 12])
% xticklabels({ ' 0', ' 1', ' 2', ' 3', ' 4', ' 5', ' 6', '7', '8', '9', '10', '11', '12+' });
% grid on
% grid minor
% % print -depsc 'ce_plot.eps'
% 
% clear tmp1 tmp2 y isx
