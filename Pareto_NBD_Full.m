%% Basic assumptions for calibration duration and cohort size

weeks_in_cohort = 10;
weeks_in_study = 80;
weeks_total = 80;
% cohort_n = floor((80-weeks_in_study)/weeks_in_cohort +1);
cohort_n = 1;
charts_hor = 1;
charts_ver = 1;
first_date = min(patients.first_visit);

set(0,'DefaultFigureWindowStyle','docked');

cumrptslsfig = figure('Name','Cumulative repeat sales');
rptslsfig = figure('Name','Repeat sales');
spendmodel = figure('Name','Spend Model');
pareto_nbd_ce = figure('Name','Pareto/NBD Conditional Expectations');

CLVheatmap = figure('Name','CLV Heatmap');
ScatterCLV = figure('Name','3D Scatter CLV');

DETchart = figure('Name','DET Surface Plot');


hold

%% Building patient study set

patients.cohort = ceil((patients.first_visit-first_date+1)/(weeks_in_cohort*7));

for i = 1:size(patients)
    if patients.cohort(i) <= cohort_n
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
        
        patients.study_visits(i) = height(visits_of_patient_in_study)-1;
        patients.holdout_visits(i) = height(visits_of_patient_in_holdout);
        patients.observation_time(i) = (study_end-patients.first_visit(i) )/7;
        if patients.study_visits(i) ~= 0
            patients(i,'last_repeat') = num2cell((visits_of_patient_in_study.date_completed_serial( ...
                height(visits_of_patient_in_study)) - patients.first_visit(i)  )/7);
        end
    end
end

visits = innerjoin(visits,patients,'RightVariables',{'cohort','first_visit'});
visits.Properties.VariableNames([8 9]) = {'cohort' 'first_visit'};
visits.week = ceil((visits.date_completed_serial-first_date+1)/7);
visits.is_first = visits.date_completed_serial== visits.first_visit;

clear i visits_of_patient study_begin study_end holdout_end studyCondition visits_of_patient_in_study ...
    holdoutCondition visits_of_patient_in_holdout;

%% Building cumulative purchases variables for each cohort

for i = 1:max(patients.cohort)
    cohortCondition = visits.cohort == i & visits.is_first ==0;
    cohort_visits = accumarray(table2array(visits(cohortCondition,'week')), 1);
    weekly_visits(:,i) = cohort_visits(:,1);
    weekly_cumul_visits(:,i) = cumsum(cohort_visits(:,1));
end

clear i cohortCondition cohort_visits

%% Log likelihood minimzation parameter determination
% Requires pareto_nbd_ll.m and gamma_gamma_ll.m
%

global p1x p2x tx T zbar
cohort_results = table;

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
%%
% Script to compute conditional expectations
%


for i = 1:size(patients)
    cohort = patients.cohort(i);
    study_visits = patients.study_visits(i);
    obs_time = patients.observation_time(i);
    lst_rpt = patients.last_repeat(i);
    avg_purchase = patients.mean_rev(i);
    
    if cohort <=cohort_n
        
        r = cohort_results.r(cohort);
        alpha = cohort_results.alpha(cohort);
        s = cohort_results.s(cohort);
        beta = cohort_results.beta(cohort);
        p = cohort_results.p(cohort);
        q = cohort_results.q(cohort);
        gam = cohort_results.gamma(cohort);
        
        % period for which conditional expectations are to be computed 
        holdout_end=  ((cohort-1) * weeks_in_cohort + weeks_in_study*2);
        study_end =  ((cohort-1) * weeks_in_cohort + weeks_in_study) ;
        
        t = min(holdout_end,weeks_total)-study_end;  
        
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

clear i cohort study_visits  obs_time lst_rpt avg_purchase r alpha s beta ...
    p q gam F0 F1 F2 A0 pactive delta L



%% Script to create the cumulative and incremental tracking plots, 
%  and the spend model plots
% 
% Also compute average E[Y(t)|p1x,tx,T] and average actual number of
% transactions in the second 39 weeks for each level of p1x



days_in_cohort = weeks_in_cohort * 7;
days_in_study = weeks_in_study * 7;

x_axis_max = 80;
y_axis_max_rpt=700;
y_axis_max_cum=20000;
x_spend = 1:500;

for i = 1:cohort_n
    p1x = table2array(patients(patients.cohort == i,'study_visits'));
    tx = table2array(patients(patients.cohort == i,'last_repeat'));
    T = table2array(patients(patients.cohort == i,'observation_time'));
    p2x = table2array(patients(patients.cohort == i,'holdout_visits'));
    zbar = table2array(patients(patients.cohort == i,'mean_rev'));
    ce = table2array(patients(patients.cohort == i,'exp_tr'));

    % determine cohort size by day of trial
    ns = [];
    for j = 1:days_in_cohort
        ns(j) = sum( round(T,2) == round((days_in_study+1-j)/7,2));
    end
    
    % generate sales cumulative forecast
    endwk = 80 ; % MODIFIED FOR FULL
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
    subplot(charts_hor,charts_ver,i)
    
    plot(1:endwk,weekly_cumul_visits(1:endwk,i),'k', ...
        1+chart_offset:endwk+chart_offset,cumrptsls,'b--');
    axis([0 x_axis_max 0 y_axis_max_cum]);
    xlabel('Week'); ylabel('Cumulative Repeat Transactions');
    legend('Actual','Model','Location','SouthEast');
    title(['Cohort ',num2str(i), ...
        ' - Frequency: ', num2str(round(cohort_results.alpha(i)/cohort_results.r(i),1)), ...
        ' - Lifetime: ', num2str(round(cohort_results.beta(i)/cohort_results.s(i),1))]);
    % print -depsc 'cumrptsls.eps'
    
    % create tracking plot of weekly repeat sales (pred. vs actual)
    figure(rptslsfig);
    subplot(charts_hor,charts_ver,i)
    
    incrptsls = [ cumrptsls(1) diff(cumrptsls) ];
    plot(1:endwk,weekly_visits(1:endwk,i),'k', ...
        1+chart_offset:endwk+chart_offset,incrptsls,'b--');
    axis([0 x_axis_max 0 y_axis_max_rpt]);
    xlabel('Week'); ylabel('Weekly Repeat Transactions');
    legend('Actual','Model','Location','SouthEast');
    title(['Cohort ',num2str(i), ...
        ' - Frequency: ', num2str(round(cohort_results.alpha(i)/cohort_results.r(i),1)), ...
        ' - Lifetime: ', num2str(round(cohort_results.beta(i)/cohort_results.s(i),1))]);
    % print -depsc 'incrptsls.eps'
    
    % create  plot of spending model
    figure(spendmodel);
    
    subplot(charts_hor,charts_ver,i)
    histogram(zbar,'Normalization','pdf', 'BinLimits',[0,500],'BinWidth',10);
    hold on;
    plot(x_spend,(cohort_results.p(i).*cohort_results.gamma(i)).^cohort_results.q(i) ...
        .*x_spend.^(-cohort_results.q(i)-1).*exp(-cohort_results.p(i).*cohort_results.gamma(i)./x_spend) ...
        ./gamma(cohort_results.q(i)));
    title(['Cohort ',num2str(i), ...
        ' - Mode: $', num2str(round(cohort_results.mode_spend(i),0)), ...
        ' - Average: $', num2str(round(cohort_results.average_spend(i),0))]);
    
    % create plots of conditional expectations
    ce_act = zeros(max(p1x)+1,1);
    ce_est = zeros(max(p1x)+1,1);
    np1x = zeros(max(p1x)+1,1);
    for y = unique(p1x)'
        isx = find(p1x==y);
        np1x(y+1) = length(isx);
        ce_act(y+1) = sum(p2x(isx))/np1x(y+1);
        ce_est(y+1) = sum(ce(isx))/np1x(y+1);
    end
    
    % create right-censored version for plot
%     censor = 12;  % right-censor at 12+
%     denom = sum(np1x(censor+1:length(np1x)));
%     
%     ce_act_cen = ce_act(1:censor);
%     ce_est_cen = ce_est(1:censor);
%     
%     ce_act_cen(censor+1) = (np1x(censor+1:length(np1x))*ce_act(censor+1:length(np1x)))/denom;
%     ce_est_cen(censor+1) = (np1x(censor+1:length(np1x))*ce_est(censor+1:length(np1x)))/denom;

    grouped = [[5;6],[7;9],[10;15],[16;20], [21;length(np1x)]];
    ce_act_cen = ce_act(1:4);
    ce_est_cen = ce_est(1:4);
    for j=1:length(grouped)
        denom = sum(np1x(grouped(:,j)));
        ce_act_cen(4+j) = sum(np1x(grouped(:,j)).*ce_act(grouped(:,j)))/denom;
        ce_est_cen(4+j) = sum(np1x(grouped(:,j)).*ce_est(grouped(:,j)))/denom;
    end
            % period for which conditional expectations are to be computed 
        holdout_end=  ((i-1) * weeks_in_cohort + weeks_in_study*2);
        study_end =  ((i-1) * weeks_in_cohort + weeks_in_study) ;
        weeks_in_holdout = min(holdout_end,weeks_total)-study_end;
    
    figure(pareto_nbd_ce);
    subplot(charts_hor,charts_ver,i)
    
    plot([0:4+length(grouped)-1],ce_act_cen,'kp-',[0:4+length(grouped)-1],ce_est_cen,'bp--');

    xlabel(strcat('During Study Period (',string(weeks_in_study),' weeks)'));
    ylabel(strcat('During Holdout Period (',string(weeks_in_holdout),' weeks)'));
    xticks([0:4+length(grouped)-1])
    xticklabels({ ' 0', ' 1', ' 2', ' 3', ' 4-5', ' 6-8', '9-14', '15-19', '20+'});
    grid on
    grid minor
    
    yyaxis right
    bar(histc(ce,[0 1 2 3 5 8 14 19 Inf])/length(ce),'FaceColor',[0.1 0.8 0.3], ...
        'XData',[0:4+length(grouped)-1],'FaceAlpha',0.3);
    axis([0 8 0 1]);
    ylabel('Proportion of sub-group in population');
    
    legend('Actual','Model','Representation','Location','North');

end

hold off

clear ns tmp1 tmp2 tmpcumsls1 tmpcumsls2 chart_offset i j x_axis_max x_spend ...
    y_axis_max_rpt y_axis_max_cum cumrptsls incrptsls dailysls;

%%
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
% Script to compute RFM Matrix with 3 groups
%

groups = 3;
group_results = table;

cohort_data = patients(patients.study_visits==0 | patients.cohort>cohort_n, ...
    {'study_visits','last_repeat','mean_rev'});
cohort_data.study_visits_rankings=zeros(size(cohort_data.study_visits));
cohort_data.last_repeat_rankings=zeros(size(cohort_data.study_visits));
cohort_data.mean_rev_rankings=zeros(size(cohort_data.study_visits));
cohort_data.study_visits_group=zeros(size(cohort_data.study_visits));
cohort_data.last_repeat_group=zeros(size(cohort_data.study_visits));
cohort_data.mean_rev_group=zeros(size(cohort_data.study_visits));

group_results = [group_results;cohort_data];


for i = 1:cohort_n
    cohort_data = patients(patients.cohort==i & patients.study_visits>0, ...
        {'study_visits','last_repeat','mean_rev'});
    
    [~,~,cohort_data.study_visits_rankings]=unique(cohort_data.study_visits);
    [~,~,cohort_data.last_repeat_rankings]=unique(cohort_data.last_repeat);
    [~,~,cohort_data.mean_rev_rankings]=unique(cohort_data.mean_rev);
    
    q = quantile(cohort_data.study_visits, [1:groups-1]/groups);
    cohort_data.study_visits_group = sum(cohort_data.study_visits > q,2)+1;
    q = quantile(cohort_data.last_repeat, [1:groups-1]/groups);
    cohort_data.last_repeat_group = sum(cohort_data.last_repeat > q,2)+1;
    q = quantile(cohort_data.mean_rev, [1:groups-1]/groups);
    cohort_data.mean_rev_group = sum(cohort_data.mean_rev > q,2)+1;
    
    group_results = [group_results;cohort_data];
    
end

patients_with_groups = join(patients, group_results, 'Keys','RowNames', ...
    'RightVariables',{'study_visits_rankings','last_repeat_rankings','mean_rev_rankings', ...
    'study_visits_group','last_repeat_group','mean_rev_group'});

RFM = grpstats(patients_with_groups, {'cohort','study_visits_group','last_repeat_group','mean_rev_group'},'sum', ...
    'Datavars',{'CLV'});

RFM.sum_CLV = RFM.sum_CLV/1000;

figure(CLVheatmap);
clf(CLVheatmap);
load cmap.mat;

for i = 1:cohort_n
    
    subplot(groups,cohort_n,i)
    heatmap(RFM((RFM.cohort==i & RFM.mean_rev_group==1) | (RFM.cohort==i & RFM.mean_rev_group==0) , ...
        {'study_visits_group','last_repeat_group','sum_CLV'}), ...
        'study_visits_group','last_repeat_group','ColorVariable','sum_CLV', ...
        'ColorMethod', 'sum', 'Colormap',mycmap,'ColorLimits',[0 500], ...
        'CellLabelFormat','$%4.0f k');
    xlabel('Purchases (3 is more)'); ylabel('Recency (3 is more recent)');
    title(['Cohort ',num2str(i),': Low "M" CLV']);
    
    subplot(groups,cohort_n,cohort_n+i)
    heatmap(RFM((RFM.cohort==i & RFM.mean_rev_group==2) | (RFM.cohort==i & RFM.mean_rev_group==0) , ...
        {'study_visits_group','last_repeat_group','sum_CLV'}), ...
        'study_visits_group','last_repeat_group','ColorVariable','sum_CLV', ...
        'ColorMethod', 'sum', 'Colormap',mycmap,'ColorLimits',[0 500], ...
        'CellLabelFormat','$%4.0f k');
    xlabel('Purchases (3 is more)'); ylabel('Recency (3 is more recent)');
    title(['Cohort ',num2str(i),': Mid "M" CLV']);
    
    subplot(groups,cohort_n,2*cohort_n+i)
    heatmap(RFM((RFM.cohort==i & RFM.mean_rev_group==3) | (RFM.cohort==i & RFM.mean_rev_group==0) , ...
        {'study_visits_group','last_repeat_group','sum_CLV'}), ...
        'study_visits_group','last_repeat_group','ColorVariable','sum_CLV', ...
        'ColorMethod', 'sum', 'Colormap',mycmap,'ColorLimits',[0 500], ...
        'CellLabelFormat','$%4.0f k');
    xlabel('Purchases (3 is more)'); ylabel('Recency (3 is more recent)');
    title(['Cohort ',num2str(i),': Top "M" CLV']);
    
end

clear i group_results cohot_data q

%%
% Script to compute RFM Matrix with 10 group for one cohort
%

figure(ScatterCLV);
clf(ScatterCLV);

z = table2array(patients(patients.cohort<=4,'mean_rev'));
y = table2array(patients(patients.cohort<=4,'study_visits'));
x = table2array(patients(patients.cohort<=4,'last_repeat'));
c = table2array(patients(patients.cohort<=4,'CLV'));
scatter3(x,y,z,30,c,'filled')
zlabel('Average Purchase size (log)')
ylabel('Number of purchases during study period (log)')
xlabel('Recency of the last purchase')
colorbar('Location', 'EastOutside')

% 
%     heatmap(RFM(RFM.cohort==cohort , {'study_visits_group','last_repeat_group','sum_CLV'}), ...
%         'study_visits_group','last_repeat_group','ColorVariable','sum_CLV', ...
%         'ColorMethod', 'sum', 'Colormap',jet, ...
%         'ColorLimits',[1 100000]);
%     xlabel('Purchases (5 is more)'); ylabel('Recency (5 is more recent)');
%     title(['Cohort ',num2str(cohort)]);
    
clear i group_results cohot_data q z y x c

%%
% Script to compute DET chart
%



cohort = 1
obs_time = 80

r = cohort_results.r(cohort);
alpha = cohort_results.alpha(cohort);
s = cohort_results.s(cohort);
beta = cohort_results.beta(cohort);

freq_points = linspace(0, 20, 100);
rec_points = linspace(0, 80, 100);

[study_visits, lst_rpt] = meshgrid(freq_points, rec_points);

delta = log(1+0.1)/52;

maxab = max(alpha,beta);
absab = abs(alpha-beta);
param2 = s+1;
if alpha < beta
    param2 = r+study_visits;
end

F0 = (alpha+obs_time).^(r+study_visits)*(beta+obs_time)^s;
F1 = h2f1(r+s+study_visits,param2,r+s+study_visits+1,absab./(maxab+lst_rpt))./...
    ((maxab+lst_rpt).^(r+s+study_visits));
F2 = h2f1(r+s+study_visits,param2,r+s+study_visits+1,absab./(maxab+obs_time))./...
    ((maxab+obs_time).^(r+s+study_visits));
A0 = F1-F2;
L = gamma(r+study_visits)*alpha^r*beta^s./gamma(r) .* ...
    (1./F0 + (s./(r+s+study_visits)) .* A0);
DET = alpha^r * beta^s * delta^(s-1) * gamma(r+study_visits+1) * ...
    kummerU(s,s,delta*(beta+obs_time)) ./ ...
    (gamma(r)*(alpha+obs_time).^(r+study_visits+1) .* L);

figure(DETchart);
subplot(1,2,1);
surf( lst_rpt,study_visits, DET)
zlabel('Discounted Expected Transactions')
ylabel('Frequency')
xlabel('Recency')
colorbar('Location', 'EastOutside')

subplot(1,2,2);
contour( lst_rpt,study_visits, DET, [0 1 2 3 4 5 7 10 15 20 30 40 50],'ShowText','on', ...
    'Fill','on', 'LineStyle','-', 'LineWidth', 1, 'LineColor',[0.3 0.3 0.3])
ylabel('Frequency')
xlabel('Recency')

clear cohort obs_time r alpha s beta freq_points rec_points delta maxab absab ...
    param2 F0 F1 F2 A0 L DET