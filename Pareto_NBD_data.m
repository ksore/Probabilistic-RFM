
formatSpec = '%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%s%s%s%[^\n\r]';

fileID = fopen('orders_subset.csv','r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false, ...
    'HeaderLines',1);

fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,6,7]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{8} = datetime(dataArray{8}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{8} = cellfun(@(x) x(2:end-1), dataArray{8}, 'UniformOutput', false);
        dates{8} = datetime(dataArray{8}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
    catch
        dates{8} = repmat(datetime([NaN NaN NaN]), size(dataArray{8}));
    end
end

anyBlankDates = cellfun(@isempty, dataArray{8});
anyInvalidDates = isnan(dates{8}.Hour) - anyBlankDates;
dates = dates(:,8);

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [1,2,6,7]);
rawCellColumns = raw(:, [3,4,5]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
orders = table;
orders.order_number = cell2mat(rawNumericColumns(:, 1));
orders.patient_number = cell2mat(rawNumericColumns(:, 2));
orders.patient_group = rawCellColumns(:, 1);
orders.location_name = rawCellColumns(:, 2);
orders.order_source = rawCellColumns(:, 3);
orders.net_revenue = cell2mat(rawNumericColumns(:, 3));
orders.cost = cell2mat(rawNumericColumns(:, 4));
orders.date_completed = dates{:, 1};

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% orders.date_completed=datenum(orders.date_completed);

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me dates blankDates anyBlankDates invalidDates anyInvalidDates rawNumericColumns rawCellColumns R;

%% Apply filters for Pareto-NBD
orders_no_anon = orders(~any(isnan(orders.patient_number),2),:);

visits = grpstats(orders_no_anon, {'patient_number','date_completed'},'sum', ...
    'Datavars',{'net_revenue','cost'});