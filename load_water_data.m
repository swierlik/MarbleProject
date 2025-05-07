filename = 'Data/water_data_45.txt';

% Read file, skipping comment lines
opts = delimitedTextImportOptions("NumVariables", 20, "Delimiter", " ", ...
    "ConsecutiveDelimitersRule", "join", "LeadingDelimitersRule", "ignore");

opts.DataLines = [3 Inf]; % Skip the first two header lines
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Define variable names manually (from header lines)
opts.VariableNames = {'YY','MM','DD','hh','mm','WDIR','WSPD','GST','WVHT','DPD','APD','MWD',...
                      'PRES','ATMP','WTMP','DEWP','VIS','PTDY','TIDE'};

opts.VariableTypes = repmat("double", 1, 19); % default to double

% Treat 'MM' as missing
opts = setvaropts(opts, opts.VariableNames, "TreatAsMissing", "MM");

% Read as table
T = readtable(filename, opts);

% Convert to struct array (if desired)
dataStruct = table2struct(T);

% Optional: view one row
disp(dataStruct(1))

save('Data/weather_table_45day.mat', 'T');          % Saves table