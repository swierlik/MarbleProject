opts = detectImportOptions('Data/dst_labels.csv', 'NumHeaderLines', 1);
opts.VariableNames = {'period', 'timedelta', 'dst'};
T = readtable('Data/dst_labels.csv', opts);

save('Data/MAG_data.mat', 'T');  % Save the data and column names