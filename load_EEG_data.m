opts = detectImportOptions('Data/EEG_raw.csv', 'NumHeaderLines', 0);
opts.VariableNames = {'Fp1','Fp2','F3','F4','F7','F8','T3','T4','C3','C4','T5','T6','P3','P4','O1','O2','Fz','Cz','Pz'};
T = readtable('Data/EEG_raw.csv', opts);

save('Data/EEG_data.mat', 'T');  % Save the data and column names