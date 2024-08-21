function [results_mean_before, results_sd_before, results_mean_after, results_sd_after] = mean_sd(wreg, flag_transf, sp , ftype, fdomain, T0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%wreg:          WHO region for which the data is analysed    
%flag_transf:   flag for data transformation
%               0=no transf. is applied to the data
%               1=transf. from GP Nason, Scientific Reports, 2020
%               2=transf. which is usually applied to daily stock markets indices at closing time
%               3=first order differences (computed in time domain)
%sp             'day'=daily data, 'week'=weekly data
%ftype          'stop'=stop band filter, 'high'=high pass filter
%fdomain        'time'=filtering in time domain, 'freq'=filtering in freq. domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the corresponding .mat file.
fina = strcat('Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, '_rect', '.mat');
load(fina);
[len,T] = size(t_before);
%Create two new cells to store the mean and standard deviation (before filtering).
results_mean_before = zeros(1, T-1);
results_sd_before = zeros(1, T-1);
nearest_number = cell(len, 1);
for i = 1:len
    % Retrieve the data from the second column of the current row.
    data = cell2mat(t_before{i, 2});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0));
    nearest_number{i, 1} = data(idx)-T0;
end
results_mean_before(1, 1) = mean(cell2mat(nearest_number(:, 1)));
results_sd_before(1, 1) = std(cell2mat(nearest_number(:, 1)));

for i = 1:len
    % Retrieve the data from the third column of the current row.
    data = cell2mat(t_before{i, 3});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/2));
    nearest_number{i, 2} = data(idx)-T0/2;
end
results_mean_before(1, 2) = mean(cell2mat(nearest_number(:, 2)));
results_sd_before(1, 2) = std(cell2mat(nearest_number(:, 2)));

for i = 1:len
    % Retrieve the data from the fourth column of the current row.
    data = cell2mat(t_before{i, 4});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/3));
    nearest_number{i, 3} = data(idx)-T0/3;
end
results_mean_before(1, 3) = mean(cell2mat(nearest_number(:, 3)));
results_sd_before(1, 3) = std(cell2mat(nearest_number(:, 3)));


for i = 1:len
    % Retrieve the data from the fifthcolumn of the current row.
    data = cell2mat(t_before{i, 5});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/4));
    nearest_number{i, 4} = data(idx)-T0/4;
end
results_mean_before(1, 4) = mean(cell2mat(nearest_number(:, 4)));
results_sd_before(1, 4) = std(cell2mat(nearest_number(:, 4)));

for i = 1:len
    % Retrieve the data from the sixth column of the current row.
    data = cell2mat(t_before{i, 6});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/5));
    nearest_number{i, 5} = data(idx)-T0/5;
end
results_mean_before(1, 5) = mean(cell2mat(nearest_number(:, 5)));
results_sd_before(1, 5) = std(cell2mat(nearest_number(:, 5)));

%Create two new cells to store the mean and standard deviation (after filtering).
[len,T] = size(t_after);
results_mean_after = zeros(1, T-1);
results_sd_after = zeros(1, T-1);
nearest_number = cell(len, 1);

for i = 1:len
    % Retrieve the data from the second column of the current row.
    data = cell2mat(t_after{i, 2});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0));
    nearest_number{i, 1} = data(idx)-T0;
end
results_mean_after(1, 1) = mean(cell2mat(nearest_number(:, 1)));
results_sd_after(1, 1) = std(cell2mat(nearest_number(:, 1)));

for i = 1:len
    % Retrieve the data from the third column of the current row.
    data = cell2mat(t_after{i, 3});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/2));
    nearest_number{i, 2} = data(idx)-T0/2;
end
results_mean_after(1, 2) = mean(cell2mat(nearest_number(:, 2)));
results_sd_after(1, 2) = std(cell2mat(nearest_number(:, 2)));

for i = 1:len
    % Retrieve the data from the fourth column of the current row.
    data = cell2mat(t_after{i, 4});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/3));
    nearest_number{i, 3} = data(idx) - T0/3;
end
results_mean_after(1, 3) = mean(cell2mat(nearest_number(:, 3)));
results_sd_after(1, 3) = std(cell2mat(nearest_number(:, 3)));

for i = 1:len
    % Retrieve the data from the fifth column of the current row.
    data = cell2mat(t_after{i, 5});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/4));
    nearest_number{i, 4} = data(idx)-T0/4;
end
results_mean_after(1, 4) = mean(cell2mat(nearest_number(:, 4)));
results_sd_after(1, 4) = std(cell2mat(nearest_number(:, 4)));

for i = 1:len
    % Retrieve the data from the sixth column of the current row.
    data = cell2mat(t_after{i, 6});

    % Find the number closest to T0.
    [~, idx] = min(abs(data - T0/5));
    nearest_number{i, 5} = data(idx)-T0/5;
end
results_mean_after(1, 5) = mean(cell2mat(nearest_number(:, 5)));
results_sd_after(1, 5) = std(cell2mat(nearest_number(:, 5)));
end