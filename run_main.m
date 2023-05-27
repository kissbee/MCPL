
clc;
clear memory;
clear all;
warning('off')
addpath('data')
addpath('func')

% make dir to save results
file_path_name_ = 'result';
if exist(file_path_name_,'file')==0   %if not exist,make one
    mkdir(file_path_name_);
end

dataname= 'MSRC';  % select data name
load(dataname);
disp(['--test data:'  dataname ])

for i = 1:size(X,2)
    X{i} = X{i}';
end
Test_MCPL(dataname,X,Y); % test data dimension: dv * n


