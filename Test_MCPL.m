function  Test_MCPL(dataname, X, Y )
% Input£º     X£ºdata matrix
%             Y£ºtruth label
% Output:     result£ºclustering result

% CESC	        5	50  500  5  0.1  0.0001				
% HW2sources	24	100  100  1  0.1  100				
% MSRC_v1	    15	0.005  100  50  0.1  5				
% MSRC	        8	1000  1000  0.05  0.1  0.0005  				
% Prokaryotic	35	5  50  1  50  5				
% Reuters	    5	1000  10  0.001  1  0.0001 				
% WikipediaArticles	99	0.5  10  0.05  1  1				
% Youtube	    5	5000  100  0.0005  0.5  0.00005				
% YaleB	        6	0.0005  1  0.0001  100  0.00005				

k = 8;
alpha=1000;
lambda=1000;
beta=0.05;
gamma=0.1;
miu=0.0005;

row1=string(zeros(1,13));
row1(1)='alpha';
row1(2)='lambda';
row1(3)='beta';
row1(4)='gamma';
row1(5)='miu';
row1(6)='ACC';
row1(7)='NMI';
row1(8)='Purity';
row1(9)='Fscore';
row1(10)='Precision';
row1(11)='Recall';
row1(12)='AR';
row1(13)='Entropy';

resultfile = strcat('result/',dataname,'.csv'); % record result(ACC NMI Purity...)
objfile = strcat('result/',dataname,'Obj.csv'); % record objvalue

f1 = fopen(resultfile, 'w+', 'n', 'utf8');      % create table title
fprintf(f1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',row1(1,1),row1(1,2),row1(1,3),row1(1,4),row1(1,5),row1(1,6),row1(1,7),row1(1,8),row1(1,9),row1(1,10),row1(1,11),row1(1,12),row1(1,13));
fclose(f1);

f2 = fopen(objfile, 'w+', 'n', 'utf8');         % create table title
fprintf(f2, '%s,%s,%s,%s,%s,%s\n', row1(1,1),row1(1,2),row1(1,3),row1(1,4),row1(1,5),'obj');
fclose(f2);

f3 = fopen(resultfile, 'a', 'n', 'utf8');       % add running results to file

V = size(X,2);      % view number
c = max(Y);         % class number
n = length(Y);      % object number

tic;
F = MCPL(dataname,X,V,n,c,k,alpha,lambda,beta,gamma,miu);

% Repeat kmeans 20 times.
for idx = 1:20
    pre_Y = kmeans(F,c,'maxiter',100,'replicates',10,'emptyaction','singleton');
    Res(idx,:) = Clustering8Measure(Y, pre_Y);
end
res = [mean(Res);std(Res)];
disp(res);

% write mean data to file
row = [alpha,lambda,beta,gamma,miu,res(1,1),res(1,2),res(1,3),res(1,4),res(1,5),res(1,6),res(1,7),res(1,8)];
fprintf(f3, '%.8f,%.8f,%.8f,%.8f,%.8f,%f,%f,%f,%f,%f,%f,%f,%f\n', row(1),row(2),row(3),row(4),row(5),row(6),row(7),row(8),row(9),row(10),row(11),row(12),row(13));

% write std data to file
fprintf(f3, '%s,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f\n',' ',' ',' ',' ',' ',res(2,1),res(2,2),res(2,3),res(2,4),res(2,5),res(2,6),res(2,7),res(2,8));

toc;

fclose(f3);

end

