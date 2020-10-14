% This code runs the LMVSC algorithm on Caltech7 dataset and records its
% performance. The results are output in *.txt format.

% Notice: The dataset is organized in a cell array with each element being
% a view. Each view is represented by a matrix, each row of which is a
% sample.

% The core of LMVSC is encapsulated in an independent matlab function.
% Visit lmv.m directly, if you want to learn the details of its
% implementation.
clc;
clear;
addpath('./measure');
addpath('./datasets');
dataname = {'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
for di = 1 : 10
    path = '../0datasets/smalldata/';
load([path, dataname{di},'.mat']);
dlmwrite('result-0929.txt',[dataname{di}],'-append','delimiter','\t','newline','pc');

y=Y;
nv=length(X);
ns=length(unique(y));

% Parameter 1: number of anchors (tunable)
numanchor=[ns 50 100];

% Parameter 2: alpha (tunable)
alpha=[.01 1 100]; 

for j=1:length(numanchor)
    
    % Perform K-Means on each view
    parfor i=1:nv
        rand('twister',5489);
        [~, H{i}] = litekmeans(X{i},numanchor(j),'MaxIter', 100,'Replicates',10);
    end

    for i=1:length(alpha)
        fprintf('params:\t numanchor=%d\t\talpha=%f\n',numanchor(j),alpha(i));
        
        for repi = 1:1
            tic;
            % Core part of this code (LMVSC)
            [F,ids] = lmv(X',y,H,alpha(i));

            % Performance evaluation of clustering result
            result(repi,:)=Clustering8Measure(ids,y);
            % result = [Fscore nmi ACC Purity];
            time(repi)=toc;
        end
        meanresult = mean(result,1);
        t=mean(time);
        fprintf('result:\t%12.6f %12.6f %12.6f %12.6f %12.6f\n',[meanresult t]);
        
        % Write the evaluation results to a text file
        dlmwrite('result-0929.txt',[alpha(i) numanchor(j) meanresult t],'-append','delimiter','\t','newline','pc');
        
    end
end
end