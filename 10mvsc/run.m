clear all;
clc;
tic;
% load bbcsport_seg14of1.mat  (data,labels)
addpath('./measure');
dataname = {'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};

for di =1:10
    daname = ['../0datasets/smalldata/',dataname{di}]
f=load(daname);
dlmwrite('result-0929.txt',[dataname{di}],'-append','delimiter','\t','newline','pc');
X = f.X;
Y = f.Y;
C = length(unique(Y)); %% Clusters of samples
N = length(Y);

maxIters = 100;

options = [];
options.NeighborMode = 'KNN';
options.k = 3;
% options.WeightMode = 'Binary';
Zall = zeros(N,N);
for iv = 1:length(X)
    X{iv} = NormalizeFea(X{iv},0);
    Z = constructW(X{iv},options);% input N*d dim X
    Z_ini{iv} = full(Z);
    Zall = Zall + Z_ini{iv};
end
    D   = diag(sum(Zall,2));
    D2 = diag(1./sqrt(diag(D)+eps));
    L   = D2 - Zall;
    Lnor = D2*L*D2;
    Lnor = (Lnor + Lnor')/2;

    options.tol = 1e-8;
    options.maxit = 30000;
    [Finit, ~] = eigs(Lnor, C, 'sa', options);



lambda1 = 10.^[-5:2:5];
lambda2 = 10.^[-5:2:5];
alpha = 10.^[1:1:2];

lambda1 = 10^-5;
lambda2 = 10^-3;
alpha=10;

for i = 1:length(lambda1)
    for j = 1:length(lambda2)
        for k = 1:length(alpha)
            disp(['lmd1 = ',num2str(lambda1(i)),'  lmd2 = ',num2str(lambda2(j)),...
                '  alpha = ',num2str(alpha(k))]);
            tic;
            [Z, F, E] = mvcsolver( X, C, Finit, lambda1(i), lambda2(j), alpha(k), maxIters );
            Ypred = kmeans(F,C);
            result = Clustering8Measure(Y, Ypred);
            t=toc;
            finres = [lambda1(i), lambda2(j), alpha(k), result, t];
            % disp(finres);
            dlmwrite('result-0929.txt', finres ,'-append','delimiter','\t','newline','pc');
        end
    end
end
%X: 1 by n cell, each cell corresponds to each view
%C: number of classes
%alpha: defined in the paper, tune it according to the performance

toc;
end
