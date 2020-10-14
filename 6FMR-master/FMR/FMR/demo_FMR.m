clc;clear;
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
for di = 1:10
    path = '/home/zpp/NewDisk/2020/2020TKDE/major_revision/AllMethod/0datasets/smalldata/';
    datapath = [path, dataname{di},'.mat'];
    disp(dataname{di});
    dlmwrite('result0929.txt', dataname{di} ,'-append','delimiter','\t','newline','pc');

    f = load(datapath);
    X=f.X;
    Y=f.Y;
    for iv = 1:length(X)
        X{iv} = X{iv}';
    end
%% flexible representation learning method
fprintf('Latent representation multiview subspace clustering\n');
numClust = size(unique(Y),1);

%%   
lambda1 = 1e-4;   lambda2 = 1e-8;  epsilon =  -5e10; alpha = 100; K = 200; %yale
%lambda1 = 1e-4;   lambda2 = 1e-8;  epsilon =  -5e8; alpha = 100; K = 300; %ORL
%lambda1 = 1e-3;   lambda2 = 1e-6;  epsilon =  -5e8; alpha = 10; K = 200; %BBCSport
%lambda1 = 1e-4;   lambda2 = 1e-7;  epsilon =  -5e5; alpha = 10; K = 800; %Notting Hill
%lambda1 = 1e-4;   lambda2 = 1e-6;  epsilon =  -5e7; alpha = 100; K = 200; %MSRC-v1

lambda1 = 1;   lambda2 = 0.001;       % MSRCV1
% lambda1 = 1;   lambda2 = 1e-9;         % Wiki
% lambda1 = 0.01;   lambda2 = 0.0001; % Caltech 7
% lambda1 = 1e-5;   lambda2 = 1e-9;     % Caltech 20
for it = 1:1
    tic;
    [nmi, ACC, f, RI, H] = HSIC(X , Y, numClust,lambda1,lambda2,K,epsilon,alpha);           
    t(it)= toc;
    result(it,:) = [nmi, ACC, f, RI];
end
meanresult= mean(result,1);
meant = mean(t);
dlmwrite('result0929.txt',[lambda1,lambda2 meanresult meant],'-append','delimiter','\t','newline','pc');
 

end

