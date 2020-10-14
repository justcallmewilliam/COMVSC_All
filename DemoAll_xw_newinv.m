clc;
clear;
%path = '.\'; % Windows
%datapath = '.\0datasets\'; % Windows

path = './'; % Ubuntu
datapath = './0datasets/'; % Ubuntu

addpath(genpath(path));
warning off;

datasetname = cell(6 , 1);

%% These datasets d<N
datasetname{1}='proteinFold';
datasetname{2}='Caltech101-7';

datasetname{3}='Caltech101-20'; % hzc
datasetname{4}='WikipediaArticles'; % jiyuan

datasetname{5}='yaleA_3view'; %
datasetname{6}='Handwritten_numerals'; %

Allresult = zeros(13,8);
% MethodOrder: Ours / FeatConcate / Co-regularized / MLRSSC / LMSC / 
%                                       RMKM / mPAC / FMR / GMC / LMVSC / PMSC / MVSC / 
% -8-Measure : Fscore Precision Recall nmi AR Entropy ACC Purity

runtimes = 1; 
%%
for di = 2:2  % Datasets
    dataName = datasetname{di};
    % load([datapath,'/',dataName,'.mat'],'X','Y'); % Windows
    load([datapath,'/',dataName,'.mat'],'X','Y');   % Ubuntu
    disp(['\n\n Current dataset  : ',dataName]);
    
    viewN = length(X);
    k = length(unique(Y));
    N = length(Y);
    
    for iv = 1:viewN
        X{iv} = mapstd(X{iv}',0,1); 
        data{iv} = X{iv}';
    end
    % X{i} is d * N
    % data{i} is N *d
    
    FV = cell(size(X));
    WV= cell(size(X));
    LV = cell(size(X));
    for i =1:viewN
        W = constructW_PKN(X{i}, 5, 1);  %   X : d*n  / 5 neighbors /  is symmetric
        D = diag(sum(W));
        L = D-W;
        [Fv, ~, ~]=eig1(L,k,0); 
        FV{i} = Fv;
        LV{i} = L;
        WV{i} = W;
    end

    %% **************************** Our Proposed Method ****************************
    lmd = 2.^[1:1:13]; %%2.^[7:2:13]
    gma = [1.1:0.1:3];
    para3=1;
    para4=1; % Fixed

    B=cell(size(X)); % B=inv(X'X+alpha*I)^(-1)
    idx = 1; 
    for i=1:length(lmd)
        for vnum=1:viewN 
            % B{vnum}= inv(X{vnum}'*X{vnum}+lmd(i)*eye(size(X{1},2)));
            I_d{vnum} = eye(size(X{vnum},1)); % X : d*N
            I_n{vnum} = eye(size(X{vnum},2));
            B{vnum} = (1/lmd(i))*(I_n{vnum}-(1/lmd(i))*X{vnum}'*inv(I_d{vnum}+(1/lmd(i))*X{vnum}*X{vnum}')*X{vnum});
        end
        for j=1:length(gma)
            tic;
            [ress,Fv,Fstar,Lv,Wv,Yres,obj]=CLOMV(X,Y,B,lmd(i),gma(j),para3,para4,FV,WV,LV);
            t(idx)=toc;
            disp(['Para1 = ',num2str(lmd(i)),'  Para2 = ',num2str(gma(j)), '  Time = ',num2str(t(idx))]);
            res_comvsc = ress(end,:);
            Result_11COMVSC(idx,:) = [lmd(i) gma(j) res_comvsc];
            
            ourFscore(i,j)  = res_comvsc(1); ourPrecision(i,j) = res_comvsc(2); ourRecall(i,j) = res_comvsc(3); 
            ourNmi(i,j) = res_comvsc(4);ourAR(i,j) = res_comvsc(5); ourEntropy(i,j) = res_comvsc(6); 
            ourACC(i,j) = res_comvsc(7); ourPurity(i,j) = res_comvsc(8);
            dlmwrite(['TxtResultOfOur_', dataName, '.txt'], Result_11COMVSC(idx,:),'-append','delimiter','\t','newline','pc');
            idx = idx+1;
        end
    end
    maxresult = max(Result_11COMVSC,[],1);
    Allresult(1,:) = maxresult(3:10);
    fprintf('- Finish Proposed Method \n ');

    save(['./' , dataName ,'_1013_Our.mat'], 't','Allresult', 'Result_11COMVSC',...
        'ourFscore', 'ourPrecision', 'ourRecall' ,'ourNmi', 'ourAR', 'ourEntropy', 'ourACC', 'ourPurity');


%     %% **************************** 0 Concatenate Features ****************************
%     allX = [];
%     for i=1:viewN % view num
%         allX = [allX data{i}];
%     end
% 
%     idx_cc = litekmeans(allX, k, 'MaxIter', 100, 'Replicates',20);
%     Allresult(2,:) = Clustering8Measure(Y, idx_cc);
%     fprintf('- Finish Concatenate Features \n ');
% 
%     %% **************************** 1Co-regularized ****************************
%     % 1.1multiview spectral (pairwise): more than 2 views
%     clear sigma;
%     for i = 1:viewN
%         sigma(i) = optSigma(data{i}'); % Input : N*d
%     end
%     lambda = [0.01 0.02 0.03 0.04 0.05]; 
%     numiter = 50;
%     for li =1:length(lambda)
%         % fprintf('1Co-regularized-Pairwise: Multiview spectral with 2 views \n ');
%         [Result_1Co_pw] = spectral_pairwise_multview(data,viewN,k,sigma,lambda(li),Y,numiter);
%     end
%     Allresult(3,:) = Result_1Co_pw;
% 
%     % 1.2multiview spectral (centroid): more than 2 views
%     % fprintf('1Co-regularized-Centroid:Multiview spectral with 2 views \n ');
%     lambda2 = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
%     [Result_1Co_ct] = spectral_centroid_multiview(data,viewN,k,sigma,lambda2,Y,numiter);
%     Allresult(4,:) = Result_1Co_ct;
%     fprintf('- Finish Co-regularized \n ');
% 
%     %% **************************** 2MLRSSC-Centroid ****************************
%     opts.mu = [10, 10^2, 10^3, 10^4];
%     lambda1 = [0.1 0.3 0.5 0.7 0.9];
%     lambda3 = [0.3 0.5 0.7 0.9];
%     Result_2MLRSSC = zeros(length(opts.mu)*length(lambda1)*length(lambda3), 3+8);
%     idx = 1;
%     for i=1:length(opts.mu)
%         for j = 1:length(lambda1)
%             for q = 1:length(lambda3)
%                 opts.lambda = [lambda1(j) (1-lambda1(j)) lambda3(q)];
%                 A = centroid_MLRSSC(data, opts); % joint affinity matrix % Input:N*d
%                 ress = spectral_clustering(A, k, Y);
%                 Result_2MLRSSC(idx,:) = [opts.mu(i), lambda1(j), lambda3(q), ress];
%                 idx = idx +1;
%             end
%         end
%     end
%     maxress = max(Result_2MLRSSC,[],1);
%     Allresult(5,:) = maxress(4:11);
%     fprintf('- Finish MLRSSC \n ');
%     
%     %% **************************** 3 LMSC ****************************
%     lambda = 10.^(-3 :2:3);    
%     K = 100; 
%     Result_3LMSC =  zeros(length(lambda), 8);
%     for i=1:length(lambda)
%     	[Result_3LMSC(i,:), H, cov_val] = LRMSC(X,Y,k,lambda(i), K); % Input:d*N
%     end
%     Allresult(6,:)= max(Result_3LMSC,[],1);
%  
%     
%     fprintf('- Finish LMSC \n ');
%     %% **************************** 4 RMKM ****************************
%     % init G0
%     Ik = eye(k);
%     randorder = randperm(size(Ik,1));
%     numceil = ceil(N/k);
%     largeG = repmat(Ik(randorder,:),numceil,1);
%     inG0= largeG(1:N,:);% N*k
% 
%     gma = linspace(1.1,0.2,20);
%     Result_4RMKM = zeros(length(gma), 9);
%     for i = 1:length(gma)
%             % Input X : d*N
%             [ outG0, outFCell, outAlpha, outAlpha_r, outObj, outNumIter ] = weighted_robust_multi_kmeans(X, gma(i), inG0);
%             [~,res_label]=max(outG0,[],2);
%             ress= Clustering8Measure(res_label,Y);
%             Result_4RMKM(i,:) = [gma(i), ress];
%     end
%     maxresult = max(Result_4RMKM,[],1);
%     Allresult(7,:) = maxresult(2:9);
%     fprintf('- Finish RMKM \n ');
%     %% **************************** 5 mPAC ****************************
%     para1=[10 20 30];
%     para2=[.0001 .001 .01];
%     para3=[.000001 .00001 .0001];
%     Result_5mPAC = zeros(length(para1)*length(para2)*length(para3), 3+8);
% 
%     Bn=cell(size(data));
%     idx = 1;
%     for i=1:length(para1)
%         for ii=1:viewN
%             Bn{ii}= inv(data{ii}*data{ii}'+para1(i)*eye(size(data{1},1)));
%         end
%         for j=1:length(para2)
%             for q=1:length(para3)
%                 disp(['para1=',num2str(para1(i)),'  Para2=',num2str(para2(j)),...
%                     '  Para3=',num2str(para3(q))]);
%             	[ress_mpac, FVret ,obj,Wv]=UMVSC(data,Y,Bn,para1(i),para2(j),para3(q),FV,WV,LV);
%                 Result_5mPAC(idx,:) = [para1(i), para2(j), para3(q), ress_mpac];
%                 idx = idx +1;
%             end
%         end
%     end
%     maxress = max(Result_5mPAC,[],1);
%     Allresult(8,:) = maxress(4:11);
%     fprintf('- Finish mPAC \n ');
% %  save(['./Result1012_' , dataName , '_mPAC.mat'], 'Allresult', 'Result_5mPAC');

%     %% **************************** 6 FMR ****************************
%     K = 200;
%     lambda1 = 10.^(-5 : 1 : 0);
%     lambda2 = 10.^(-10 : 1 : -2);
%     epsilon = -5e10;
%     alpha = 100;
% 
%     Result_6FMR = zeros(length(lambda1)*length(lambda2), 2+8);
%     
%     idx=1;
%     for i = 1:length(lambda1)
%         for j = 1:length(lambda2)
%             [ress_fmr , H] = HSIC(X , Y, k,lambda1(i),lambda2(j),K,epsilon,alpha);
%             Result_6FMR(idx,:) = [lambda1(i), lambda2(j), ress_fmr];
%             idx=idx+1;
%         end
%     end
%     maxress = max(Result_6FMR,[],1);
%     Allresult(9,:) = maxress(3:10);
%     
%     fprintf('- Finish FMR \n ');
%     
%     %% **************************** 7 GMC ****************************
%     [y_pred, U, S0, S0_initial, F, evs] = GMC(X, k);
%     Allresult(10,:) = Clustering8Measure(Y, y_pred);
%     fprintf('- Finish GMC \n ');
%     
%     %% **************************** 8 LMVSC ****************************
%     % Parameter 1: number of anchors (tunable)
%     numanchor=[k 50 100];
%     % Parameter 2: alpha (tunable)
%     alpha=[.01 1 100];
% 
%     Result_8LMVSC = zeros(length(numanchor)*length(alpha), 2+8);
%     idx=1;
%     for j=1:length(numanchor)
%     % Perform K-Means on each view
%         parfor i=1:viewN
%             rand('twister',5489);
%             [~, H_lm{i}] = litekmeans(data{i}, numanchor(j),'MaxIter', 100,'Replicates',10);
%         end
%         for i=1:length(alpha)
%             % Core part of this code (LMVSC)
%             [F,ids] = lmv(data,Y,H_lm,alpha(i)); % X:N*d
%             % Performance evaluation of clustering result
%             ress_lmvsc=Clustering8Measure(ids,Y);
%             Result_8LMVSC(idx,:) = [alpha(i) numanchor(j) ress_lmvsc];
%             idx=idx+1;
%         end
%     end
%     maxress = max(Result_8LMVSC,[],1);
%     Allresult(11,:) = maxress(3:10);
% 
%     fprintf('- Finish LMVSC \n ');
% 
%     %% **************************** 9 PMSC  ****************************
%     para1= 50 : 50 : 250;
%     para2= 10.^[-2 -1 0];
%     para3= 10.^( -5 : 1: -2 );
% 
%     Result_9PMSC = zeros(length(para1)*length(para2)*length(para3),3+8);
% 
%     idx=1;
%     Bp=cell(size(data));
%     for i=1:length(para1)
%         for ii=1:size(data,1)
%             Bp{ii}= inv(data{ii}*data{ii}'+para1(i)*eye(size(data{1},1)));
%         end
%         for j=1:length(para2)
%             for q=1:length(para3)
%                 [result,YY,Fv,S]=secondnorm22(data,Y,Bp,para1(i),para2(j),para3(q));% Input : N*d
%                 Result_9PMSC(idx,:) = [para1(i),para2(j),para3(q),result];
%                 idx = idx+1;
%             end
%         end
%     end
%     maxress = max(Result_9PMSC,[],1);
%     Allresult(12,:) = maxress(4:11);
%     fprintf('- Finish PMSC \n ');
%  
%     %% **************************** 10 MVSC ****************************
%     Wall=zeros(N,N);
%     for iv = 1:viewN
%         Wall = Wall + WV{iv};
%     end
%     D = diag(sum(Wall,2));
%     D2 = diag(1./sqrt(diag(D)+eps));
%     L   = D2 - Wall;
%     Lnor = D2*L*D2;
%     Lnor = (Lnor + Lnor')/2;
% 
%     options.tol = 1e-8;
%     options.maxit = 30000;
%     [Finit, ~] = eigs(Lnor, k, 'sa', options);
% 
%     maxIters = 100;
%     lambda1 = 10.^(-5:2:5);
%     lambda2 = 10.^(-5:2:5);
%     alpha = 10.^(1:1:2);
%     Result_10MVSC=zeros(length(lambda1)*length(lambda2)*length(alpha),3+8);
% 
%     idx=1;
%     for i = 1:length(lambda1)
%         for j = 1:length(lambda2)
%             for q = 1:length(alpha)
%                 [Z, F, E] = mvcsolver( data, k, Finit, lambda1(i), lambda2(j), alpha(q), maxIters );% Input :N*d
%                 Ypred = kmeans(F,k);
%                 result = Clustering8Measure(Y, Ypred);
%                 Result_10MVSC(idx,:) = [lambda1(i), lambda2(j), alpha(q), result];
%                 idx=idx+1;
%             end
%         end
%     end
%     maxress = max(Result_10MVSC,[],1);
%     Allresult(13,:) = maxress(4:11);
% 
%     fprintf('- Finish MVSC \n '); 
% 
%     save(['.\Result0930_' , dataName , '.mat'], 'Allresult',...
%                 'Result_11COMVSC', 'Result_1Co_pw', 'Result_1Co_ct', 'Result_2MLRSSC', 'Result_3LMSC',...
%                 'Result_4RMKM', 'Result_5mPAC', 'Result_6FMR', 'Result_8LMVSC', 'Result_9PMSC', 'Result_10MVSC');
end % Finish datasets
