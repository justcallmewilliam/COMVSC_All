clc;clear;
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};

for di = 1:10
    path = '/home/zpp/NewDisk/2020/2020TKDE/major_revision/AllMethod/0datasets/smalldata/';
    datapath = [path,dataname{di},'.mat'];
    disp(dataname{di});
    f = load(datapath);
    X=f.X;
    Y=f.Y;
    clear sigma;
    num_views = length(X);
    numClust = length(unique(Y));
    for i = 1:num_views
        sigma(i) = optSigma(X{i});
    end
    lambda = [0.01 0.02 0.03 0.04 0.05]; 
    numiter = 50;

    %% single best view
    % fprintf('Running with view1\n');
    % [E F P R nmi avgent] = baseline_spectral(X1,numClust,sigma1,truth);
    % fprintf('Running with view2\n');
    % [E F P R nmi avgent] = baseline_spectral(X2,numClust,sigma2,truth);
    % fprintf('Running with view3\n');
    % [E F P R nmi avgent] = baseline_spectral(X3,numClust,sigma3,truth);
%     for i = 1:num_views
%         [E F P R nmi avgent AR ress] = baseline_spectral(X{i},numClust,sigma(i),Y);
%         disp([num2str(i),'-view: ',num2str(ress)]);
%     end
    %% two views
    % % feature concat
    % fprintf('Running with the feature concatenation of two views\n');
    % [E F P R nmi avgent] = baseline_spectral([X1 X2],numClust,optSigma([X1 X2]),truth);
    % %de sa mindis
    % fprintf('Running with the De Sa MinDis\n');
    % [F P R nmi avgent AR] = spectral_mindis(X1, X2, numClust,sigma1,sigma2, truth);
    % fprintf('Our approach with 2 views (pairwise)\n');
    % [U U1 F P R nmi avgent AR] = spectral_pairwise(X1,X2,numClust,sigma1,sigma2,lambda,truth,numiter);
    % fprintf('Our approach with 2 views (centroid)\n');
    % lambda1 = 0.5; lambda2 = 0.5;
    % [U U1 F P R nmi avgent AR] = spectral_centroid(X1,X2,numClust,sigma1,sigma2,lambda1,lambda2,truth,numiter);
    % 
    % 
    % %% three views
    % fprintf('Running with the feature concatenation of three views\n');
    % [E F P R nmi avgent] = baseline_spectral([X1 X2 X3],numClust,optSigma([X1 X2 X3]),truth);
    
    %% multiview spectral (pairwise): more than 2 views
    lambda = [0.04 0.05];
    for li =1:length(lambda)
        fprintf('pairwise:Multiview spectral with 2 views \n ');
        disp(['lambda:' num2str(lambda(li))]);
        for it = 1:1
        tic;
        [F P R nmi avgent AR ress] = spectral_pairwise_multview(X,num_views,numClust,sigma,lambda(li),Y,numiter);
        t1(it)=toc;
        disp(ress);
        end
        disp(t1);
        tmean1 = mean(t1);
        disp(tmean1);
%         filename=['./result0912/',datasetname{di},'_pairwise0113.mat'];
%         save(filename,'F','P','R','nmi','avgent','AR','ress');
         savetxt = ('./major_result/coreg_pairwise0930.txt');
         % Fscore nmi ACC Purity
         dlmwrite(savetxt,[char(dataname{di}) ': ' num2str(lambda(li)),':  ',num2str(ress), '-',num2str(t1)],'-append','delimiter','\t','newline','pc');
    end
    %% multiview spectral (centroid): more than 2 views
    fprintf('centroid:Multiview spectral with 2 views \n ');
    lambda2 = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    for it = 1:1
        tic;
        [F P R nmi avgent AR ress] = spectral_centroid_multiview(X,num_views,numClust,sigma,lambda2,Y,numiter);
        t2(it)=toc;
        disp(ress);
    end
    disp(t2);
    tmean2 = mean(t2);
    disp(tmean2);
%     filename=['../result0912/',datasetname{di},'_centroid.mat'];
%     save(filename,'F','P','R','nmi','avgent','AR','ress');
     savetxt = ('./major_result/coreg_centroid0930.txt');
     dlmwrite(savetxt,[char(dataname{di}) ': ' num2str(ress), '-', num2str(t2)],'-append','delimiter','\t','newline','pc');
    
end