% 
% Run MLRSSC on UCI digit dataset. Parameters are optimized over NMI
% measure.
%
%-------------------------------------------------------
addpath ./clustering/ ./performance/
clc;clear;
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};

for di = 3:10
    path = '/home/zpp/NewDisk/2020/2020TKDE/major_revision/AllMethod/0datasets/smalldata/';
    datapath = [path,dataname{di},'.mat'];
    disp(dataname{di});
    f = load(datapath);
    X      = f.X;
    truth = f.Y; 

    k = length(unique(truth));
    num_views = length(X);
    num_iter = 100;
    %% Linear kernel multi-view LRSSC

%    fprintf('\nPairwise multiview LRSSC\n');
%     opts.mu = 10^2;
%     lambda1 = 0.5;
%     lambda3 = 0.7;
%     opts.mu = 10;
%     lambda1 = [0.1 0.3 0.5 0.7 0.9];
%     lambda3 = [0.5 0.7, 0.9];
%     for i=1:length(opts.mu)
%         for j = 1:length(lambda1)
%             for k = 1:length(lambda3)
%                 disp([num2str(opts.mu(i)),'_',num2str(lambda1(j)),'_',num2str(lambda3(k))]);
%                 opts.lambda = [lambda1 (1-lambda1) lambda3];
%                 opts.noisy = false;
%                 for it = 1:1
%                     tic;
%                     A = pairwise_MLRSSC(X, opts); % joint affinity matrix
%                     [result] = spectral_clustering(A, k, truth);
%                     t1(it) = toc;
%                 end
%                 disp(result);
%                 disp(t1);
%                 tmean1 = mean(t1);
%                 disp(tmean1);

% %                 folder = ['../result/LRSSC/pairwise/',datasetname{di},'/'];
% %                 if exist(folder)==0   %该文件夹不存在，则直接创�?                    
% %                     mkdir(folder);
% %                 else
% %                     dir=1;
% %                 end
% %                 filename=[folder,num2str(opts.mu(i)),'_',num2str(lambda1(j)),'_',num2str(lambda3(k)),'_Pairwise','.mat'];
% %                 save(filename,'result');
% %                 savetxt = ['../result/LRSSC/pairwise/',datasetname{di},'_','Pairwise_0115.txt'];
% %                 dlmwrite(savetxt,[opts.mu(i) lambda1(j) lambda3(k)  result],'-append','delimiter','\t','newline','pc');
%             end
%         end
%     end


    fprintf('\nCentroid multiview LRSSC\n');
    opts.mu = [10, 10^2, 10^3, 10^4];
    opts.mu = 10;
    lambda1 = [0.1 0.3 0.5 0.7 0.9];
    lambda3 = [0.3 0.5 0.7 0.9];

    for i=1:length(opts.mu)
        for j = 1:length(lambda1)
            for k = 1:length(lambda3)
                opts.lambda = [lambda1 (1-lambda1) lambda3];
                disp([num2str(opts.mu(i)),'_',num2str(lambda1(j)),'_',num2str(lambda3(k))]);
                for it = 1:1
                    tic;
                    A = centroid_MLRSSC(X, opts); % joint affinity matrix
                    result = spectral_clustering(A, k, truth);
                    t(it) = toc;
                    %best
                end
                % disp(result);
                % disp(t);
                tmean = mean(t);
                % disp(tmean);
%                 folder = ['../result/LRSSC/centroid/',datasetname{di},'/'];
%                 if exist(folder)==0   %该文件夹不存在，则直接chuangjian                    
%                     mkdir(folder);
%                 else
%                     dir=1;
%                 end
%                 filename=[folder,num2str(opts.mu(i)),'_',num2str(lambda1(j)),'_',num2str(lambda3(k)),'_Centroid','.mat'];
%                 save(filename,'result');
                 savetxt = ['./major_result/',dataname{di},'_','Centroid_0929.txt'];
                 dlmwrite(savetxt,[opts.mu(i) lambda1(j) lambda3(k)  result tmean],'-append','delimiter','\t','newline','pc');
            end
        end
    end
% 
%     disp(['===== Finish ', dataname{di},' =====']);
%     %% Gaussian kernel multi-view LRSSC
% 
%     opts.kernel = 'Gaussian';
%     opts.err_thr = 10^(-5);
%     for v=1:num_views
%        sigma(v) = opt_sigma(X{v});
%     end
%     opts.sigma = [5*sigma(1) 0.5*sigma(2) 0.5*sigma(3)]; 
% 
%     fprintf('\nKernel pairwise multiview LRSSC\n');
%     opts.mu = 10^4;
%     lambda1 = 0.7;
%     lambda3 = 0.7;
%     opts.lambda = [lambda1 (1-lambda1) lambda3];
% 
%     A = pairwise_MLRSSC(X, opts); % joint affinity matrix
%     [best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
%     best
% 
%     fprintf('\nKernel centroid multiview LRSSC\n');
%     opts.mu = 10^4;
%     lambda1 = 0.7;
%     lambda3 = 0.5;
%     opts.lambda = [lambda1 (1-lambda1) lambda3];
% 
%     A = centroid_MLRSSC(X, opts); % joint affinity matrix
%     [best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
%     best

end
