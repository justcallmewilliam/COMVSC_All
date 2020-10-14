%% LMSC (CVPR-17)
clc;clear;
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};

for di = 9:10
    path = '/home/zpp/NewDisk/2020/2020TKDE/major_revision/AllMethod/0datasets/smalldata/';
    datapath = [path,dataname{di},'.mat'];
    disp(dataname{di});
    load(datapath);
    dlmwrite('./LMSC_0929.txt',dataname{di},'-append','delimiter','\t','newline','pc');

    fprintf('Latent representation multiview subspace clustering\n');
    num_views = length(X);
    numClust = size(unique(Y),1);
    %Y=Y+1;% if Y has 0
    for i = 1:num_views
        X{i}=X{i}';
    end
    lambda = 10.^[-3 :2:3];
    
    lambda = 0.1;
    
    K = 100; 

    for i=1:length(lambda)
        disp(['lmd = ', num2str(lambda(i))]);
        for it = 1:1
            tic;
            [nmi(it), ACC(it), f(it), RI(it), H, cov_val] = LRMSC(X,Y,numClust,lambda(i),K);
            t(it) = toc;
            % disp([nmi(it), ACC(it), f(it), t(it)]);
        end
        nmimean = mean(nmi);
        accmean = mean(ACC);
        fmean = mean(f);
        RImean = mean(RI);
        tmean = mean(t);
        % disp(tmean);
%         savefile=['H:\\zp\\var\\LMSC\\caltech7\\','lmd=',num2str(lambda(i)),'.mat'];
%         save(savefile,'nmi','ACC','f','RI','H','cov_val');
        
        dlmwrite('./LMSC_0929.txt',[num2str(lambda(i)),' - ', num2str(nmimean),...
            num2str(accmean),num2str(fmean),'-',num2str(tmean)],'-append','delimiter','\t','newline','pc');
    end
end