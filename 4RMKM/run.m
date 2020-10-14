clc;clear;
% dataname = {'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
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
    N = length(Y); %% Number of samples
    K = length(unique(Y)); %% Clusters of samples
    ViewN = size(X,2); 
    inXCell = cell(ViewN,1);
    for iv =1:ViewN
        inXCell{iv}=X{iv}';  
    end
    for vIndex=1:ViewN
        Temp=[];Temp=inXCell{vIndex};
        inXCell{vIndex}=minmaxNormalization(Temp); %Normalize into [-1,1]
    end
    % init G0
    Ik = eye(K);
    randorder = randperm(size(Ik,1));
    numceil = ceil(N/K);
    largeG = repmat(Ik(randorder,:),numceil,1);
    inG0= largeG(1:N,:);% N*k

    %%
    gma = linspace(1.1,10,20);
    
%    gma = 4.3789; % MSRCV1
%     gma = 8.5947; % Wiki
%     gma = 9.5316; % Caltech-7
     gma = 6.7211; % Caltech-20
    
    for i = 1:length(gma)
        disp(['gma = ', num2str(gma(i))]);
        for it = 1:1
            tic;
            [ outG0, outFCell, outAlpha, outAlpha_r, outObj, outNumIter ] = weighted_robust_multi_kmeans( inXCell,gma(i) ,inG0 );
            [~,res_label]=max(outG0,[],2);
            result(it,:) = Clustering8Measure(res_label,Y);
            t(it) = toc;
            disp(result(it,:));
            disp(t(it));
        end
        resultmean = mean(result,1);
        tmean = mean(t);
%        disp(tmean);
%         save(['H:\\zp\\var\\RMKM\\MSRCV1\\','GMA=',num2str(gma(i)),'.mat']);
        rmresult = [gma(i), resultmean, tmean];
         dlmwrite('RMKM_0929.txt', rmresult, '-append','delimiter','\t','newline','pc');
    end
end

