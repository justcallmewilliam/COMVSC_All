clear;
dataname = {'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
for di = 9:10
    path = '/home/zpp/NewDisk/2020/2020TKDE/major_revision/AllMethod/0datasets/';
    datapath = [path,dataname{di},'.mat'];
    disp(dataname{di});
    dlmwrite('result0930.txt',dataname{di},'-append','delimiter','\t','newline','pc');
    
    f = load(datapath);
    data=f.X;
    label=f.Y;
    viewN = length(data);
    
    para1=[10 20 30];
    para2=[.001 .01 .1];
    para3=[.000001 .00001 .0001];

%     para1 = 30; para2= 0.01;    para3=0.000001; % MSRCV1
%     para1 = 20;para2= 0.001; para3=0.000001; % Wiki
%     para1 = 10;para2= 0.001; para3=0.00001; % Caltech-7
%     para1 = 10;para2= 0.1; para3=0.000001; % Caltech-20 

    for i=1:viewN
        dist = max(max(data{i})) - min(min(data{i}));
        m01 = (data{i} - min(min(data{i})))/dist;
        data{i} = 2 * m01 - 1;
    end
    [v1,v2]=size(data);
    FV=cell(v1,v2);
    SV=cell(v1,v2);
    LV=cell(v1,v2);
    for i=1:viewN
        Sv= constructW_PKN(data{i}', 5, 1);
        D = diag(sum(Sv));
        Lv = D-Sv;
        [Fv, ~, ev]=eig1(Lv,length(unique(label)),0);
        FV{i}=Fv;
        SV{i}=Sv;
        LV{i}=Lv;
    end

    B=cell(size(data));
    for i=1:length(para1)
        for ii=1:viewN
            B{ii}= inv(data{ii}*data{ii}'+para1(i)*eye(size(data{1},1)));
        end
        for j=1:length(para2)
            for k=1:length(para3)
                disp([num2str(para1(i)),'_',num2str(para2(j)),'_',num2str(para3(k))]);
                for it = 1:1
                    tic;
                    [res(it,:),FVret,obj,Wv]=UMVSC(data,label,B,para1(i),para2(j),para3(k),FV,SV,LV);
                    t(it) = toc;
                    % disp(res(it,:));
                    % disp(t(it));
                end
                meanres = mean(res,1);
                meant = mean(t,2);
                disp(['meanresult ', num2str(meanres)]);
                disp(['mean time  ', num2str(meant)]);
%                 folder = ['./result/',newdata{di},'/'];
%                 if exist(folder)==0   %该文件夹不存在，则直接创�?
%                     mkdir(folder);
%                 else
%                     dir=1;
%                 end
%                 filename=[folder,num2str(para1(i)),'_',num2str(para2(j)),'_',num2str(para3(k)),'.mat'];
%                 save(filename,'res','FVret','obj','Wv');
                 
                 dlmwrite('result0930.txt',[para1(i) para2(j) para3(k) meanres(1) meanres(2)...
                     meanres(3) meanres(4) meanres(5) meanres(6) meanres(7) meanres(8) meant],...
                     '-append','delimiter','\t','newline','pc');

            end
        end
    end
end
