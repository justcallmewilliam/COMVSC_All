clc;clear;
dataname = {'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};
dataname = {'3sources', 'ORL_mtv', 'proteinFold','WebKB_cor2views',...
'WebKB_Wisconsin2views', 'yaleA_3view','WebKB', 'WebKB_2views',...
'bbcsport_seg14of4', 'Handwritten_numerals',...
'MSRCV1','WikipediaArticles','Caltech101-7','Caltech101-20'};

for di = 1:10

    f = load(['../0datasets/smalldata/',dataname{di}]);
    disp(dataname{di});
    dlmwrite('result-0929.txt',[dataname{di}],'-append','delimiter','\t','newline','pc');

% f=load('bbc_seg14of4.mat');%=load('C:\Users\User\Desktop\research\multiviewlearning\shiguoxin\bbc_seg14of4.mat');
data=f.X;
label=f.Y;
% addpath('C:\Users\User\Desktop\research\kernelclusteringexp')
para1= [50 : 50 : 250];
para2= 10.^[-2 -1 0];
para3= 10.^[ -5 : 1: -2 ];

para1=250;
para2=0.1;
para3=0.001;

% for i=1:size(data,1)
% dist = max(max(data{i})) - min(min(data{i}));
% m01 = (data{i} - min(min(data{i})))/dist;
% data{i} = 2 * m01 - 1;
% end
B=cell(size(data));
for i=1:length(para1)
   for ii=1:size(data,1)
           tic;
           B{ii}= inv(data{ii}*data{ii}'+para1(i)*eye(size(data{1},1)));
   end
    for j=1:length(para2)
        for k=1:length(para3)
            fprintf('params%12.6f%12.6f%12.6f\n',para1(i),para2(j),para3(k));
            [result,YY,Fv,S]=secondnorm22(data,label,B,para1(i),para2(j),para3(k));
            t=toc;
            dlmwrite('result-0929.txt',[para1(i) para2(j) para3(k)...
                result(1,1) result(2,1) result(3,1) result(4,1) t],'-append','delimiter','\t','newline','pc');
        end
    end
end
end
