function [ress]=clustering(S, cls_num, gt)

[C] = SpectralClustering(S,cls_num);
% [A nmi avgent] = compute_nmi(gt,C);
% ACC = Accuracy(C,double(gt));
% [f,p,r] = compute_f(gt,C);
% [AR,RI,MI,HI]=RandIndex(gt,C);

 ress = Clustering8Measure(gt,C);
end
