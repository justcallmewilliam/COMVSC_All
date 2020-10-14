function [ Xnor ] = minmaxNormalization(X)
%[-1 , 1]
% Xnormalization = (x - xmin)/(xmax- xmin)
Xnor = (X - min(min(X)))./(max(max(X)) - min(min(X)));
Xnor = 2*Xnor-1;
end

