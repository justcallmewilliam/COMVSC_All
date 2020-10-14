function [p] = calculate(F)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
[sample,classnumber]= size(F);
p =zeros(sample,sample);
for i=1:sample
    for j=i+1:sample
        p(i,j) = norm(F(i,:)- F(j,:),2);
    end
end
p = p + p';

end

