function [p] = calculate(F)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[sample,classnumber]= size(F);
p =zeros(sample,sample);
for i=1:sample
    for j=i+1:sample
        p(i,j) = norm(F(i,:)- F(j,:),2);
    end
end
p = p + p';

end

