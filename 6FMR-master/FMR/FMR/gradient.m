function [H,val,k]=gradient(fun,gfun,X0,epsilon,alpha) %x0��H(:,m)�ĳ�ʼֵ
% ����: �������½��������Լ������:  min f(x)
%����:  x0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�����ݶ�
%���:  x, val�ֱ��ǽ������ŵ������ֵ,  k�ǵ�������.
maxk=20000;   %����������
% % % rho=0.5;sigma=0.4;
k=0;%  epsilon=-1.32e+3;


while(k<maxk)
    g=feval(gfun,X0);  %�����ݶ�
    d=-g;    %������������
     val0=feval(fun,X0);
    
     k=k+1;
    if(val0<epsilon)
        break; 
    end
      X0=X0+alpha*d;  
      % X0=X0./repmat(max(1,sqrt(sum(X0.*X0,1)+eps)),size(X0,1),1);
      val0=feval(fun,X0);   
end

H=X0;
val=feval(fun,X0);
end