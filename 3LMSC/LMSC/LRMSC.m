 function [ress,H,cov_val] = LRMSC(X,gt,cls_num,lambda,K)
%% Initialize P,H,Z,E_h,E_r,J,F and Y_1,Y_2,Y_3,Y_4
% X is d*N
V = length(X);
N = size(X{1},2); % number of data points
%K = 100; % dimension of H

for v=1:V
    alpha{v} = 1/V;
end

for i=1:V
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
    %X{i} = X{i}(:,1:100);
end

for i=1:V
    D{i} = size(X{i},1); % dimension of each view
end
SD = 0;
M = [];
for i=1:V
    SD = SD + D{i};
    M = [M;X{i}];
end
%MX = M;
P = zeros(SD,K);
H = zeros(K,N); Er = zeros(K,N); 
Z = zeros(N,N);J = zeros(N,N);
Eh = zeros(SD,N);
Y1 = zeros(SD,N);Y2 = zeros(K,N);Y3 = zeros(N,N);
% [COEFF,SCORE] = pca(M');
% T = SCORE(:,1:K);
% H = T';
% H= mypca(M',K);
% H=H';
H = rand(K,N);
%% updating variables...
IsConverge = 0;
mu = 1e-4;
%lambda = 100;
pho = 1.2;
max_mu = 1e6;
max_iter = 100;
iter = 1;
% a = sum(M,2)/size(M,2);
% b = sum(H,2)/size(H,2);
% M= a./repmat(sqrt(sum(a.*a,1)+eps),size(a,1),1);
% H =b./repmat(sqrt(sum(b.*b,1)+eps),size(b,1),1);

% M = M./repmat(sqrt(sum(M.^2,1)),size(M,1),1);
% H = H./repmat(sqrt(sum(H.^2,1)),size(H,1),1);
while (IsConverge == 0&&iter<max_iter+1)
    % update P
    P=updateP2(Y1,mu,M,H,Eh);
    % update H
    A = mu*(P'*P); B = mu*(Z*Z'-Z-Z'+eye(N))+eye(N)*1e-6;
    C = P'*Y1+Y2*Z'-Y2+mu*(P'*M+Er-P'*Eh-Er*Z');
    H = lyap(A,B,C);
    % update J
    %[J] = solve_l1l2((Z-Y3/mu)',lambda/mu);
    %J = J';
    J = softth((Z-Y3/mu)+eye(N)*1e-8,lambda/mu);
    % update multipliers
    % update Z
    Z = inv(H'*H+eye(N))*((J+H'*H-H'*Er)+(Y3+H'*Y2)/mu);
    % update E
    G = [M-P*H+Y1/mu; H-H*Z+Y2/mu];
    [E] = solve_l1l2(G,1/mu);
    Eh = E(1:SD,:); 
    Er = E(1+SD:SD+K,:);
   % updata multipliers
    Y1 = Y1+ mu*(M-P*H-Eh);
    Y2 = Y2+ mu*(H-H*Z-Er);
    Y3 = Y3+ mu*(J-Z);
    mu = min(pho*mu, max_mu);
    % convergence conditions
    thrsh = 0.000001;
    if(norm(M-P*H-Eh,inf)<thrsh && norm(H-H*Z-Er,inf)<thrsh && norm(J-Z,inf)<thrsh)
        IsConverge = 1;
    end
    cov_val(iter) = norm(H-H*Z-Er,inf);
    %norm(H-H*Z-Er,inf)
    %norm(J-Z,inf)
    
    iter = iter + 1;
end
[ress]=clustering(abs(Z)+abs(Z'), cls_num, gt);