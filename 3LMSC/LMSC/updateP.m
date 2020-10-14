function P=updateP(Y1,K,mu,P,M,H,Eh)
Q = zeros(size(P));
S = zeros(size(P));
beta = 1e-4;
rho = 1.2;
max_beta = 1e6;
epson = 1e-4;
max_iter = 50;
IsConverge = 0;
iter = 1;
while (IsConverge == 0&&iter<max_iter)
   % old_cst = norm(M-P*H-Eh,'fro');
    P = (Y1*H'+mu*(M*H'-Eh*H')+beta*(Q-S))*inv(mu*(H*H')+beta*eye(K));
    %P = (M*H')*(inv(H*H'+10e1*eye(K)));
    Q = S + P; 
    %Q = Q./repmat(sqrt(sum(Q.^2,1)),size(Q,1),1)+epson;
    Q = Q./repmat(max(1,sqrt(sum(Q.*Q,1)+eps)),size(Q,1),1);
    %Q = T';
    S = S + (P - Q);
    beta = min(rho*beta, max_beta);
    
   % new_cst = norm(M-P*H-Eh,'fro');
    if(norm(P-Q,inf)<epson)
        IsConverge = 1;
        a = 111
    end
    %if iter == 190
        %norm(P-Q,inf)
    %end
    iter = iter+1;
end

end
