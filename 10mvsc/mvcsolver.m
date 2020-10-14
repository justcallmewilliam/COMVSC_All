function [Z, F, E] = mvcsolver(X, C, Finit, lambda1, lambda2, alpha, maxIters)

%: min \sum_{v} ||X_v - X_vZ_v - E_v||^2 + \lambda_1 \sum_{v}Tr(F^TL_vF) 
%                                                             + \lambda_2 \sum_{v} ||E_v||_1

% X: 1 by n cell, each cell corresponds to each view
% C: number of classes
% alpha: defined in the paper, tune it according to the performance


numOfViews = length(X);

for iv = 1:numOfViews
    X{iv} = X{iv}';
end
numOfSamples = size(X{1}, 2);


for i=1:numOfViews
    Z{i} = zeros(numOfSamples, numOfSamples);
    E{i} = zeros(size(X{i}));
end

    
F = Finit;


for it=1:maxIters
   
    %-------update Z---------
    for j=1:numOfViews
        Z{j} = solveZ(X{j}, F, Z{j}, E{j}, lambda1, alpha);
    end
    
    
    %------update E---------
    for j=1:numOfViews
        E{j} = solveE(X{j}, Z{j}, lambda2);
    end
    
 
    
    %------update F--------
    F = solveF(Z, C);
    
   
    %------print OBJ-------
    % fprintf('iters = %d\n', it);
    obj = 0.0;
    for j=1:numOfViews
        W = (abs(Z{j})+abs(Z{j}'))/2;
        L = diag(sum(W, 2)) - W;
        obj = obj + norm((X{j} - X{j}*Z{j}))^2 + lambda1*sum(diag(F'*L*F)) + lambda2*norm(E{j}, 1);
    end
    % fprintf('obj = %f\n', obj);
    
end





