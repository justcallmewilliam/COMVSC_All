function E = solveE(X, Z, lambda)

D = X - X*Z;

[row, col] = size(D);
for i=1:row
    for j=1:col
        if D(i, j)>lambda/2
            E(i, j) = D(i, j) - lambda/2;
        elseif D(i, j) < -lambda/2
            E(i, j) = D(i, j) + lambda/2;
        else
            E(i, j) = 0;
        end       
    end
end
