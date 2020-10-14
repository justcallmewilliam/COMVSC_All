function Z = solveZ(X, F, Z, E, lamda, alpha)

[row, col] = size(F);


for i=1:row
    for j=1:row
        deltaF = F(i,:) - F(j,:);
        P(i, j) = deltaF*deltaF';
    end
end

numOfSamples = size(X, 2);
%X = matrixNormalize(X);
X = [X; alpha*ones(1, numOfSamples)]; 

E = [E; zeros(1, numOfSamples)];


for i=1:row
    X_1 = X - (X*Z - X(:, i)*Z(i, :)) - E;
    v = X_1'*X(:, i)/(X(:, i)'*X(:, i)+eps);
    
    for j=1:row
        temp = lamda * P(j, i) / 4;
        if v(j) > temp
            Z(i, j) = v(j) - temp;
        elseif v(j) < -temp
            Z(i, j) = v(j) + temp;
        else
            Z(i, j) = 0;
        end
    end
    
    Z(i, i) = 0;
end


