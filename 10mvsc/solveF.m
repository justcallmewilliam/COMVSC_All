function F = solveF(Z, numOfClasses)

numOfViews = length(Z);
numOfSamples = size(Z{1}, 1);
M = zeros(numOfSamples, numOfSamples);

for i=1:numOfViews
    W = (abs(Z{i})+abs(Z{i}'))/2;
  

    DN = diag( 1./sqrt(sum(W, 2)+eps) );
    LapN = speye(numOfSamples) - DN * W * DN;
    M = M + LapN;
    
end
M = M./numOfViews;

[V,D] = eig(M);
[D_sort, ind] = sort(diag(D));
ind2 = find(D_sort>1e-6);
F = V(:, ind2(1:numOfClasses));



