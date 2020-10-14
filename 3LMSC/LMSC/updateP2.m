function P=updateP2(Y1,mu,M,H,Eh)
G = H';
Q = (1/mu*Y1+M-Eh)';
W = G'*Q;
[U,S,V] = svd (W,'econ'); 

PT = U*V';
P = PT';
end
