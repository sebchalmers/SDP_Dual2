function [ W, G, IG ] = computeWandG( X, Z )
%Compute the matrices W and G of the NT method
    
    L = chol(X,'lower'); R = chol(Z,'lower');
    [U,D,V] = svd(R.'*L);
    G = L*V*diag(1./sqrt(diag(D)));
    W = 0;%G*G.';
    IG = sqrt(D)*V.'*inv(L);
    %check = norm(X - W*Z*W)
end

% check = norm(X - W*Z*W)
