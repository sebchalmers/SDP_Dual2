function [ W, G, IG ] = computeWandG( X, Z )
%Compute the matrices W and G of the NT method
    
    try
        L = chol(X,'lower'); 
    catch
        display('Cholesky factorization of X failed')
        display('Eig of X')
        min(eig(X))
        keyboard
    end
    
    try
        R = chol(Z,'lower');
    catch
        display('Cholesky factorization of Z failed')
        display('Eig of Z')
        min(eig(Z))
        keyboard
    end
    
    [U,D,V] = svd(R.'*L);
    G = L*V*diag(1./sqrt(diag(D)));
    W = 0;%G*G.';
    IG = sqrt(D)*V.'*inv(L);
    %check = norm(X - W*Z*W)
end

% check = norm(X - W*Z*W)
