function [ GskronK ] = skron( G, K )
%Compute the skron product

    n = size(G,1);
    
    %Build U
    U = zeros(n*(n+1)/2,n^2);
    
    IndicesLines = [];
    for j = 1:n
        for i = j:n
            IndicesLines = [IndicesLines;[i j]];
        end
    end
    %IndicesLines
    
    IndicesCols = [];
    for l = 1:n
        for k = 1:n
            IndicesCols = [IndicesCols;[k l]];
        end
    end
    %IndicesCols
    
    for line = 1:size(U,1)
        for col = 1:size(U,2)
            i = IndicesLines(line,1);
            j = IndicesLines(line,2);
            k = IndicesCols(col,1);
            l = IndicesCols(col,2);
            
            if (i==j) && (j==k) && (k==l)
                U(line,col) = 1;
            end
            if ((i==k) && (j==l) && not(k==j)) || ((i==l) && (j==k) && not(l==j))
                U(line,col) = 1/sqrt(2);
            end
        end
    end


    GskronK = 0.5*U*(kron(G,K) + kron(K,G))*U.';

end



% %check skron
% skron(X,Z)*svec(A) - 0.5*svec(Z*A*X.' + X*A*Z.')
% skron(X,Z) - skron(Z,X)
% skron(X,Z).' - skron(X.',Z.')
% skron(X,eye(Nvar)).' - skron(X,eye(Nvar))
% norm(inv(skron(X,X)) - skron(inv(X),inv(X)))
