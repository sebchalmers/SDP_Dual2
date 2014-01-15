function [ alpha ] = alphaPos( X, dX )
%Find step-size to maintain positivity

    alpha = -0.98/min(eig(X\dX));
    
    if alpha < 0
        alpha = 1;
    else
        alpha = min([alpha,1]);
    end
    
end

