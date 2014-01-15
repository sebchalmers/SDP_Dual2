function [ check ] = PosCheck( X, Z )
%Check the positivity of the primal-dual solution (for debugging purposes)
    
    if (min(eig(X)) < 0) || (min(eig(Z)) < 0) 
        check = 0;
        keyboard
    else
        check = 1;
    end


end

