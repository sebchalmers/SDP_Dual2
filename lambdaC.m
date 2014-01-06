function [ DC ] = lambdaC( lambda, C, P )
%Build \sum lambda_i C_i

    Nvar = size(C,1);
    
    %Build dualization
    DC = zeros(Nvar,Nvar);
    for i = 1:length(P)
        const = P(i);
        DC = DC + lambda(const)*C(:,:,const);
    end
end

