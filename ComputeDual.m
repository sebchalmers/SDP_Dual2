function [ D, res ] = ComputeDual( X, lambda, Q, C, P, Cset, WeightNuclear, tau )

    %Some dimensions
    Nagent = length(X);

  
    res = 0;D = 0;
    for agent = 1:Nagent

        %Dual value
        DC = lambdaC( lambda, C{agent}, P{agent} );
        D = D + trace((Q{agent}  + WeightNuclear*eye(size(Q{agent}))  + DC)*X{agent}) - tau*log(det(X{agent}));

        %Residual (to be checked)
        for const = 1:length(Cset) %walk through the C-sets
            res(const,1) = 0; %res for constraint 'const'
            for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
                participating_agent = Cset{const}(index);
                res(const,1) = res(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
                %display(['Agent = ',num2str(participating_agent),' / Const = ',num2str(const)]);
            end
        end
    end

end

