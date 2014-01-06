clear all
close all
clc

% run /Users/sebastien/Desktop/cvx/cvx_setup

% Test on 2 agents

%First or Second order Dual update
Order = 'First';
Lipschitz = 10;

Nvar   = 3;
Nconst = 1;

% % Three agents
Nagent = 3;

%Define coupling constraints (C-set)
Cset = {[1 2],
        [2 3]};

% res = [trace(C{1}(:,:,1)*X{1} + C{2}(:,:,1)*X{2});
%        trace(C{2}(:,:,2)*X{1} + C{3}(:,:,2)*X{2})]

% % Two agents
% Nagent = 2;
% 
% %Define coupling constraints (C-set)
% Cset = {[1 2]};

%% res = trace(C{1}*X{1} + C{2}*X{2})

%Participating set
for k = 1:Nagent
    P{k} = []; %create set of empty sets
end
for const = 1:length(Cset)
    for k = 1:length(Cset{const})
        P{Cset{const}(k)} = [P{Cset{const}(k)} const];
    end
end

WeightNuclear = 10;

for agent = 1:Nagent
    Q{agent} = random('norm',0,1,Nvar,Nvar);
    Q{agent} = Q{agent} + Q{agent}.';
    while (min(real(eig(Q{agent}))) < 0)
        Q{agent} = Q{agent} + (1-min(real(eig(Q{agent}))))*eye(Nvar);
    end
    

    % C{agent}(:,:,index of coupling constraint)
    for const = 1:length(P{agent})
        C{agent}(:,:,P{agent}(const)) = random('norm',0,1,Nvar,Nvar);
        C{agent}(:,:,P{agent}(const)) = C{agent}(:,:,P{agent}(const)) + C{agent}(:,:,P{agent}(const)).';
    end


    for k = 1:Nconst
        A{agent}(:,:,k) = random('norm',0,1,Nvar,Nvar);
        A{agent}(:,:,k) = A{agent}(:,:,k) + A{agent}(:,:,k).';        
        a{agent}(k)     = random('norm',0,1);
    end

end

lambda = zeros(length(Cset),1);
% for agent = 1:Nagent
%     DC(:,:,agent) = lambdaC( lambda, C{agent}, P{agent} );
% end

for agent = 1:Nagent
    cvx_begin

        cvx_precision( 1e-1 );

        variable Xk(Nvar,Nvar);% symmetric;

        dual variables y{1+Nconst}  

        %Dualize constraints
        DC_agent = lambdaC( lambda, C{agent}, P{agent} );
        
        minimize( trace((Q{agent} + WeightNuclear*eye(Nvar) + DC_agent)*Xk) );
        subject to
            Xk == semidefinite(Nvar) : y{1};
            for k = 1:Nconst
                trace(A{agent}(:,:,k)*Xk) == a{agent}(k) : y{k+1} ;
            end

        cvx_problem
    cvx_end

    Z0{agent} = y{1};
    for k = 1:Nconst
        mu0{agent}(k) = y{k+1};
    end
    X0{agent} = Xk;
end

tol = 1e-6;
tau_table = logspace(0,-6,10);

X = X0; Z = Z0; mu = mu0;

for tau_index = 1:length(tau_table)
    tau = tau_table(tau_index);

    iter = 1;


    res = 1;DualIncrease = 1;
    while (norm(res) > 1e-3) && (DualIncrease == 1)
        D = 0;res = 0;
        for agent = 1:Nagent
            [ X{agent}, Z{agent}, mu{agent}, X_sens{agent} ] = NTSolve(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau, tol, X{agent}, Z{agent}, mu{agent}, P{agent});

            %Dual value
            DC = lambdaC( lambda, C{agent}, P{agent} );
            D = D + trace((Q{agent}  + WeightNuclear*eye(Nvar)  + DC)*X{agent});

            %Residual (to be checked)
            clc
            for const = 1:length(Cset) %walk through the C-sets
                res(const,1) = 0; %res for constraint 'const'
                for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
                    participating_agent = Cset{const}(index);
                    res(const,1) = res(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
                    %display(['Agent = ',num2str(participating_agent),' / Const = ',num2str(const)]);
                end
            end
            
%           %Check residual 2 agents
            res
            res2 = [trace(C{1}(:,:,1)*X{1}) + trace(C{2}(:,:,1)*X{2});
                    trace(C{2}(:,:,2)*X{2}) + trace(C{3}(:,:,2)*X{3})]
            %pause
            

        end

        %Compute dual Hessian
        DH = zeros(length(Cset),length(Cset));
        for m = 1:length(Cset)
            for l = 1:length(Cset)
                for i = 1:length(Cset{l})
                    k = Cset{l}(i);
                    DH(m,l) = DH(m,l) + trace(C{k}(:,:,l)*X_sens{agent}(:,:,m));
                end
            end
        end


        
        %Store
        lambda_store{tau_index}(iter,:) = lambda;
        D_store{tau_index}(iter)        = D;
        res_store{tau_index}(iter)      = norm(res);
        DH_store{tau_index}(iter)       = norm(DH);
        
        %%%%%%%%%%%%%%%%
        %  Dual Update %
        %%%%%%%%%%%%%%%%
        
        switch Order
            case 'First'
                lambda = lambda + res/Lipschitz;
            case 'Second'     
                lambda = lambda - DH\res;
        end

        if iter > 1
            %Check Dual Increase
            if (D_store{tau_index}(iter) < D_store{tau_index}(iter-1))
                DualIncrease = 0;
            end
        end
        
        %Next iterate
        iter = iter + 1;

        if (length(lambda) == 1)
            %Two agents
            figure(2);clf
            subplot(1,2,1)
            plot(lambda_store{tau_index},D_store{tau_index},'linestyle','none','marker','.');hold on
            title('Dual Value')
            subplot(1,2,2)
            semilogy(res_store{tau_index});hold on
            title('Residual')
            grid on
        end
        
        if (length(lambda) == 2)
            %Three agents
            figure(2);clf
            subplot(1,2,1)
            plot3(lambda_store{tau_index}(:,1),lambda_store{tau_index}(:,2),D_store{tau_index},'linestyle','none','marker','.');hold on
            title('Dual Value')
            grid on
            subplot(1,2,2)
            semilogy(res_store{tau_index});hold on
            title('Residual')
            grid on
        end
        
    end
    iter_store(tau_index) = iter;

    Central_path(tau_index,:) = [lambda_store{tau_index}(end,:) D_store{tau_index}(end)];

    
    if (length(lambda) == 1)
        figure(3);
        plot(lambda_store{tau_index},D_store{tau_index},'linestyle',':','marker','.');hold on
        plot(lambda_store{tau_index}(end),D_store{tau_index}(end),'linestyle','none','marker','.','color','r');
        plot(Central_path(:,1),Central_path(:,2))
        title('Central path & with Dual values over local updates')
        xlabel('Dual variable');ylabel('Dual Value')
        grid on
    end
    
    if (length(lambda) == 2)
        figure(3);
        plot3(lambda_store{tau_index}(:,1),lambda_store{tau_index}(:,2),D_store{tau_index},'linestyle','none','marker','.');hold on
        plot3(Central_path(:,1),Central_path(:,2),Central_path(:,3))
        title('Central path & with Dual values over local updates')
        xlabel('Dual variable');ylabel('Dual Value')
        grid on
    end
    
    figure(4)
    subplot(2,1,1)
    plot3([1:1:length(DH_store{tau_index})],-log10(tau_table(tau_index))*ones(length(DH_store{tau_index}),1),DH_store{tau_index});hold on;grid on
    title('Dual curvature at each barrier value')
    xlabel('Dual iteration');ylabel('-log10(Barrier)');zlabel('||Dual Hessian||')

    subplot(2,1,2)
    bar(-log10(tau_table(1:tau_index)),iter_store)
    xlabel('-log10(Barrier)');ylabel('#Dual iter')
    grid on    
    title('# of iteration at each barrier value')
    %pause
end

if (length(lambda) == 2)
    Check_residual = [trace(C{1}(:,:,1)*X{1} + C{2}(:,:,1)*X{2});
                      trace(C{2}(:,:,2)*X{2} + C{3}(:,:,2)*X{3})]
end
