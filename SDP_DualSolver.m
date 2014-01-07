clear all
close all
clc

% run /Users/sebastien/Desktop/cvx/cvx_setup

% Test on 2 agents

%First or Second order Dual update
%Method    = 'GA';     % Gradient ascent
Method    = 'FGM';    % Fast gradient 
%Method    = 'Second'; % 2nd-order dual update

Lipschitz = 10; % For 'GA' and 'FGM'

display_subproblems = 0;

Nvar   = 8;
Nconst = 2;

% % Four agents
% Nagent = 4;
% 
% %Define coupling constraints (C-set)
% Cset = {[1 2],
%         [2 3],
%         [3 4]};

% % Three agents
Nagent = 3;

%Define coupling constraints (C-set)
Cset = {[1 2],
        [2 3]};

% res = [trace(C{1}(:,:,1)*X{1} + C{2}(:,:,1)*X{2});
%        trace(C{2}(:,:,2)*X{1} + C{3}(:,:,2)*X{2})]

% %% Two agents
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
%         while (min(real(eig(C{agent}(:,:,P{agent}(const))))) < 0)
%             C{agent}(:,:,P{agent}(const)) = C{agent}(:,:,P{agent}(const)) + (1-min(real(eig(C{agent}(:,:,P{agent}(const))))))*eye(Nvar);
%         end
    end


    for k = 1:Nconst
        A{agent}(:,:,k) = random('norm',0,1,Nvar,Nvar);
        A{agent}(:,:,k) = A{agent}(:,:,k) + A{agent}(:,:,k).';        
        a{agent}(k,1)   = random('norm',0,1);
    end

end

lambda    = zeros(length(Cset),1);
xFGM      = lambda; %for Fast Gradient Method
xFGM_prev = xFGM;

% Initiate subproblems
% for agent = 1:Nagent
%     cvx_begin
% 
%         cvx_precision( 1e-1 );
% 
%         variable Xk(Nvar,Nvar);% symmetric;
% 
%         dual variables y{1+Nconst}  
% 
%         %Dualize constraints
%         DC_agent = lambdaC( lambda, C{agent}, P{agent} );
%         
%         minimize( trace((Q{agent} + WeightNuclear*eye(Nvar) + DC_agent)*Xk) );
%         subject to
%             Xk == semidefinite(Nvar) : y{1};
%             for k = 1:Nconst
%                 trace(A{agent}(:,:,k)*Xk) == a{agent}(k) : y{k+1} ;
%             end
% 
%         cvx_problem
%     cvx_end
% 
%     Z0{agent} = y{1};
%     for k = 1:Nconst
%         mu0{agent}(k) = y{k+1};
%     end
%     X0{agent} = Xk;
% end
% X = X0;Z = Z0;mu = mu0;

tolIP = 1e-6;tauEnd = -6; %(order)
tau_table = logspace(0,tauEnd,20);

%Naive guess
for agent = 1:Nagent
    X{agent} = eye(Nvar); Z{agent} = eye(Nvar); 
    mu{agent} = zeros(Nconst,1);
end
    

for tau_index = 1:length(tau_table)
    tau     = tau_table(tau_index);
    tolDual = tau;
    
    iter = 1;


    res = 1e6;DualIncrease = 1;
    while (norm(res) > tolDual) && (DualIncrease == 1)
        D = 0;res = 0;
        for agent = 1:Nagent
            [ X{agent}, Z{agent}, mu{agent}, X_sens{agent} ] = NTSolve(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau, tolIP, X{agent}, Z{agent}, mu{agent}, P{agent}, display_subproblems);

            %Dual value
            DC = lambdaC( lambda, C{agent}, P{agent} );
            D = D + trace((Q{agent}  + WeightNuclear*eye(Nvar)  + DC)*X{agent}) - tau*log(det(X{agent}));

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
        
        %Compute dual Hessian
        DH = zeros(length(Cset),length(Cset));
        for m = 1:length(Cset)
            for l = 1:length(Cset)
                for i = 1:length(Cset{l})
                    k = Cset{l}(i);
                    if max(m == P{k})                      
                        DH(m,l) = DH(m,l) + trace(C{k}(:,:,l)*X_sens{k}(:,:,m))   ;                     
                    end
                end
            end
        end
        
        Lipschitz = max(abs(eig(DH)));

        %Store
        lambda_store{tau_index}(iter,:) = lambda;
        D_store{tau_index}(iter)        = D;
        res_store{tau_index}(iter)      = norm(res);
        DH_store{tau_index}(iter)       = Lipschitz;
        
        %%%%%%%%%%%%%%%%
        %  Dual Update %
        %%%%%%%%%%%%%%%%
        
        switch Method
            case 'GA'
                lambda      = lambda + res/Lipschitz;
            case 'FGM'
                xFGM      = lambda + res/Lipschitz;
                lambda    = xFGM + (iter-1)/(iter+2)*(xFGM - xFGM_prev);
                xFGM_prev = xFGM;
            case 'Second'     
                lambda = lambda - DH\res;
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
            ylim([tolDual max(10*tolDual,res_store{tau_index}(1))])
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
            ylim([tolDual max(10*tolDual,res_store{tau_index}(1))])
            title('Residual')
            grid on
        end
        
        if (length(lambda) > 2)
            figure(2);clf
            semilogy(res_store{tau_index});hold on
            title('Residual')
            grid on
            ylim([tolDual max(10*tolDual,res_store{tau_index}(1))])
        end
        
    end
    iter_store(tau_index) = iter-1;

    Central_path(tau_index,:) = [lambda_store{tau_index}(end,:) D_store{tau_index}(end)];

    
    if (length(lambda) == 1)
        figure(3);
        plot(lambda_store{tau_index},D_store{tau_index},'linestyle',':','marker','.');hold on
        plot(lambda_store{tau_index}(end),D_store{tau_index}(end),'linestyle',':','marker','.','color','r');
        plot(Central_path(:,1),Central_path(:,2))
        title('Central path & with Dual values over local updates')
        xlabel('Dual variable');ylabel('Dual Value')
        grid on
    end
    
    if (length(lambda) == 2)
        figure(3);
        plot3(lambda_store{tau_index}(:,1),lambda_store{tau_index}(:,2),D_store{tau_index},'linestyle',':','marker','.');hold on
        plot3(Central_path(:,1),Central_path(:,2),Central_path(:,3))
        title('Central path & with Dual values over local updates')
        xlabel('Dual variable');ylabel('Dual Value')
        grid on
    end
    

    
    figure(4)
    barh(log10(tau_table(1:tau_index)),iter_store)
    ylabel('log10(Barrier)');xlabel('#Dual iter')
    grid on    
    title('# of iteration at each barrier value')
    %pause
end

if (length(lambda) == 2)
    Check_residual = [trace(C{1}(:,:,1)*X{1} + C{2}(:,:,1)*X{2});
                      trace(C{2}(:,:,2)*X{2} + C{3}(:,:,2)*X{3})]
end

if (length(lambda) == 1)
    Check_residual = [trace(C{1}(:,:,1)*X{1} + C{2}(:,:,1)*X{2})]
end