clear all
close all
clc

List_of_Methods = {'Second_Predictor','Second_NoPredictor','FGM'};
Marker = {'o','*','x'};

label = '';

save = 1;
display_subproblems = 0;

Nvar   = 8;
Nconst = 2;

%Data = 4;
Data = 'Data1';

if isstr(Data)
    %Load predefined scenario
    load(Data)
else
    %Create random scenario
    Nagent = Data;
    
    %Define coupling constraints (C-set)
    switch Data
        
        case 4 %Four agents
            Cset = {[1 2],
                    [2 3],
                    [3 4]};
                
        case 3 %Three agents
            Cset = {[1 2],
                    [2 3]};

        case 2 %Two agents
            Cset = {[1 2]};
    end

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
            a{agent}(k,1)   = random('norm',0,1);
        end

    end
end

%return

%%%% Start Algorithm %%%%
for method_number = 1:length(List_of_Methods)
    Method = List_of_Methods{method_number}

    lambda    = zeros(length(Cset),1);
    xFGM      = lambda; %for Fast Gradient Method
    xFGM_prev = xFGM;

    toltau    = 1e-6;
    
    %Naive guess
    for agent = 1:Nagent
        X{agent} = eye(Nvar); Z{agent} = eye(Nvar); 
        mu{agent} = zeros(Nconst,1);
    end    
            
    %Storage for method compare
    all_res_store{method_number}     = [];
    tau_vs_iter_total{method_number} = [];
    iterNT_total{method_number}      = [];
    iterDual_store{method_number}    = [];
    lambda_store{method_number}      = [];
    
    iter_Dual_total = 1; %total number of dual iteration of the algorithm
    dlambda_dtau = 0;tau = 1;
    while tau > toltau

        tolDual = tau;
        tolNT   = tau^2;
        
        tau_vs_iter_total{method_number} = [tau_vs_iter_total{method_number};
                                            iter_Dual_total-1   tau  tolDual tolNT ];
        
        iter_Dual = 1;   %number of dual iterations for the current value of the barrier parameter

       
        
        res = 1e6;lambda_old = lambda;
        while (norm(res,inf) > tolDual)  
            D = 0;res = 0;
            tau_total(iter_Dual_total) = tau;
            for agent = 1:Nagent
                [ X{agent}, Z{agent}, mu{agent}, X_sens{agent}, X_sens_tau{agent}, Z_sens{agent}, Z_sens_tau{agent}, mu_sens{agent}, mu_sens_tau{agent}, iterNT ] = NTSolve(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau, tolNT, X{agent}, Z{agent}, mu{agent}, P{agent}, display_subproblems);
                
                iterNT_total{method_number}(agent,iter_Dual_total) = iterNT;
                
                
                %Dual value
                DC = lambdaC( lambda, C{agent}, P{agent} );
                D = D + trace((Q{agent}  + WeightNuclear*eye(Nvar)  + DC)*X{agent}) - tau*log(det(X{agent}));

                %Dual Residual 
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
                            DH(m,l) = DH(m,l) + trace(C{k}(:,:,l)*X_sens{k}(:,:,m));                            
                        end
                    end
                end
            end

            %Compute D_lambda_tau
            Dtau = zeros(length(Cset),1);
            for l = 1:length(Cset)
                for i = 1:length(Cset{l})
                    k = Cset{l}(i);                     
                    Dtau(l) = Dtau(l) + trace(C{k}(:,:,l)*X_sens_tau{k});                            
                end
            end
            dlambda_dtau = -DH\Dtau;

            Lipschitz = max(abs(eig(DH)));

            %Store
            res_store(iter_Dual)                          = norm(res,inf);
            all_res_store{method_number}(iter_Dual_total) = norm(res,inf);


            %%%%%%%%%%%%%%%%
            %  Dual Update %
            %%%%%%%%%%%%%%%%

            switch Method
                case 'GA' %Gradient ascent
                    lambda  = lambda + res/Lipschitz;
                    
                case 'FGM' %Fast Gradient with restart
                    restart = 0;lambda_old = lambda;
                    if (iter_Dual > 1)
                        if (res_store(iter_Dual) > res_store(iter_Dual-1))
                            display('restart')
                            restart   = 1;
                            xFGM_prev = xFGM; % restart inertia
                            lambda    = xFGM; % update lambda according to gradient ascent on the previous lambda
                        else                            
                            xFGM    = lambda + res/Lipschitz;
                        end
                    else
                        xFGM        = lambda + res/Lipschitz;
                    end
                    
                    if (restart == 0)
                        lambda    = xFGM + (iter_Dual-1)/(iter_Dual+2)*(xFGM - xFGM_prev);
                        xFGM_prev = xFGM;
                    end
                    
                case 'Second_Predictor'   %Second order update 
                    dlambda = - DH\res;
                    lambda  = lambda + dlambda;
                    delta_func = sqrt(res.'*dlambda)/sqrt(tau);
                    
                    %Update local variables based on dlambda
                    for agent = 1:Nagent
                        dmu = 0;dX  = 0;dZ  = 0;

                        for k = 1:length(P{agent})
                            dmu = dmu + mu_sens{agent}(:,P{agent}(k))*dlambda(P{agent}(k)) ;
                            dX  = dX  + X_sens{agent}(:,:,P{agent}(k))*dlambda(P{agent}(k));
                            dZ  = dZ  + Z_sens{agent}(:,:,P{agent}(k))*dlambda(P{agent}(k));
                        end

                        alphaX = alphaPos( X{agent}, dX );
                        alphaZ = alphaPos( Z{agent}, dZ );

                        mu{agent} = mu{agent} + dmu;
                        X{agent}  = X{agent}  + alphaX*dX;
                        Z{agent}  = Z{agent}  + alphaZ*dZ;

                        PosCheck( X{agent}, Z{agent} );
                    end
                    
                case 'Second_NoPredictor'   %Second order update
                    dlambda = - DH\res;
                    lambda  = lambda + dlambda;
                    delta_func = sqrt(res.'*dlambda)/sqrt(tau);
            end
            lambda_store{method_number} = [lambda_store{method_number} lambda];
                       
            
            
            %Next iterate
            iter_Dual       = iter_Dual + 1;
            iter_Dual_total = iter_Dual_total + 1;


            
          
        end
        
        iterDual_store{method_number} = [iterDual_store{method_number};tau  iter_Dual-1];
        
        tau_vs_iter_total{method_number} = [tau_vs_iter_total{method_number};
                                            iter_Dual_total-1   tau  tolDual tolNT];
       
        
        %Update tau
        dtau = -0.5*tau;
        
        
        tau  = tau + dtau;

      
        
        
        %Update of the dual variable for the coming tau update
        if (tau > toltau) 
            switch Method
                case 'Second_Predictor'   
                    %Update primal-dual after barrier update using predictors
                    dlambda = dlambda_dtau*dtau;
                    lambda = lambda + dlambda;   
                    
                    %Update local variables based on dtau AND dlambda
                    for agent = 1:Nagent
                        dmu = mu_sens_tau{agent}*dtau;
                        dX  = X_sens_tau{agent}*dtau;
                        dZ  = Z_sens_tau{agent}*dtau;

                        for k = 1:length(P{agent})
                            dmu = dmu + mu_sens{agent}(:,P{agent}(k))*dlambda(P{agent}(k)) ;
                            dX  = dX  + X_sens{agent}(:,:,P{agent}(k))*dlambda(P{agent}(k));
                            dZ  = dZ  + Z_sens{agent}(:,:,P{agent}(k))*dlambda(P{agent}(k));
                        end

                        alphaX = alphaPos( X{agent}, dX );
                        alphaZ = alphaPos( Z{agent}, dZ );

                        mu{agent} = mu{agent} + dmu;
                        X{agent}  = X{agent}  + alphaX*dX;
                        Z{agent}  = Z{agent}  + alphaZ*dZ;

                        PosCheck( X{agent}, Z{agent} );

                    end 
                case 'Second_NoPredictor'
                    
                    %Update local variables based on dtau only (classical
                    %predictor-corrector)
                    for agent = 1:Nagent
                        dmu = mu_sens_tau{agent}*dtau;
                        dX  = X_sens_tau{agent}*dtau;
                        dZ  = Z_sens_tau{agent}*dtau;

                        alphaX = alphaPos( X{agent}, dX );
                        alphaZ = alphaPos( Z{agent}, dZ );

                        mu{agent} = mu{agent} + dmu;
                        X{agent}  = X{agent}  + alphaX*dX;
                        Z{agent}  = Z{agent}  + alphaZ*dZ;
                        
                        PosCheck( X{agent}, Z{agent} );

                    end
                    
            end
        end
            
        
     

        
    end



end

%%%%%   Check solution   %%%%%%

%Coupling Residual 
for const = 1:length(Cset) %walk through the C-sets
    resC2(const,1) = 0; %res for constraint 'const'
    for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
        participating_agent = Cset{const}(index);
        resC2(const,1) = resC2(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
        %display(['Agent = ',num2str(participating_agent),' / Const = ',num2str(const)]);
    end
end

%Primal residual
for agent = 1:Nagent
    for k = 1:Nconst
        resP2{agent} = trace(A{agent}(:,:,k)*X{agent}) - a{agent}(k)  ;
    end        
end

%Dual residual 
DC = lambdaC( lambda, C{agent}, P{agent} );

resD2 = Q{agent} + WeightNuclear*eye(Nvar) + DC - Z{agent};
for k = 1:Nconst
    resD2 = resD2 - mu{agent}(k)*A{agent}(:,:,k);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



LS = '-';FS = 16;


%%%% Compare all methods
fig = figure(2);clf
List_titles = {'NDSDP with predictor','NDSDP w/o predictor','rFGM'};
for method_number = 1:length(List_of_Methods)   
    subplot(3,1,method_number)
    semilogy(all_res_store{method_number},'linestyle',LS,'marker',Marker{method_number},'color','k');hold on
    plot(tau_vs_iter_total{method_number}(1:end,1),tau_vs_iter_total{method_number}(1:end,2),'linestyle','-','color','k','linewidth',1);hold on
    plot(tau_vs_iter_total{method_number}(1:end,1),tau_vs_iter_total{method_number}(1:end,3),'linestyle','--','color','k','linewidth',1);hold on
    ylabel('$$\|\nabla_\lambda D\|\quad\mathrm{and}\quad \tau$$','interpreter','latex','fontsize',FS)
    set(gca,'fontsize',FS)
    title(List_titles{method_number})
    grid on
    axis tight
end

%legend('N-D with predictor','N-D w/o predictor','rFGM','location','best')
xlabel('Dual iteration','fontsize',FS)

if save
    PapPos = get(gcf,'PaperPosition');PapPos(4) = 1.5*PapPos(4);
    set(gcf,'PaperPosition',PapPos)

    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/CompareSolvers',label];
    exportfig(fig, FileName,'color','cmyk')
end


%%%% Compare number of iteration
Marker = {'o','*','x'};
fig = figure(4);clf
%subplot(2,1,1)
for method_number = 1:length(List_of_Methods)
    loglog(iterDual_store{method_number}(:,1),iterDual_store{method_number}(:,2),'linestyle','none','marker',Marker{method_number},'color','k');hold on
end
axis tight
xlabel('$$\tau$$','interpreter','latex','fontsize',FS);
ylabel(['#Dual iteration'],'fontsize',FS)
grid on    
legend('NDSDP with predictor','NDSDP w/o predictor','rFGM','location','east')

axis tight
set(gca,'XDir','reverse','fontsize',FS);


if save
    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/CompareSolvers_iterations',label];
    exportfig(fig, FileName,'color','cmyk')
end



fig = figure(6);clf
for method_number = 1:length(List_of_Methods)
    sumIterIP(method_number) = sum(max(iterNT_total{method_number}));
end  
barh(sumIterIP,'k')
set(gca,'YTickLabel',{'NDSDP with predictor','NDSDP w/o predictor','rFGM'},'fontsize',2*FS)
xlabel('# Nesterov-Todd steps','fontsize',2*FS)
%ylabel('Method')
grid on
axis tight
if save
    PapPos = get(gcf,'PaperPosition');PapPos(3) = 2*PapPos(3);
    set(gcf,'PaperPosition',PapPos)
    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/CompareSolvers_NTsteps',label];
    exportfig(fig, FileName,'color','cmyk')
end