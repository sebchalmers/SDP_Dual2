clear all
close all
clc

List_of_Methods = {'Second_Dual_Predictor','Second_All_Predictor','Adaptive'};
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
            


    gamma = 0.7;
      
    dlambda_dtau = 0;
    tau = 1;
    res = 1;
    tau_index = 1;
    
    %%%% Initialize on the central path (step 1)
    for agent = 1:Nagent    
        [ X{agent}, Z{agent}, mu{agent}, X_sens{agent}, X_sens_tau{agent}, Z_sens{agent}, Z_sens_tau{agent}, mu_sens{agent}, mu_sens_tau{agent}, Store{method_number}.iterNT(agent) ] = NTSolveMehrotra(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau, tau, X{agent}, Z{agent}, mu{agent}, P{agent}, display_subproblems);
    end
    
    while tau > toltau || norm(res) > toltau
        
        Store{method_number}.tau(tau_index) = tau;
   
        %%%% Step NT system (step 2)
            
        %Initiate res and D for summation
        res = 0;D = 0;
        for agent = 1:Nagent
            [ X{agent}, Z{agent}, mu{agent}, X_sens{agent}, X_sens_tau{agent}, Z_sens{agent}, Z_sens_tau{agent}, mu_sens{agent}, mu_sens_tau{agent} ] = NTStepMehrotra(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau,  X{agent}, Z{agent}, mu{agent}, P{agent});
            


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
            
        Store{method_number}.res(tau_index) = norm(res);
        Store{method_number}.lambda(:,tau_index) = lambda;
            
        %%%% Dual gradient and Hessian (step 3)
        
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

        %Compute gradient of D w.r.t barrier tau
        Dtau = zeros(length(Cset),1);
        for l = 1:length(Cset)
            for i = 1:length(Cset{l})
                k = Cset{l}(i);                     
                Dtau(l) = Dtau(l) + trace(C{k}(:,:,l)*X_sens_tau{k});                            
            end
        end

        %%%% Update barrier tau (step 4)  
        if (tau <= toltau)
            gamma = 1;
        end
        
                %Update tau

        switch Method
            case 'Adaptive'
                %Feedback on gamma
                if (tau < 1)
                    Kp = 0.02;Kd = -0.01;
                else
                    Kp = 0;Kd = 0;
                end
                
                Error = (norm(res)-sqrt(tau))/sqrt(tau);
                
                Pcorr = Kp*Error;
                if (tau < 1)
                    delta_Error = Error-Error_prev;
                    Dcorr = Kd*delta_Error;
                else
                    Dcorr = 0;
                end
                gamma = max(0.5,min(1,gamma + Pcorr + Dcorr));

                Error_prev = Error;

        end
        
        
        tau_new = gamma*tau;
        dtau = tau_new - tau;
        tau = tau_new;
  


        %%%%  Dual Update  (step 5)  

        dlambda = - DH\(res + Dtau*dtau);

        %Update dual variables with predictor
        lambda = lambda + dlambda;
        
        %%%% Predictor on the local solutions (step 6)
        switch Method
            case 'Second_Dual_Predictor'
                % No local update
                
            case 'Second_All_Predictor'

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
                
        end

       
                                           
    
        tau_index = tau_index + 1;
        
            
    end
    Store{method_number}.iter = tau_index-1;


end

LS = '-';FS = 16;



%%%% Compare all methods, version 1
fig = figure(1);clf
 
for method_number = 1:length(List_of_Methods)
    loglog(Store{method_number}.tau,Store{method_number}.res,'linestyle',LS,'marker',Marker{method_number},'color','k');hold on
end
legend({'Algorithm 2, w/o step 7','Algorithm 2'},'location','southwest')
loglog(Store{method_number}.tau,Store{method_number}.tau,'linestyle','--','color','k');hold on

line([10 tau],[toltau toltau],'linestyle','--','color','k')
line([toltau toltau],[1 toltau],'linestyle','--','color','k')

xlabel('$$\tau$$','interpreter','latex','fontsize',FS)
ylabel('$$\|\nabla_\lambda D\|$$','interpreter','latex','fontsize',FS)
set(gca,'fontsize',FS)
%title(List_titles{method_number})
grid on
axis tight
xlim([tau/2 1])
ylim([toltau/2 2])

set(gca,'XDir','reverse','fontsize',FS);

title('NDSDP Path-Following, barrier and dual residual')

if save
    PapPos = get(gcf,'PaperPosition');%PapPos(4) = 1.5*PapPos(4);
    set(gcf,'PaperPosition',PapPos)

    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/PathFollowing_tau_res',label];
    exportfig(fig, FileName,'color','cmyk')
end

%%%% Compare all methods, version 2
fig = figure(2);clf
 
for method_number = 1:length(List_of_Methods)
    semilogy(Store{method_number}.res,'linestyle',LS,'marker',Marker{method_number},'color','k');hold on
end
semilogy(Store{method_number}.tau,'linestyle','-','color','k');hold on
legend({'Algorithm 2, w/o step 7\qquad','Algorithm 2','Barrier parameter $$\tau$$'},'interpreter','latex','fontsize',FS,'location','southwest')


line([1 Store{method_number}.iter],[toltau toltau],'linestyle','--','color','k')
%ylabel('$$\|\nabla_\lambda D\|$$ and $$\tau$$','interpreter','latex','fontsize',FS)

%title(List_titles{method_number})
grid on
axis tight

xlabel('Iteration','fontsize',FS)


title('NDSDP Path-Following, barrier and dual residual')
set(gca,'fontsize',FS);

if save
    PapPos = get(gcf,'PaperPosition');%PapPos(4) = 1.5*PapPos(4);
    set(gcf,'PaperPosition',PapPos)

    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/PathFollowing_tau_and_res',label];
    exportfig(fig, FileName,'color','cmyk')
end



%%% Dual variables
fig = figure(4);clf
for method_number = 1:length(List_of_Methods)
        plot(Store{method_number}.lambda(1,:),Store{method_number}.lambda(2,:),'marker',Marker{method_number},'color','k','linestyle',LS);hold on
end
grid on
xlabel('$$\lambda_1$$','interpreter','latex','fontsize',FS)
ylabel('$$\lambda_2$$','interpreter','latex','fontsize',FS)
set(gca,'fontsize',FS)
legend({'Algorithm 2, w/o step 7\qquad','Algorithm 2'},'interpreter','latex','location','southeast')
title('NDSDP Path-Following, trajectory of the dual variables')

if save
    PapPos = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',PapPos)
    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/PathFollowing_Dual',label];
    exportfig(fig, FileName,'color','cmyk')
end

