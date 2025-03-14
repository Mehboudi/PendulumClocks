%% Load and Plot the asymptotic limit cycles
figure
myVars0={'M_vec'};
sub_folder_name='M1';
load([sub_folder_name,'/in_cond_M1','traj1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
iMmax=length(M_vec);
for iM=1:iMmax
    M=M_vec(1,iM);
    sub_folder_name=['M',num2str(iM)];
    myVars = {'x_m_vec','p_m_vec'};
    load([sub_folder_name,'/in_cond_M',num2str(iM)],myVars{:})
    subplot(2,4,iM);
    plot(-x_m_vec/sqrt(2),1i*p_m_vec/sqrt(2),'LineWidth',2)
    fontsize(20,"points")
    set(gca,'linewidth',1)
    %ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
    %'Interpreter','latex','FontSize', 20);
    %xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 20);
    title(['$M~ = ~$',num2str(iM)],'FontSize', 20,'Interpreter','latex');
    %xlim([-150,150]/sqrt(2))
    %ylim([-150,150]/sqrt(2))
    [iM iMmax]
end
%% %% Load and Plot the asymptotic limit cycles conditional
figure
myVars0={'M_vec'};
sub_folder_name='M1';
load([sub_folder_name,'/in_cond_M1','traj1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
iMmax=length(M_vec);
for iM=1:iMmax
    M=M_vec(1,iM);
    sub_folder_name=['M',num2str(iM)];
    myVars = {'x_m_vec','p_m_vec'};
    load([sub_folder_name,'/in_cond_M',num2str(iM),'traj1'],myVars{:})
    subplot(3,3,iM);
    plot(-x_m_vec/sqrt(2),1i*p_m_vec/sqrt(2),'LineWidth',2)
    fontsize(20,"points")
    set(gca,'linewidth',1)
    %ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
    %'Interpreter','latex','FontSize', 20);
    %xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 20);
    title(['$M~ = ~$',num2str(iM)],'FontSize', 20,'Interpreter','latex');
    %xlim([-150,150]/sqrt(2))
    %ylim([-150,150]/sqrt(2))
    [iM iMmax]
end
%%%
%%%
%% The real deal (histogram etc)
figure
myVars0={'M_vec'};
sub_folder_name='M2';
load([sub_folder_name,'/in_cond_M2','traj1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
iMmax=length(M_vec);
%You only plot the following indices (i.e., number of copies). You can modify them.
plotindices=[1,2,4,6];splot=length(plotindices);
scount=0;
t = tiledlayout(1, splot, 'TileSpacing', 'compact', 'Padding', 'compact');
for iM=1:iMmax
    if ismember(iM, plotindices)
        scount=scount+1;
        nexttile(scount);
        M=M_vec(1,iM);
        sub_folder_name=['M',num2str(iM)];
        myVars = {'x_m_vec','p_m_vec'};
        load([sub_folder_name,'/in_cond_M',num2str(iM),'traj1'],myVars{:})
        %subplot(1,splot,scount);
        histogram2(-[110,-110,x_m_vec/sqrt(2),-110,110],...
            [110,-110,1i*p_m_vec/sqrt(2),-110,110],[200,200],...
            'DisplayStyle', 'tile', 'ShowEmptyBins', 'on','EdgeColor','None'...
            )%I add the new points so that the axis limits are equal everywhere
        %histogram2(-x_m_vec/sqrt(2),...
        %1i*p_m_vec/sqrt(2),...
        %'DisplayStyle', 'tile', 'ShowEmptyBins', 'on')%I add the new points so that the axis limits are equal everywhere
        fontsize(20,"points")
        set(gca,'linewidth',1)
        colormap('parula')
        %ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
        %'Interpreter','latex','FontSize', 20);
        %xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 20);
        title(['$M~ = ~$',num2str(iM)],'FontSize', 20,'Interpreter','latex');
        %colorbar;
        %clim([0 .5]);
        %xLimits{1,iM} = xlim;
        %yLimits{1,iM} = ylim;
        hold on
        %xlim([-100,100])
        %ylim([-100,100])
        [iM iMmax]
    end
end
%%
%%%
%%%
% plot the limit cycle too
%%%
%%%
myVars0={'M_vec'};
sub_folder_name='M2';
%load([sub_folder_name,'/in_cond_M2','traj1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
scount=0;
for iM=1:iMmax
    if ismember(iM, plotindices)
        scount=scount+1;
        nexttile(scount);
        M=M_vec(1,iM);
        sub_folder_name=['M',num2str(iM)];
        myVars = {'x_m_vec','p_m_vec'};
        load([sub_folder_name,'/in_cond_M',num2str(iM)],myVars{:})
        %subplot(1,splot,scount);
        plot(-x_m_vec/sqrt(2),1i*p_m_vec/sqrt(2),'LineWidth',2)
        %xlim(xLimits{1,iM});
        %ylim(yLimits{1,iM});
        %fontsize(20,"points")
        set(gca,'linewidth',1)
        xlabel('$x_{\rm m}$','Interpreter','latex')
        if iM<2
            ylabel('$p_{\rm m}$','Interpreter','latex')
        end
        [iM iMmax]
    end
end
%%%
%%%
% Use tightfig to remove extra white space
tightfig;
