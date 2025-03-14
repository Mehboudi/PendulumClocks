%%%Warning, this is a slow code, only for phase space limit cycle!!
%% Simulate trajectories.
% % If you have done it already, comment for further analysis
clear all
imin=1;
imax=1;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
n_c=0;
sub_folder_name='Data';
mkdir(sub_folder_name)
for ur=0:1
    if ur==0
        Factorisation;
        myVars = {'x_m_vec','p_m_vec','w_hot','w_cold','w_cav','n_h','n_c','w_m','f','g'...
            ,'k','g_h','g_c','g_m','dt'};
        save([sub_folder_name,'/unconditional'],myVars{:});
        % plot(x_m_vec,1i*p_m_vec,'LineWidth',2);
        % xlim([-40 40])
        % ylim([-50 50])
    else
        for i1=imin:imax
            [i1,imax]
            Factorisation;
            tvec_dN1=jump_times;
            myVars2={'x_m_vec','p_m_vec','w_hot','w_cold','w_cav',...
                'n_h','n_c','w_m','f','g','k','g_h','g_c','g_m','dt'};
            save([sub_folder_name,'/conditional_traj',num2str(i1)],myVars2{:});
        end
    end
end
%% Load and Plot the asymptotic limit cycles
figure
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/unconditional'],myVars{:})
%subplot(1,4,iM);
plot(-x_m_vec/sqrt(2),1i*p_m_vec/sqrt(2),'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('$p_{\rm m}$','Interpreter','latex','FontSize', 20);
xlabel('$x_{\rm m}$','Interpreter','latex','FontSize', 20);
%% %% Load and Plot the asymptotic limit cycles conditional
figure
M=M_vec(1,iM);
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/conditional_traj1'],myVars{:})
endlim=floor(.01*length(x_m_vec));
%subplot(1,4,iM);
plot(-x_m_vec(1:endlim)/sqrt(2),1i*p_m_vec(1:endlim)/sqrt(2),'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('$p_{\rm m}$','Interpreter','latex','FontSize', 20);
xlabel('$x_{\rm m}$','Interpreter','latex','FontSize', 20);
%%%
%%%
%% The real deal (histogram etc)
figure
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/conditional_traj1'],myVars{:})
%subplot(1,6,iM);
histogram2(-[30,-30,x_m_vec/sqrt(2),-30,30],...
    [30,-30,1i*p_m_vec/sqrt(2),030,30],[100,100],...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'on','EdgeColor','None')%I add the new points so that the axis limits are equal everywhere
%histogram2(-x_m_vec/sqrt(2),...
%1i*p_m_vec/sqrt(2),...
%'DisplayStyle', 'tile', 'ShowEmptyBins', 'on')%I add the new points so that the axis limits are equal everywhere
fontsize(20,"points")
set(gca,'linewidth',1)
colormap('parula')
ylabel('$p_{\rm m}$','Interpreter','latex','FontSize', 20);
xlabel('$x_{\rm m}$','Interpreter','latex','FontSize', 20);
xLimits{1,iM} = xlim;
yLimits{1,iM} = ylim;
hold on
%%
%%%
%%%
% plot the limit cycle too
%%%
%%%
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/unconditional'],myVars{:})
%subplot(1,6,iM);
plot(-x_m_vec/sqrt(2),1i*p_m_vec/sqrt(2),'LineWidth',2)
%xlim(xLimits{1,iM});
%ylim(yLimits{1,iM});
%fontsize(20,"points")
set(gca,'linewidth',1)
%%%
%%%
% Use tightfig to remove extra white space
xlim([-30,30])
ylim([-30,30])
tightfig;
