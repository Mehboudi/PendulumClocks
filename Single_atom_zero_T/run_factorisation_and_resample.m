%% Simulate trajectories.
% % If you have done it already, comment for further analysis below
clear all
iM=1;
imin=1;
imax=10;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
n_c=0;T_c=0;
sub_folder_name='Data';
mkdir(sub_folder_name)
for ur=0:1
    if ur==0
        i1=1;%don't touch, this is called in Factorisation; 
        Factorisation;
        myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
            'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_hot','w_cold','w_cav','n_h','n_c','w_m','f','g'...
            ,'k','g_h','g_c','g_m','dt'};
        save([sub_folder_name,'/unconditional'],myVars{:});
        % plot(x_m_vec,1i*p_m_vec,'LineWidth',2);
        % xlim([-40 40])
        % ylim([-50 50])
    else
        for i1=imin:imax
            Factorisation;
            tvec_dN1=jump_times;
            myVars2={"tvec_dN1","p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
                'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_hot','w_cold','w_cav',...
                'n_h','n_c','w_m','f','g','k','g_h','g_c','g_m','dt'};
            save([sub_folder_name,'/conditional_traj',num2str(i1)],myVars2{:});
            [i1,imax]
        end
    end
end
%%%
iM=1;%This is the number of atoms, here it is one, cannot be changed. 
% See multiple copies folder if you wanna change it. It shows up in some dependent codes
%%%
%% Load and Plot the asymptotic limit cycles; the trajectories for the populations and ticks (vlines)
%%%%The uncond
myVars = {'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_m'};
sub_folder_name='Data';
mkdir('Data/Pics')
load([sub_folder_name,'/unconditional'],myVars{:})
%%%%reduce size
ds=1e3;%downsample
x_m_vec=x_m_vec(1:ds:end);
p_m_vec=p_m_vec(1:ds:end);
p1_vec=p1_vec(1:ds:end);
p2_vec=p2_vec(1:ds:end);
na_vec=na_vec(1:ds:end);
t_vec_i1=t_vec_i1(1:ds:end);
%%%%Now start plotting
figure
xlim([-30 30])
ylim([-30 30])
plot(x_m_vec,1i*p_m_vec,'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
    'Interpreter','latex','FontSize', 20);
xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 20);
title('\rm Unconditional evolution','FontSize', 22);
saveas(gcf,[pwd '/Data/Pics/Phase_space_uncond.png'])
saveas(gcf,[pwd '/Data/Pics/Phase_space_uncond.fig'])
%%%
figure
xlim([-30 30])
ylim([-30 30])
plot(t_vec_i1(2:end)*w_m/pi,[x_m_vec(2:end);1i*p_m_vec(2:end)]*w_m/pi,'LineWidth',2)
xlim([t_vec_i1(1,2)*w_m/pi,t_vec_i1(1,end)*w_m/pi])
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('Mech. osc. quadratures','FontSize', 20);
xlabel('$\omega_m t/\pi$','Interpreter','latex','FontSize', 20);
title('\rm Unconditional evolution','FontSize', 22);
legend('$\left\langle b+b^{\dagger} \right\rangle$',...
    '$i\left\langle b-b^{\dagger} \right\rangle$','Interpreter','latex')
saveas(gcf,[pwd '/Data/Pics/quadratures_uncond.png'])
saveas(gcf,[pwd '/Data/Pics/quadratures_uncond.fig'])
%%%
figure
p3_vec=1-p1_vec-p2_vec;
plot(t_vec_i1(2:end)*w_m/pi,[p1_vec(2:end);p2_vec(2:end);p3_vec(2:end)],'LineWidth',2)
xlim([t_vec_i1(1,2)*w_m/pi,t_vec_i1(1,end)*w_m/pi])
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('Populations','FontSize', 20);
xlabel('$\omega_m t/\pi$','Interpreter','latex','FontSize', 20);
title('\rm Unconditional evolution','FontSize', 22);
legend('$p_1$','$p_2$','$p_3$','Interpreter','latex')
saveas(gcf,[pwd '/Data/Pics/populations_uncond.png'])
saveas(gcf,[pwd '/Data/Pics/populations_uncond.fig'])
%%%
figure
plot(t_vec_i1(2:end)*w_m/pi,na_vec(2:end),'LineWidth',2)
xlim([t_vec_i1(1,2)*w_m/pi,t_vec_i1(1,end)*w_m/pi])
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('Cavity occupation','FontSize', 20);
xlabel('$\omega_m t/\pi$','Interpreter','latex','FontSize', 20);
title('\rm Unconditional evolution','FontSize', 22);
saveas(gcf,[pwd '/Data/Pics/cav_occu_uncond.png'])
saveas(gcf,[pwd '/Data/Pics/cav_occu_uncond.fig'])
%%
%%%%The cond
%%%%
myVars = {'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_m','tvec_dN1'};
load([sub_folder_name,'/conditional_traj',num2str(1)],myVars{:})
%%%%reduce size
ds=1e3;%downsample
x_m_vec=x_m_vec(1:ds:end);
p_m_vec=p_m_vec(1:ds:end);
p1_vec=p1_vec(1:ds:end);
p2_vec=p2_vec(1:ds:end);
na_vec=na_vec(1:ds:end);
t_vec_i1=t_vec_i1(1:ds:end);
%this next one is tricky, as it is all the lines.
tvec_dN1=tvec_dN1(tvec_dN1>t_vec_i1(1,2));
tvec_dN1=tvec_dN1(tvec_dN1<t_vec_i1(1,end));
%%%%Now start plotting
figure
xlim([-30 30])
ylim([-30 30])
plot(-x_m_vec/sqrt(2),1i*p_m_vec(sqrt(2)),'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
    'Interpreter','latex','FontSize', 20);
xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 20);
title('\rm Conditional evolution','FontSize', 22);
saveas(gcf,[pwd '/Data/Pics/Phase_space_cond.png'])
saveas(gcf,[pwd '/Data/Pics/Phase_space_cond.fig'])
%%%
figure
xlim([-30 30])
ylim([-30 30])
plot(t_vec_i1(2:end)*w_m/pi,[-x_m_vec(2:end)/sqrt(2);...
    1i*p_m_vec(2:end)/sqrt(2)]*w_m/pi,'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
hold on
xline(tvec_dN1(2:end)*w_m/pi,'-k','FontSize',2,'HandleVisibility','off')
xlim([t_vec_i1(1,2)*w_m/pi,t_vec_i1(1,end)*w_m/pi])
ylabel('Mech. osc. quadratures','FontSize', 20);
xlabel('$\omega_m t/\pi$','Interpreter','latex','FontSize', 20);
title('\rm Conditional evolution','FontSize', 22);
legend('$\left\langle b+b^{\dagger} \right\rangle$',...
    '$i\left\langle b-b^{\dagger} \right\rangle$','Interpreter','latex')
saveas(gcf,[pwd '/Data/Pics/quadratures_cond.png'])
saveas(gcf,[pwd '/Data/Pics/quadratures_cond.fig'])
%%%
figure
p3_vec=1-p1_vec-p2_vec;
plot(t_vec_i1(2:end)*w_m/pi,[p1_vec(2:end);p2_vec(2:end);p3_vec(2:end)],'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
hold on
xline(tvec_dN1(2:end)*w_m/pi,'-k','FontSize',2,'HandleVisibility','off')
xlim([t_vec_i1(1,2)*w_m/pi,t_vec_i1(1,end)*w_m/pi])
ylabel('Populations','FontSize', 20);
xlabel('$\omega_m t/\pi$','Interpreter','latex','FontSize', 20);
title('\rm Conditional evolution','FontSize', 22);
legend('$p_1$','$p_2$','$p_3$','Interpreter','latex')
saveas(gcf,[pwd '/Data/Pics/populations_cond.png'])
saveas(gcf,[pwd '/Data/Pics/populations_cond.fig'])
%%%
figure
plot(t_vec_i1(2:end)*w_m/pi,na_vec(2:end),'LineWidth',2)
fontsize(20,"points")
set(gca,'linewidth',1)
hold on
xline(tvec_dN1*w_m/pi,'-k','FontSize',2,'HandleVisibility','off')
xlim([t_vec_i1(1,2)*w_m/pi,t_vec_i1(1,end)*w_m/pi])
ylabel('Cavity occupation','FontSize', 20);
xlabel('$\omega_m t/\pi$','Interpreter','latex','FontSize', 20);
title('\rm Conditional evolution','FontSize', 22);
saveas(gcf,[pwd '/Data/Pics/cav_occu_cond.png'])
saveas(gcf,[pwd '/Data/Pics/cav_occu_cond.fig'])
%%%
% myVars = {'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_m','tvec_dN1'};
% load([sub_folder_name,'/conditional_traj',num2str(1)],myVars{:})
plot_filter=1;
Detector_Filter_saturation;
s1=size(tvec_dN1)
s2=size(tvec_dN1_I2)
%%%
plot_filter=0;%do not coment this
%% Load and analyze the data: Tick states
% %-----------------
% %   Tick stats; No filter
% %-----------------
imin=1;
%imax=1;
det_filt=0;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
sub_folder_name='Data';
for i1=imin:1:imax
    myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','n_c'};
    load([sub_folder_name,'/conditional_traj',num2str(i1)],myVars{:})
    %Let's renormalise everything!
    tvec_dN1=tvec_dN1*w_m/pi;
    %%%%This line will be passed only if you want to filter (detector dead time)
    if det_filt==1
        Detector_Filter_saturation;
        tvec_dN1=tvec_dN1_I2(1:end);
    end
    %%%%Otherwise carryout as usual
    dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
    size(dtjump);
    dtj=[dtj,dtjump];
    muvec(1,i1)=mean(dtjump);
    varvec(1,i1)=std(dtjump)^2;
end
figure
mu_=mean(muvec)
var_=mean(varvec)%Note we take mean of the var over different rounds.
N=mu_.^2./var_
% % %%
bin=200;
histogram(dtj(2:end),bin,Normalization="probability",EdgeColor="none")
% Create xline
xline([0 1 2 3 4],'-k','FontSize',2,'HandleVisibility','off');
tname=(['~~${\cal N}=$',num2str(N),'$~~~\nu=$',num2str(1/mu_),'$\Omega_{\rm m}/\pi$' ]);
title(tname,'Interpreter','latex','FontSize',18)
fontsize(20,"points")
set(gca,'linewidth',1)
xlim([-0.1,3])
xlabel('$\Omega_{\rm m}\tau /\pi$','FontSize', 20,'Interpreter','latex');
%saveas(gcf,[pwd '/Data/Pics/Histogram_no_filter.png'])
%saveas(gcf,[pwd '/Data/Pics/Histogram_no_filter.fig'])
% Allan
% xlim([1,2e3])
% saveas(gcf,[pwd '/Data/Pics/Allan_no_filter_scrach.png'])
% saveas(gcf,[pwd '/Data/Pics//Allan_no_filter_scrach.fig'])
%% Now, Take the ticks, and consider the dead-time of the detectors
% %-----------------
% %   Tick stats; Filter
% %-----------------
plot_filter=1;
det_filt=1;
imin=1;
%imax=100;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
sub_folder_name='Data';
for i1=imin:1:imax
    myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav'};
    load([sub_folder_name,'/conditional_traj',num2str(i1)],myVars{:})
    %%%%This line will be passed only if you want to filter (detector dead time)
    %%%The detector parameters are set on the other code, check it out
    %%%there.
    if det_filt==1
        Detector_Filter_saturation;
        tvec_dN1=tvec_dN1_I2(1:end);
    end
    %Let's renormalise everything!
    tvec_dN1=tvec_dN1*w_m/pi;
    %%%%Otherwise carryout as usual
    dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
    dtj=[dtj,dtjump(2:end)];
    muvec(1,i1)=mean(dtjump(2:end));
    varvec(1,i1)=std(dtjump(2:end))^2;
    [i1 imax]
end
% dtj=dtj(2:end);
% res=mean(dtj)
% res_stable=mean(dtj_stable)
% std_=std(dtj)
% std_stable=std(dtj_stable)
%%%IMPORTANT: This is a test, change back to the immediate lower two lines!!!
mu_=mean(muvec,'omitnan')
var_=mean(varvec,'omitnan')
N=mu_.^2./var_
bin=200;
figure
histogram(dtj(2:end),bin,Normalization="probability",EdgeColor="none")
% Create xline
xline([0 1 2 3 4],'-k','FontSize',2,'HandleVisibility','off');
tname=(['~~${\cal N}=$',num2str(N),'$~~~\nu=$',num2str(1/mu_),'$\Omega_{\rm m}/\pi$' ]);
title(tname,'Interpreter','latex','FontSize',18)
fontsize(20,"points")
set(gca,'linewidth',1)
xlabel('$\Omega_{\rm m}\tau/ \pi$','FontSize', 20,'Interpreter','latex');
xlim([-.1,3])
%saveas(gcf,[pwd '/Data/Pics/Histogram_filter.png'])
%saveas(gcf,[pwd '/Data/Pics/Histogram_filter.fig'])
Allan
xlim([1,1e3])
%saveas(gcf,[pwd '/Data/Pics/Allan_filter_scrach.png'])
%saveas(gcf,[pwd '/Data/Pics//Allan_filter_scrach.fig'])
