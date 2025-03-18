%Run this code section by section if you know what you are doing. That
%would allow you to control, chose, what you are doing. Specially, the
%first section below is the most time consuming, generates data, and saves it into your computer. Later
%sections are analysing this data. You can come back to your saved data
%anytime, just make sure you don't get any errors when saving. A priori,
%you should not get errors, but you may need to change file directories
%etc.
%% Simulate trajectories.
% This simulates trajectories (imax) for each cold bath occupation (n_c_vec)
% It first runs the code for a longer time, without unraveling. This will
% be used for (1) the steady state and limit cycle properties, heat
% currents etc. (2) the start point (initial conditions) of the unravelled trajectories. 
% % If you have done it already, comment for further analysis
clear all
iM=1;
imin=1;
imax=1;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
%iTmax=12;n_c_vec=[logspace(-5,0,iTmax/2),linspace(1,n_h,iTmax/2)];
iTmax=20;
n_c_vec=[logspace(-5,-2,iTmax)];
n_c_vec=unique(n_c_vec);
iTmax=length(n_c_vec);
figure
for iT=1:iTmax
    n_c=n_c_vec(1,iT);
    T_c=w_cold/(log((n_c+1)/n_c));
    sub_folder_name=['n_c',num2str(iT)];
    mkdir(sub_folder_name)
    for ur=0:1
        if ur==0
            Factorisation;
            myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
                'x_m_vec','p_m_vec','J_h','J_m','J_cold','J_cav','w_hot','w_cold','w_cav','n_c_vec'};
            save([sub_folder_name,'/in_cond_n_c',num2str(iT)],myVars{:});
            subplot(9,5,iT);
            plot(x_m_vec,1i*p_m_vec,'-*','LineWidth',1);
            xlim([-40 40])
            ylim([-50 50])
            title('nc=',num2str(n_c));
            [iT iTmax]
        else
            for i1=imin:imax
                Factorisation;
                [i1,imax; iT, iTmax]
                tvec_dN1=jump_times;
                myVars2={"tvec_dN1","w_m",'w_hot','w_cold','w_cav',"n_h", 'Q_h','Q_h_f',...
                    'J_h','J_m','J_cold','J_cav','n_c_vec'};
                save([sub_folder_name,'/in_cond_n_c',num2str(iT),'traj',num2str(i1)],myVars2{:});
            end
        end
    end
end
%% Load and Plot the asymptotic limit cycles
% This plots limit cycles of the mechanical resonator.
figure
myVars0={'n_c_vec'};
sub_folder_name='n_c1';
load([sub_folder_name,'/in_cond_n_c1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
iTmax=length(n_c_vec);
icntr=1;
for iT=1:1:iTmax
    n_c=n_c_vec(1,iT);
    sub_folder_name=['n_c',num2str(iT)];
    myVars = {'x_m_vec','p_m_vec'};
    load([sub_folder_name,'/in_cond_n_c',num2str(iT)],myVars{:})
    subplot(5,4,iT);icntr=icntr+1;
    %plot(x_m_vec,1i*p_m_vec,'-*','LineWidth',1);
    plot(-x_m_vec/sqrt(2),1i*p_m_vec/sqrt(2),'linewidth',2);
    xlim([-35 35])
    ylim([-35 35])
    title(['$\bar n_c=$',num2str(n_c)],'Interpreter','latex');
    fontsize(20,"points")
    set(gca,'linewidth',1)
    %ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
        %'Interpreter','latex','FontSize', 10);
    %xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 10);
    [iT iTmax]
end
figure
plot(x_m_vec)
%% Mechanical vs hot bath heat currents
% Title is self explanatory. It specifically can be checked, for the
% parameters in the paper, that J_m<<J_hot
figure
imin=1;
imax=1;
myVars0={'n_c_vec'};
myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
               'x_m_vec','p_m_vec','J_h','J_m','J_cold','J_cav','w_hot','w_cold','w_cav','n_c_vec'};
sub_folder_name='n_c1';
load([sub_folder_name,'/in_cond_n_c1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
sTmax=length(n_c_vec);
J_hot_vec=zeros(1,iTmax);
J_m_vec=zeros(1,iTmax);
J_cold_vec=zeros(1,iTmax);
J_cav_vec=zeros(1,iTmax);
for iT=1:iTmax
    n_c=n_c_vec(1,iT);
    sub_folder_name=['n_c',num2str(iT)];
    myVars = {'J_h','J_m','J_cold','J_cav','w_hot','w_cold','w_cav'};
    load([sub_folder_name,'/in_cond_n_c',num2str(iT)],myVars{:})
    J_hot_vec(1,iT)=J_h;
    J_m_vec(1,iT)=J_m;
    J_cold_vec(1,iT)=J_cold;
    J_cav_vec(1,iT)=J_cav;
end
%plot(n_c_vec,[-J_hot_vec/w_hot;J_cold_vec/w_cold;J_cav_vec/w_cav])
figure
plot(n_c_vec,[-J_hot_vec;J_m_vec]/(w_m*w_hot))
fontsize(10,"points")
set(gca,'linewidth',1)
box on
xlabel('$n_{\rm c}$','Interpreter','latex')
xscale log
legend('$J_{\rm h}/\Omega_{\rm m} \omega_{13}$','$J_{\rm m}/ \Omega_{\rm m} \omega_{13}$','interpreter','latex')
%% Load and analyze the data
% %-----------------
% %   Tick stats; No filter
% %-----------------
% We look at the tick stats, heat currents, entropy production etc.
det_filt=0;
imin=1;
%imax=1;
iTmax=length(n_c_vec);
N=zeros(1,iTmax);
mu_=zeros(1,iTmax);
var_=zeros(1,iTmax);
% Jhmat=zeros(iTmax,imax);
Jmmat=zeros(iTmax,imax);
% Jcoldmat=zeros(iTmax,imax);
% Jcavmat=zeros(iTmax,imax);
% click_num=zeros(iTmax,imax);
figure
scount=0;%counts sub-plot number for histograms.
for iT=1:iTmax
    sub_folder_name=['n_c',num2str(iT)];
    dtj=[];
    muvec=zeros(1,imax);
    varvec=zeros(1,imax);
    %for i1=imin:1:100
    for i1=imin:1:imax
        myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','J_h','J_m','J_cold','J_cav','n_c_vec','Q_h','Q_h_f'};
        load([sub_folder_name,'/in_cond_n_c',num2str(iT),'traj',num2str(i1)],myVars{:})
        %%%%This line will be passed only if you want to filter (detector dead time)
        if det_filt==1
            Detector_Filter_saturation;
            tvec_dN1=tvec_dN1_I2(1:end);
        end
        %Let's renormalise everything!
        tvec_dN1=tvec_dN1*w_m/pi;
        %%%%Otherwise carryout as usual
        dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
        size(dtjump);
        dtj=[dtj,dtjump(2:end)];
        muvec(1,i1)=mean(dtjump);
        varvec(1,i1)=std(dtjump)^2;
        %%%HEAT CURRENT
        % Jhmat(iT,i1)=J_h;
        % Jmmat(iT,i1)=J_m;
        % Jcoldmat(iT,i1)=J_cold;
        % Jcavmat(iT,i1)=J_cav;
        % click_num(iT,i1)=length(tvec_dN1);
        %%%Other currents
    end
    mu_(1,iT)=mean(muvec)
    var_(1,iT)=mean(varvec)%Note we take mean of the var over different rounds.
    N(1,iT)=mu_(1,iT).^2./var_(1,iT);
    % % %%
    bin=40;
    % hold on
    scount=scount+1;
    subplot(4,5,scount)
    histogram(dtj(2:end),bin,'LineStyle','none')
    % Create xline
    xline([0 1 2]);
    tname=(['$\omega_m =$',num2str(w_m),'$n_c = $',num2str(n_c_vec(1,iT)), ...
        '$~~~\mu=$',num2str(mu_),'~~$\sigma^2=$',num2str(var_),'~~$N=$',num2str(N) ]);
    %title(tname,'Interpreter','latex')
    % Create xlabel
    %xlabel('$\omega_m t/\pi$','Interpreter','latex');
    n_c=n_c_vec(1,iT);
    title(['$\bar n_c=$',num2str(n_c)],'Interpreter','latex');
    fontsize(10,"points")
    set(gca,'linewidth',1)
    box on
end
figure
plot(n_c_vec,mu_)
xlabel('$n_c$','Interpreter','latex')
ylabel('$\mu$','Interpreter','latex')
title('Resolution')
xscale log
xlim([n_c_vec(1,1) n_c_vec(1,end)])
fontsize(20,"points")
set(gca,'linewidth',1)
grid on
%%%
figure
plot(n_c_vec,N)
xlabel('$n_c$','Interpreter','latex')
ylabel('$N$','Interpreter','latex')
title('\rm Accuracy')
xscale log
xlim([n_c_vec(1,1) n_c_vec(1,end)])
fontsize(20,"points")
set(gca,'linewidth',1)
grid on
%%%HEAT CURRENT
% Jhmean=mean(Jhmat,2,'omitnan');%The heat current (instantaniuous)
% Jmmean=mean(Jmmat,2,'omitnan');%The heat current (instantaniuous)
% Jcoldmean=mean(Jcoldmat,2,'omitnan');%This is collective heat instead of the instantanious one
% Jcavmean=mean(Jcavmat,2,'omitnan');
% Q_click=-w_hot*mean(click_num,2,'omitnan');%This is the dissipated heat, with some approximation valid for low n_c (for some reason).
% figure
% %plot(n_c_vec,[Jhmean/w_hot,Jcoldmean/w_cold,Jcavmean/w_cav])
% plot(n_c_vec,(Jhmean-J_hot_vec')./Jhmean)
% figure
% plot(n_c_vec,(Jmmean-J_m_vec')./Jmmean)
% %%%%
% xlabel('$n_c$','Interpreter','latex')
% ylabel('$J_\alpha/\omega_\alpha$','Interpreter','latex')
% title('\rm Heat current to frequency rate')
% legend('$J_{\rm hot}/\omega_{\rm hot}$','$J_{\rm cold}/\omega_{\rm cold}$' ...
%     ,'$J_{\rm cav}/\Omega_{\rm cav}$','interpreter','latex')
% xscale log
% xlim([n_c_vec(1,1) n_c_vec(1,end)])
% fontsize(20,"points")
% set(gca,'linewidth',1)
% grid on
%%%
%Jhstd=std(Jhmat','omitnan');
T_c_vec=w_cold./(log((1+n_c_vec)./n_c_vec));
%ent_prod=Jhmean'.*(1./T_h-1./T_c_vec);
ent_prod=J_hot_vec.*(1./T_h-1./T_c_vec);%J_hot_vec comes from the previous code (unconditional evolution)
% figure
% hold on
% fill([n_c_vec, flip(n_c_vec)], [Jhmean'+Jhstd, flip(Jhmean'-Jhstd)], [0.8 0.8 0.8])
% plot(n_c_vec,Jhmean)
% xlabel('$n_c$','Interpreter','latex')
% ylabel('$J_h$','Interpreter','latex')
% title('\rm Heat current from trajectories')
% grid on
% %grid on
% xscale log
% xlim([n_c_vec(1,1) n_c_vec(1,end)])
% fontsize(20,"points")
% set(gca,'linewidth',1)
% grid on
%%%
% figure
% plot(ent_prod,'LineWidth',2)
% xlabel('$n_c$','Interpreter','latex')
% ylabel('Ent. Prod.')
% title('\rm Entropy Production')
% %grid on
% xscale log
% xlim([n_c_vec(1,1) n_c_vec(1,end)])
% fontsize(20,"points")
% set(gca,'linewidth',1)
% grid on
%Allan
% figure
% histogram(dtj_stable,bin)
%----------------
%%%Put on the same figure the accuracy and the entropy production
%Next section should be ran right after, so that the next graph is
%overlayed on the last one from this section.
%----------------
figure
hold on
yyaxis left
plot(n_c_vec,N,'LineWidth',2)
xlabel('$n_{\rm c}$','Interpreter','latex')
ylabel('${\cal N}$','Interpreter','latex')
%title('\rm Accuracy')
grid on
%
yyaxis right
plot(n_c_vec,ent_prod/(w_m),'LineWidth',2)%I am deviding by the mechanical frequency
xlabel('$n_{\rm c}$','Interpreter','latex')
ylabel('$\dot \Sigma /\Omega_{\rm m}$','Interpreter','latex')
%title('\dot \Sigma')
grid on
xscale log
yscale log
xlim([n_c_vec(1,1) n_c_vec(1,end)])
xlim([n_c_vec(1,1) .01])
fontsize(20,"points")
set(gca,'linewidth',1)
box on
%% Now, Take the ticks, and consider the dead-time of the detectors
% %-----------------
% %   Tick stats; Filter
% %-----------------
plot_filter=0;
det_filt=1;
imin=1;
%imax=10;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
N=zeros(1,iTmax);
Qhfmat=zeros(iTmax,imax);
for iT=1:1:iTmax
    sub_folder_name=['n_c',num2str(iT)];
    for i1=imin:1:imax
        myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','Q_h_f','n_c_vec'};
        load([sub_folder_name,'/in_cond_n_c',num2str(iT),'traj',num2str(i1)],myVars{:})
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
        dtj=[dtj,dtjump];
        muvec(1,i1)=mean(dtjump(2:end));
        varvec(1,i1)=std(dtjump(2:end))^2;
        [i1 imax;iT iTmax]
    end
    dtj=dtj(2:end);
    % res=mean(dtj)
    % res_stable=mean(dtj_stable)
    % std_=std(dtj)
    % std_stable=std(dtj_stable)
    %%%IMPORTANT: This is a test, change back to the immediate lower two lines!!!
    mu_=mean(muvec,'omitnan');
    var_=mean(varvec,'omitnan');
    N(1,iT)=mu_.^2./var_
    %%%Total HEAT including f
    Qhfmat(iT,i1)=Q_h_f;
    %%%
    % bin=200;
    % figure
    % hold on
    % histogram(dtj(2:end),bin)
    % % Create xline
    % xline([0 1 2]);
    % tname=(['$\omega_m =$',num2str(w_m),';~~~$\mu=$',num2str(mu_),',~~$\sigma^2=$',num2str(var_),',~~$N=$',num2str(N)]);
    % title(tname,'Interpreter','latex')
    % % Create xlabel
    % xlabel('$\omega_m t/\pi$','Interpreter','latex');
end
mu_
var_
N
yyaxis left
plot(n_c_vec,N,'LineWidth',2)
xlim([n_c_vec(1,1) n_c_vec(1,end)])
xscale log
yscale log
grid on
%%%%
% figure
% plot(n_c_vec,N)
% xlim([n_c_vec(1,1) n_c_vec(1,end)])
% xscale log
% fontsize(20,"points")
% set(gca,'linewidth',1)
% grid on
%%%
% %%%HEAT CURRENT
% Qhfmean=mean(Qhfmat,2,'omitnan');
% Qhfstd=std(Qhfmat','omitnan');
% T_c_vec=w_cold./(log((1+n_c_vec)./n_c_vec));
% ent_prod=w_hot*Qhfmean'.*(1./T_h-1/T_c_vec);
% figure
% plot(n_c_vec,Qhfmean)
% figure
% plot(n_c_vec,ent_prod)
% %Allan
% % figure
% % histogram(dtj_stable,bin)
%%
%%
%%
%%
%%%
% what if we look at nu N?
%%%
%Run this code section by section if you know what you are doing. That
%would allow you to control, chose, what you are doing. Specially, the
%first section below is the most time consuming, generates data, and saves it into your computer. Later
%sections are analysing this data. You can come back to your saved data
%anytime, just make sure you don't get any errors when saving. A priori,
%you should not get errors, but you may need to change file directories
%etc.
%% Simulate trajectories.
% This simulates trajectories (imax) for each cold bath occupation (n_c_vec)
% It first runs the code for a longer time, without unraveling. This will
% be used for (1) the steady state and limit cycle properties, heat
% currents etc. (2) the start point (initial conditions) of the unravelled trajectories. 
% % If you have done it already, comment for further analysis
imin=1;
imax=100;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
%iTmax=12;n_c_vec=[logspace(-5,0,iTmax/2),linspace(1,n_h,iTmax/2)];
iTmax=20;
n_c_vec=[logspace(-5,-2,iTmax)];
n_c_vec=unique(n_c_vec);
iTmax=length(n_c_vec);
figure
for iT=1:1
    n_c=n_c_vec(1,iT);
    sub_folder_name=['n_c',num2str(iT)];
    mkdir(sub_folder_name)
    for ur=1:1
        if ur==0
            Factorisation;
            myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
                'x_m_vec','p_m_vec','J_h','J_m','J_cold','J_cav','w_hot','w_cold','w_cav','n_c_vec'};
            save([sub_folder_name,'/in_cond_n_c',num2str(iT)],myVars{:});
            subplot(9,5,iT);
            plot(x_m_vec,1i*p_m_vec,'-*','LineWidth',1);
            xlim([-40 40])
            ylim([-50 50])
            title('nc=',num2str(n_c));
            [iT iTmax]
        else
            for i1=imin:imax
                [i1,imax; iT, iTmax]
                Factorisation;
                tvec_dN1=jump_times;
                myVars2={"tvec_dN1","w_m",'w_hot','w_cold','w_cav',"n_h", 'Q_h','Q_h_f',...
                    'J_h','J_m','J_cold','J_cav','n_c_vec'};
                save([sub_folder_name,'/in_cond_n_c',num2str(iT),'traj',num2str(i1)],myVars2{:});
            end
        end
    end
end
%% Load and Plot the asymptotic limit cycles
% This plots limit cycles of the mechanical resonator.
figure
myVars0={'n_c_vec'};
sub_folder_name='n_c1';
load([sub_folder_name,'/in_cond_n_c1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
iTmax=length(n_c_vec);
icntr=1;
for iT=1:1:iTmax
    n_c=n_c_vec(1,iT);
    sub_folder_name=['n_c',num2str(iT)];
    myVars = {'x_m_vec','p_m_vec'};
    load([sub_folder_name,'/in_cond_n_c',num2str(iT)],myVars{:})
    subplot(5,4,iT);icntr=icntr+1;
    %plot(x_m_vec,1i*p_m_vec,'-*','LineWidth',1);
    plot(x_m_vec,1i*p_m_vec,'linewidth',2);
    xlim([-35 35])
    ylim([-35 35])
    title(['$\bar n_c=$',num2str(n_c)],'Interpreter','latex');
    fontsize(20,"points")
    set(gca,'linewidth',1)
    %ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
        %'Interpreter','latex','FontSize', 10);
    %xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 10);
    [iT iTmax]
end
figure
plot(x_m_vec)
%% Mechanical vs hot bath heat currents
% Title is self explanatory. It specifically can be checked, for the
% parameters in the paper, that J_m<<J_hot
figure
imin=1;
imax=100;
myVars0={'n_c_vec'};
myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
               'x_m_vec','p_m_vec','J_h','J_m','J_cold','J_cav','w_hot','w_cold','w_cav','n_c_vec'};
sub_folder_name='n_c1';
load([sub_folder_name,'/in_cond_n_c1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
sTmax=length(n_c_vec);
J_hot_vec=zeros(1,iTmax);
J_m_vec=zeros(1,iTmax);
J_cold_vec=zeros(1,iTmax);
J_cav_vec=zeros(1,iTmax);
for iT=1:iTmax
    n_c=n_c_vec(1,iT);
    sub_folder_name=['n_c',num2str(iT)];
    myVars = {'J_h','J_m','J_cold','J_cav','w_hot','w_cold','w_cav'};
    load([sub_folder_name,'/in_cond_n_c',num2str(iT)],myVars{:})
    J_hot_vec(1,iT)=J_h;
    J_m_vec(1,iT)=J_m;
    J_cold_vec(1,iT)=J_cold;
    J_cav_vec(1,iT)=J_cav;
end
%plot(n_c_vec,[-J_hot_vec/w_hot;J_cold_vec/w_cold;J_cav_vec/w_cav])
figure
plot(n_c_vec,[-J_hot_vec;J_m_vec])
fontsize(10,"points")
set(gca,'linewidth',1)
box on
xlabel('$n_{\rm c}$','Interpreter','latex')
xscale log
legend('$J_{\rm h}$','$J_{\rm m}$','interpreter','latex')
%% Load and analyze the data
% %-----------------
% %   Tick stats; No filter
% %-----------------
% We look at the tick stats, heat currents, entropy production etc.
det_filt=0;
imin=1;
imax=100;
iTmax=length(n_c_vec);
N=zeros(1,iTmax);
mu_=zeros(1,iTmax);
Var_=zeros(1,iTmax);
Jhmat=zeros(iTmax,imax);
Jmmat=zeros(iTmax,imax);
Jcoldmat=zeros(iTmax,imax);
Jcavmat=zeros(iTmax,imax);
click_num=zeros(iTmax,imax);
figure
scount=0;%counts sub-plot number for histograms.
for iT=1:iTmax
    sub_folder_name=['n_c',num2str(iT)];
    dtj=[];
    muvec=zeros(1,imax);
    varvec=zeros(1,imax);
    %for i1=imin:1:100
    for i1=imin:1:imax
        myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','J_h','J_m','J_cold','J_cav','n_c_vec','Q_h','Q_h_f'};
        load([sub_folder_name,'/in_cond_n_c',num2str(iT),'traj',num2str(i1)],myVars{:})
        %%%%This line will be passed only if you want to filter (detector dead time)
        if det_filt==1
            Detector_Filter_saturation;
            tvec_dN1=tvec_dN1_I2(1:end);
        end
        %Let's renormalise everything!
        tvec_dN1=tvec_dN1*w_m/pi;
        %%%%Otherwise carryout as usual
        dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
        size(dtjump);
        dtj=[dtj,dtjump(2:end)];
        muvec(1,i1)=mean(dtjump);
        varvec(1,i1)=std(dtjump)^2;
        %%%HEAT CURRENT
        Jhmat(iT,i1)=J_h;
        Jmmat(iT,i1)=J_m;
        Jcoldmat(iT,i1)=J_cold;
        Jcavmat(iT,i1)=J_cav;
        click_num(iT,i1)=length(tvec_dN1);
        %%%Other currents
    end
    mu_(1,iT)=mean(muvec)
    var_(1,iT)=mean(varvec)%Note we take mean of the var over different rounds.
    N(1,iT)=mu_(1,iT).^2./var_(1,iT);
    % % %%
    bin=40;
    % hold on
    scount=scount+1;
    subplot(4,5,scount)
    histogram(dtj(2:end),bin,'LineStyle','none')
    % Create xline
    xline([0 1 2]);
    tname=(['$\omega_m =$',num2str(w_m),'$n_c = $',num2str(n_c_vec(1,iT)), ...
        '$~~~\mu=$',num2str(mu_),'~~$\sigma^2=$',num2str(var_),'~~$N=$',num2str(N) ]);
    %title(tname,'Interpreter','latex')
    % Create xlabel
    %xlabel('$\omega_m t/\pi$','Interpreter','latex');
    n_c=n_c_vec(1,iT);
    title(['$\bar n_c=$',num2str(n_c)],'Interpreter','latex');
    fontsize(10,"points")
    set(gca,'linewidth',1)
    box on
end
T_c_vec=w_cold./(log((1+n_c_vec)./n_c_vec));
%ent_prod=Jhmean'.*(1./T_h-1./T_c_vec);
ent_prod=J_hot_vec.*(1./T_h-1./T_c_vec);
figure
hold on
yyaxis left
plot(n_c_vec,N./mu_,'LineWidth',2)
xlabel('$n_{\rm c}$','Interpreter','latex')
ylabel('${\cal N}\nu$','Interpreter','latex')
%title('\rm Accuracy')
grid on
%
yyaxis right
plot(n_c_vec,ent_prod/(w_m),'LineWidth',2)%I am deviding by total number of cycles (tmax*w_m/pi)
xlabel('$n_{\rm c}$','Interpreter','latex')
ylabel('$\dot \Sigma /\Omega_{\rm m}$','Interpreter','latex')
%title('\dot \Sigma')
grid on
xscale log
yscale log
xlim([n_c_vec(1,1) n_c_vec(1,end)])
xlim([n_c_vec(1,1) .01])
fontsize(20,"points")
set(gca,'linewidth',1)
box on
%% Now, Take the ticks, and consider the dead-time of the detectors
% %-----------------
% %   Tick stats; Filter
% %-----------------
det_filt=1;
imin=1;
%imax=10;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
N=zeros(1,iTmax);
Qhfmat=zeros(iTmax,imax);
for iT=1:1:20
    sub_folder_name=['n_c',num2str(iT)];
    for i1=imin:1:imax
        myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','Q_h_f','n_c_vec'};
        load([sub_folder_name,'/in_cond_n_c',num2str(iT),'traj',num2str(i1)],myVars{:})
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
        dtj=[dtj,dtjump];
        muvec(1,i1)=mean(dtjump(2:end));
        varvec(1,i1)=std(dtjump(2:end))^2;
        [i1 imax;iT iTmax]
    end
    dtj=dtj(2:end);
    % res=mean(dtj)
    % res_stable=mean(dtj_stable)
    % std_=std(dtj)
    % std_stable=std(dtj_stable)
    %%%IMPORTANT: This is a test, change back to the immediate lower two lines!!!
    mu_=mean(muvec,'omitnan');
    var_=mean(varvec,'omitnan');
    N(1,iT)=mu_.^2./var_
    %%%Total HEAT including f
    Qhfmat(iT,i1)=Q_h_f;
    %%%
    % bin=200;
    % figure
    % hold on
    % histogram(dtj(2:end),bin)
    % % Create xline
    % xline([0 1 2]);
    % tname=(['$\omega_m =$',num2str(w_m),';~~~$\mu=$',num2str(mu_),',~~$\sigma^2=$',num2str(var_),',~~$N=$',num2str(N)]);
    % title(tname,'Interpreter','latex')
    % % Create xlabel
    % xlabel('$\omega_m t/\pi$','Interpreter','latex');
end
mu_
var_
N
yyaxis left
plot(n_c_vec,N./mu_,'LineWidth',2)
xlim([n_c_vec(1,1) n_c_vec(1,end)])
xscale log
yscale log
grid on