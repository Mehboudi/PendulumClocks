%% Simulate trajectories.
% % If you have done it already, don't run this section (or comment for
% further analysis)
plot_filter=0;
imin=1;
imax=1;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
%iTmax=12;n_c_vec=[logspace(-5,0,iTmax/2),linspace(1,n_h,iTmax/2)];
n_c=1e-5;T_c=w_cold/(log((n_c+1)/n_c));
iMmax=6;%maximum number of copies to be simulated
M_vec=1:1:iMmax;%don't change the first and second point of this vector, it should always be 1:1.
%If you want to start from a different M, set the start of the loop below
%to your desired number.
figure
for iM=1:1:iMmax
    M=M_vec(1,iM);
    sub_folder_name=['M',num2str(iM)];
    mkdir(sub_folder_name)
    for ur=0:1
        if ur==0
            Factorisation;
            myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
                'x_m_vec','p_m_vec','J_h','J_cold','J_cav','J_m','n_c','w_hot','w_cold','w_cav','M_vec'};
            save([sub_folder_name,'/in_cond_M',num2str(iM)],myVars{:});
            subplot(2,5,iM);
            plot(x_m_vec,1i*p_m_vec,'-*','LineWidth',1,'MarkerSize',2);
            %xlim([-40 40])
            %ylim([-50 50])
            title('M=',num2str(M));
            [iM iMmax]
        else
            for i1=imin:imax
                %
                Factorisation;
                tvec_dN1=jump_times;
                myVars2={"tvec_dN1","w_m",'w_hot','w_cold','w_cav',"n_h", 'n_c' 'Q_h','Q_h_f',...
                    'x_m_vec','p_m_vec','J_h','J_m','J_cold','J_cav','M_vec'};
                save([sub_folder_name,'/in_cond_M',num2str(iM),'traj',num2str(i1)],myVars2{:});
                [i1,imax; iM, iMmax]
            end
        end
    end
end
%% Load and Plot the asymptotic limit cycles
figure
myVars0={'M_vec'};
sub_folder_name='M1';
load([sub_folder_name,'/in_cond_M1','traj1'],myVars0{:});%Just pick the n_c_vec from the first available mat file
%iMmax=length(M_vec);
J_hot_vec=0*[1:iMmax];
for iM=1:iMmax
    M=M_vec(1,iM);
    sub_folder_name=['M',num2str(iM)];
    myVars = {'x_m_vec','p_m_vec','J_h'};
    load([sub_folder_name,'/in_cond_M',num2str(iM)],myVars{:})
    J_hot_vec(1,iM)=J_h;
    subplot(2,3,iM);
    plot(x_m_vec,1i*p_m_vec,'LineWidth',2)
    fontsize(20,"points")
    set(gca,'linewidth',1)
    %ylabel('$i\left\langle b - b^{\dagger} \right\rangle$',...
    %'Interpreter','latex','FontSize', 20);
    %xlabel('$\left\langle b+b^{\dagger} \right\rangle$','Interpreter','latex','FontSize', 20);
    title(['\rm M = ',num2str(iM)],'FontSize', 20);
    xlim([-150,150])
    ylim([-150,150])
    [iM iMmax]
end
%% Load and analyze the data
% %-----------------
% %   Tick stats; No filter
% %-----------------
dtjcell=cell(1,iMmax);
det_filt=0;
imin=1;
imax=1;
N=zeros(1,iMmax);
mu_=zeros(1,iMmax);
Var_=zeros(1,iMmax);
Jhmat=zeros(iMmax,imax);
Jcoldmat=zeros(iMmax,imax);
Jcavmat=zeros(iMmax,imax);
click_num=zeros(iMmax,imax);
figure
scount=0;%counts sub-plot number for histograms.
%for iT=1:iTmax
for iM=1:1:iMmax
    sub_folder_name=['M',num2str(iM)];
    dtj=[];
    muvec=zeros(1,imax);
    varvec=zeros(1,imax);
    for i1=imin:1:imax
        myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','J_h','J_cold','J_cav','n_c','Q_h','Q_h_f'};
        load([sub_folder_name,'/in_cond_M',num2str(iM),'traj',num2str(i1)],myVars{:})
        %Let's renormalise everything!
        tvec_dN1=tvec_dN1*w_m/pi;
        %%%%This line will be passed only if you want to filter (detector dead time)
        if det_filt==1
            Detector_Filter_Saturation;
            tvec_dN1=tvec_dN1_I2(1:end);
        end
        %%%%Otherwise carryout as usual
        dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
        size(dtjump);
        dtj=[dtj,dtjump];
        muvec(1,i1)=mean(dtjump);
        varvec(1,i1)=std(dtjump)^2;
        %%%HEAT CURRENT
        Jhmat(iM,i1)=J_h;
        Jcoldmat(iM,i1)=J_cold;
        Jcavmat(iM,i1)=J_cav;
        click_num(iM,i1)=length(tvec_dN1);
        %%%Other currents
    end
    %res=mean(dtj);
    % res_stable=mean(dtj_stable)
    %std_2=std(dtj)^2
    % std_stable=std(dtj_stable)
    muvec=muvec(1,imin:imax);
    varvec=varvec(1,imin:imax);
    mu_(1,iM)=mean(muvec);
    var_(1,iM)=mean(varvec);%Note we take mean of the var over different rounds.
    N(1,iM)=mu_(1,iM).^2./var_(1,iM);
    % % %%
    bin=200;
    % hold on
    scount=scount+1;
    subplot(2,3,scount)
    histogram(dtj(2:end),bin)
    xlim([.5 1.5]);
    ylim([0,250]);
    % Create xline
    xline([0 1 2]);
    %tname=(['$\omega_m =$',num2str(w_m),'$M = $',num2str(M_vec(1,iM)), ...
    %   '$~~~\mu=$',num2str(mu_),'~~$\sigma^2=$',num2str(var_),'~~$N=$',num2str(N) ]);
    tname=(['$M = $',num2str(M_vec(1,iM))]);
    title(tname,'Interpreter','latex')
    % Create xlabel
    xlabel('$\omega_m t/\pi$','Interpreter','latex');
    dtjcell{1,iM}=dtj;
    [iM sum(dtj) length(dtj); iM sum((dtj>0.5).*dtj) sum(dtj>0.5)]
end
%%%
figure
plot(M_vec,mu_)
xlabel('$M$','Interpreter','latex')
ylabel('$\mu$','Interpreter','latex')
title('Resolution')
grid on
%xscale log
xlim([M_vec(1,1) M_vec(1,end)])
figure
plot(M_vec,N)
xlabel('$M$','Interpreter','latex')
ylabel('$N$','Interpreter','latex')
title('Accuracy')
grid on
%xscale log
xlim([M_vec(1,1) M_vec(1,end)])
%%%HEAT CURRENT
Jhmean=mean(Jhmat,2,'omitnan');%The heat current (instantaniuous)
Jcoldmean=mean(Jcoldmat,2,'omitnan');%This is collective heat instead of the instantanious one
Jcavmean=mean(Jcavmat,2,'omitnan');
Q_click=-w_hot*mean(click_num,2,'omitnan');%This is the dissipated heat, with some approximation valid for low n_c (for some reason).
figure
plot(M_vec,[Jhmean/w_hot,Jcoldmean/w_cold,Jcavmean/w_cav])
xlabel('$M$','Interpreter','latex')
ylabel('$J_\alpha/\omega_\alpha$','Interpreter','latex')
title('Heat current to frequency rate')
legend('$J_{\rm hot}/\omega_{\rm hot}$','$J_{\rm cold}/\omega_{\rm cold}$' ...
    ,'$J_{\rm cav}/\Omega_{\rm cav}$','interpreter','latex')
grid on
%xscale log
xlim([M_vec(1,1) M_vec(1,end)])
%%%%%%%%%%%%%%%
Jhstd=std(Jhmat','omitnan');
T_c_vec=w_cold./(log((1+n_c)./n_c));
%ent_prod=w_hot*Jhmean'.*(1./T_h-1./T_c);
ent_prod=J_hot_vec.*(1./T_h-1./T_c)/w_m;
figure
hold on
fill([M_vec, flip(M_vec)], [Jhmean'+Jhstd, flip(Jhmean'-Jhstd)], [0.8 0.8 0.8])
plot(M_vec,Jhmean)
xlabel('$M$','Interpreter','latex')
ylabel('$J_h$','Interpreter','latex')
title('Heat current from trajectories')
grid on
%xscale log
xlim([M_vec(1,1) M_vec(1,end)])
figure
plot(M_vec,ent_prod.*M_vec)
xlabel('$M$','Interpreter','latex')
ylabel('Ent. Prod.')
title('$\dot \Sigma /\Omega_{\rm m}$','Interpreter','latex')
grid on
%xscale log
xlim([M_vec(1,1) M_vec(1,end)])
%%%%%%%%%%%%%%%
figure
hold on
yyaxis left
plot(M_vec,N,'LineWidth',2)
xlabel('$M$','Interpreter','latex')
ylabel('${\cal N}$','Interpreter','latex')
yyaxis right
plot(M_vec,ent_prod.*M_vec,'LineWidth',2)
ylabel('${\dot \Sigma}/\Omega_{\rm m}$','Interpreter','latex')
yyaxis left
yscale log
fontsize(20,"points")
set(gca,'linewidth',1)
box on
%%%%%%%%%%%%%%%
%Allan
% figure
% histogram(dtj_stable,bin)
%% Now, Take the ticks, and consider the dead-time of the detectors
% %-----------------
% %   Tick stats; Filter
% %-----------------
plot_filter=0;
det_filt=1;
imin=1;
imax=1;
scount=0;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
N=zeros(1,iMmax);
Qhfmat=zeros(iMmax,imax);
%thresh_vec=.45:0.01:.6;sthresh=length(thresh_vec);
%N_ithresh=0*thresh_vec;
%thresh_opt=zeros(1,iMmax);
for iM=1:1:iMmax
    sub_folder_name=['M',num2str(iM)];
    % for ithresh=1:sthresh
    %     threshhold=thresh_vec(1,ithresh);
        for i1=imin:1:imax
            myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','Q_h_f','n_c'};
            load([sub_folder_name,'/in_cond_M',num2str(iM),'traj',num2str(i1)],myVars{:})
            %%%%This line will be passed only if you want to filter (detector dead time)
            %%%The detector parameters are set on the other code, check it out
            %%%there.
            if det_filt==1
                Detector_Filter_Saturation;
                tvec_dN1=tvec_dN1_I2(1:end);
            end
            %Let's renormalise everything!
            tvec_dN1=tvec_dN1*w_m/pi;
            %%%%Otherwise carryout as usual
            dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
            dtj=[dtj,dtjump];
            muvec(1,i1)=mean(dtjump(2:end));
            varvec(1,i1)=std(dtjump(2:end))^2;
            [i1 imax;iM iMmax]
        end
        dtj=dtj(2:end);
        % res=mean(dtj)
        % res_stable=mean(dtj_stable)
        % std_=std(dtj)
        % std_stable=std(dtj_stable)
        %%%IMPORTANT: This is a test, change back to the immediate lower two lines!!!
        muvec=muvec(1,imin:imax);
        varvec=varvec(1,imin:imax);
        mu_=mean(muvec,'omitnan');
        var_=mean(varvec,'omitnan');
        %N_ithresh(1,ithresh)=mu_.^2./var_;
    %end
    %dummy=thresh_vec(max(N_ithresh)==N_ithresh);
    %thresh_opt(1,iM)=dummy(1,1);
    N(1,iM)=mu_.^2./var_;%max(N_ithresh);
    %%%
    bin=200;
    %figure
    % hold on
    % scount=scount+1;
    % subplot(2,3,scount)
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
%figure
yyaxis left
plot(M_vec,N,'LineWidth',2)
yscale log
xscale log
%%%
figure
plot(M_vec,N,'linewidth',1)
xlim([M_vec(1,1) M_vec(1,end)])
fontsize(20,"points")
set(gca,'linewidth',1)
ylabel('$N$','Interpreter','Latex','FontSize', 20);
xlabel('$M$','Interpreter','Latex','FontSize', 20);
