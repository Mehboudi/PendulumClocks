%% Load and analyze the data
% %-----------------
% %   Tick stats; No filter
% %-----------------
imin=1;
%imax=1;
%iTmax=length(n_c_vec);%iTmax=20;
N=zeros(1,iTmax);
mu_=zeros(1,iTmax);
Var_=zeros(1,iTmax);
Jhmat=zeros(iTmax,imax);
Jcoldmat=zeros(iTmax,imax);
Jcavmat=zeros(iTmax,imax);
click_num=zeros(iTmax,imax);
figure
scount=0;%counts sub-plot number for histograms.
for iT=[1,9,14]
    scount=scount+1;
    subplot(1,3,scount)
    hold on
    for det_filt=0:1
        sub_folder_name=['n_c',num2str(iT)];
        dtj=[];
        muvec=zeros(1,imax);
        varvec=zeros(1,imax);
        for i1=imin:1:imax
            myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav','J_h','J_cold','J_cav','n_c_vec','Q_h','Q_h_f'};
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
            Jcoldmat(iT,i1)=J_cold;
            Jcavmat(iT,i1)=J_cav;
            click_num(iT,i1)=length(tvec_dN1);
            %%%Other currents
        end
        mu_(1,iT)=mean(muvec)
        var_(1,iT)=mean(varvec)%Note we take mean of the var over different rounds.
        N(1,iT)=mu_(1,iT).^2./var_(1,iT);
        % % %%
        bin=300;
        % hold on
        if det_filt==0
            histogram(dtj(2:end),bin,'facealpha',.5,'edgecolor','none',...
                Normalization='probability')
        else
            %histogram(dtj(2:end),bin,'facealpha',.0,'DisplayStyle','stairs','LineWidth',2)
            histogram(dtj(2:end),bin,'facealpha',.7,'edgecolor','none',...
                Normalization='probability')
        end
        % Create xline
        xline([0 1 2]);
        tname=(['$\omega_m =$',num2str(w_m),'$n_c = $',num2str(n_c_vec(1,iT)), ...
            '$~~~\mu=$',num2str(mu_),'~~$\sigma^2=$',num2str(var_),'~~$N=$',num2str(N) ]);
        %title(tname,'Interpreter','latex')
        % Create xlabel
        xlabel('$\Omega_m \tau/\pi$','Interpreter','latex');
        n_c=n_c_vec(1,iT);
        title(['$\bar n_c=$',num2str(n_c)],'Interpreter','latex');
        fontsize(20,"points")
        box on
        set(gca,'linewidth',1)  
        xlim([0,2.2])
        ylim([0,.06])
        %ytickformat('%.1e'); % Set y-axis to scientific notation
        %legend('no filter','filter')
    end
end