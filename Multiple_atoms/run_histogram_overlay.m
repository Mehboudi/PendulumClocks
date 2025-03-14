%% Load and analyze the data
% %-----------------
% %   Tick stats; No filter
% %-----------------
iMmax=6;
M_vec=1:iMmax;
plot_filter=1;
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
figure(1)
scount=0;%counts sub-plot number for histograms.
%You only plot the following indices (i.e., number of copies). You can modify them. 
plotindices=[1,2,4,6];splot=length(plotindices);
for iM=1:1:iMmax
    if ismember(iM, plotindices)
            scount=scount+1;
    end
    for det_filt=0:1
        sub_folder_name=['M',num2str(iM)];
        dtj=[];
        muvec=zeros(1,imax);
        varvec=zeros(1,imax);
        for i1=imin:1:imax
            myVars = {"tvec_dN1",'w_m','n_c'};
            load([sub_folder_name,'/in_cond_M',num2str(iM),'traj',num2str(i1)],myVars{:})
            %%%%This line will be passed only if you want to filter (detector dead time)
            if det_filt==1
                Detector_Filter_Saturation;
                %Detector_Filter_old;
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
        end
        muvec=muvec(1,imin:imax);
        varvec=varvec(1,imin:imax);
        mu_(1,iM)=mean(muvec);
        var_(1,iM)=mean(varvec);%Note we take mean of the var over different rounds.
        N(1,iM)=mu_(1,iM).^2./var_(1,iM);
        % % %%
        if ismember(iM, plotindices)
            bin=200;
            figure(1)
            subplot(1,splot,scount)
            hold on
            if det_filt==0
                histogram(dtj(2:end),bin,'facealpha',.5,'edgecolor','none','Normalization', 'probability')
                % Switch the axes by rotating the view
                xline([0 1 2]);
                %tname=(['$\omega_m =$',num2str(w_m),'$M = $',num2str(M_vec(1,iM)), ...
                %   '$~~~\mu=$',num2str(mu_),'~~$\sigma^2=$',num2str(var_),'~~$N=$',num2str(N) ]);
                tname=(['$M = $',num2str(M_vec(1,iM))]);
                %title(tname,'Interpreter','latex')
            else
                %histogram(dtj(2:end),bin,'facealpha',.0,'DisplayStyle','stairs','LineWidth',2)
                histogram(dtj(2:end),bin,'facealpha',.5,'edgecolor','none',...
                    'Normalization', 'probability')
                xline([0 1 2]);
                %
                %view([90 -90])
                tname=(['$M = $',num2str(M_vec(1,iM))]);
                %title(tname,'Interpreter','latex');
                if iM<7
                    % Create xlabel
                    xlabel('$\Omega_{\rm m} \tau/\pi$','Interpreter','latex');
                end
                % Add a box around the plot
                box on;

                % Set the thickness of the box to 1
                ax = gca; % Get the current axes
                ax.LineWidth = 1; % Set the line width of the box
            end
            ylim([0,.2]);
            %yticks([0,2000]);
            xticks([0,1,2]);
            figure(1)
            [iM det_filt]
            xlim([0 2]);
            %ylim([0 100])
            fontsize(20,"points")
            box on
            set(gca,'linewidth',1)
            %set(gca, 'YScale', 'log');
        end
    end
end
%
figure(10)
N
plot(M_vec,N)
