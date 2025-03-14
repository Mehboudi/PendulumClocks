% xt=tvec_dN1;
% t_max=max(xt);
% %%
% %min number of devisions, 2;
% tauvec=0:floor(log2(t_max/2));
% %tauvec=linspace(-5,log2(10),30);
% tauvec=2.^floor(tauvec);
% stau_=length(tauvec);
% allvar=zeros(1,stau_);
% for itau=1:stau_
%     tau_=tauvec(1,itau);
%     nmax=floor(t_max/tau_);
%     yvec=zeros(1,nmax-1);
%     for n=1:nmax-1
%         xnm1=sum(xt<=(n-1)*tau_);
%         xn=sum(xt<=n*tau_);
%         xnp1=sum(xt<=(n+1)*tau_);
%         yvec(1,n)=(xnp1-2*xn+xnm1)^2;
%     end
%     allvar(1,itau)=mean(yvec)/(2*tau_^2);
% end
%%
% The overlapping allanvariance
xxt=tvec_dN1';
L = length(xxt);
t0=1;
xt=zeros(1,L);
for ik=1:L
    xt(1,ik)=sum(xxt<=ik);
end
maxNumM = 100;
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.

tau = m*t0;

avard = zeros(numel(m), 1);
for i_m = 1:numel(m)
    mi = m(i_m);
    avard(i_m,:) = sum( ...
        (xt(1:L-2*mi) - 2*xt(1+mi:L-mi) + xt(1+2*mi:L)).^2, 1);
end
allvar = avard ./ (2*tau.^2 .* (L - 2*m));allvar=transpose(allvar);
tauvec=tau;

%% 
% %% Now, Take the ticks, and consider the dead-time of the detectors
% % %-----------------
% % %   Tick stats; Filter
% % %-----------------
% det_filt=1;
% imin=1;
% %imax=10;
% dtj=[];
% muvec=zeros(1,imax);
% varvec=zeros(1,imax);
% sub_folder_name='Data';
% for i1=imin:1:imax
%     myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav'};
%     load([sub_folder_name,'/conditional_traj',num2str(i1)],myVars{:})
%     %%%%This line will be passed only if you want to filter (detector dead time)
%     %%%The detector parameters are set on the other code, check it out
%     %%%there.
%     if det_filt==1
%         Detector_Filter;
%         tvec_dN1=tvec_dN1_I2(1:end);
%     end
%     %Let's renormalise everything!
%     tvec_dN1=tvec_dN1*w_m/pi;
%     %%%%Otherwise carryout as usual
%     dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
%     dtj=[dtj,dtjump];
%     muvec(1,i1)=mean(dtjump(2:end));
%     varvec(1,i1)=std(dtjump(2:end))^2;
%     [i1 imax]
% end
% dtj=dtj(2:end);
% % res=mean(dtj)
% % res_stable=mean(dtj_stable)
% % std_=std(dtj)
% % std_stable=std(dtj_stable)
% %%%IMPORTANT: This is a test, change back to the immediate lower two lines!!!
% mu_=mean(muvec,'omitnan')
% var_=mean(varvec,'omitnan')
% N=mu_.^2./var_
% bin=200;
% histogram(dtj(2:end),bin)
% % Create xline
% xline([0 1 2 3 4],'-k','FontSize',2,'HandleVisibility','off');
% tname=(['$\omega_m =$',num2str(w_m), ...
%     '$~~~\mu=$',num2str(mu_),'~~$\sigma^2=$',num2str(var_),'~~$N=$',num2str(N) ]);
% title(tname,'Interpreter','latex','FontSize',18)
% fontsize(20,"points")
% set(gca,'linewidth',1)
% xlabel('Time between consecutive ticks (\pi/\omega_m)','FontSize', 20);
% %saveas(gcf,[pwd '/Data/Pics/Histogram_filter.png'])
% %saveas(gcf,[pwd '/Data/Pics/Histogram_filter.fig'])
% Allan
% xlim([1,1e3])
% %saveas(gcf,[pwd '/Data/Pics/Allan_filter_scrach.png'])
% %saveas(gcf,[pwd '/Data/Pics//Allan_filter_scrach.fig'])
