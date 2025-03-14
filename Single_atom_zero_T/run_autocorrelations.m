%% Load and analyze the data: Tick states
% %-----------------
% %   Tick stats; No filter
% %-----------------
imin=1;
imax=1;
det_filt=0;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
sub_folder_name='Data';
lgs=100;
acf=zeros(imax,lgs+1);
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
    dtj=[dtj,dtjump];%this is not to be used here, we don't stick ticks together, only for histogram
    muvec(1,i1)=mean(dtjump);
    varvec(1,i1)=std(dtjump)^2;
    dtjump=dtjump(2:end);
    [acf(i1,:),lags] = autocorr(dtjump,NumLags=lgs,NumSTD=0);
end
mu_=mean(muvec)
var_=mean(varvec)%Note we take mean of the var over different rounds.
N=mu_.^2./var_
% % %%
%%%%Auto_correlations
%figure
%plot(lags,mean(acf))
figure
%stem(lags, mean(acf,1), 'filled', 'LineWidth', 1.5);
%or normalise differently, multiply by var to cancel, and devide by mu^2.
%That is N actually!
stem(lags, mean(acf,1)/N, 'filled', 'LineWidth', 1.5);
% Create ylabel
ylabel('${\rm Cov}(\tau_{j+k},\tau_j)/\left\langle\tau\right\rangle^2$','Interpreter','latex');
% Create xlabel
xlabel('$k$','Interpreter','latex');
box('on');
% Set the remaining axes properties
%% Now, Take the ticks, and consider the dead-time of the detectors
% %-----------------
% %   Tick stats; Filter
% %-----------------
det_filt=1;
imin=1;
imax=1;
dtj=[];
muvec=zeros(1,imax);
varvec=zeros(1,imax);
sub_folder_name='Data';
lgs=100;
acf_filter=zeros(imax,lgs+1);
for i1=imin:1:imax
    myVars = {"tvec_dN1",'w_m','w_hot','w_cold','w_cav'};
    load([sub_folder_name,'/conditional_traj',num2str(i1)],myVars{:})
    %%%%This line will be passed only if you want to filter (detector dead time)
    %%%The detector parameters are set on the other code, check it out
    %%%there.
    if det_filt==1
        Detector_Filter_old;
        tvec_dN1=tvec_dN1_I2(1:end);
    end
    %Let's renormalise everything!
    tvec_dN1=tvec_dN1*w_m/pi;
    %%%%Otherwise carryout as usual
    dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
    dtj=[dtj,dtjump];%this is not to be used here, we don't stick ticks together, only for histogram
    muvec(1,i1)=mean(dtjump(2:end));
    varvec(1,i1)=std(dtjump(2:end))^2;
    [i1 imax]
    dtjump=dtjump(2:end);
    [acf_filter(i1,:),lags] = autocorr(dtjump,NumLags=lgs,NumSTD=0);
end
% res=mean(dtj)
% res_stable=mean(dtj_stable)
% std_=std(dtj)
% std_stable=std(dtj_stable)
%%%IMPORTANT: This is a test, change back to the immediate lower two lines!!!
mu_=mean(muvec,'omitnan')
var_=mean(varvec,'omitnan')
N=mu_.^2./var_
%%%%Auto_correlations
%figure
%plot(lags,mean(acf_filter))
%autocorr(dtj(2:end),NumLags=100,NumSTD=0)
figure
%stem(lags, mean(acf_filter,1), 'filled', 'LineWidth', 1.5);
%or normalise differently, multiply by var to cancel, and devide by mu^2.
%That is N actually!
stem(lags, mean(acf_filter,1)/N, 'filled', 'LineWidth', 1.5);
% Create ylabel
ylabel('${\rm Cov}(\tau_{j+k},\tau_j)/\left\langle\tau\right\rangle^2$','Interpreter','latex');
% Create xlabel
xlabel('$k$','Interpreter','latex');
%% 
%%%%
%Every second tic is a tic (tic-toc)
%%%%
