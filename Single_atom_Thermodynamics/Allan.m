%% Allan Variance
% Important: Allan doesn't calculate the actual variance. That is
% calculated in another code. DO NOT FORGET this and load var and mu from
% the memory; this could be fatal!
%
avar=cell(1,imax);
ttau=cell(1,imax);%There might be a tau saved in the data, that's why I changed this to ttau
for i1=1:imax
    dtj=[];
    dtj_stable=[];
    myVars = {"tvec_dN1","w_m"};
    load(['new_data_down_resloved_',num2str(i1)],myVars{:})
    %%%%This line will be passed only if you want to filter (detector dead time)
    %%%The detector parameters are set on the other code, check it out
    %%%there.
    %%%%
    if det_filt==1
        Detector_Filter;
        tvec_dN1=tvec_dN1_I2(1:end);
    end
    %Let's renormalise everything!
    tvec_dN1=tvec_dN1*w_m/pi;
    %%%%Otherwise carryout as usual
    dtjump=[diff([0,tvec_dN1])];sdtj=length(dtjump);
    dtjump_stable=dtjump(floor(sdtj/2):end);
    dtj=[dtj,dtjump];%will be used in Allan_from_scratch
    dtj_stable=[dtj_stable,dtjump_stable];
    %[avar{1,i1},ttau{1,i1}]=allanvar(dtj','octave');
    Allan_from_scratch;
    avar{1,i1}=allvar;
    ttau{1,i1}=tauvec;
end
stau=zeros(1,imax);
for i1=1:imax
        stau(1,i1)=length(avar{1,i1});
end
stau=min(stau);
avarmat=zeros(imax,stau);
for i1=1:imax
    dummy=avar{1,i1};
    avarmat(i1,:)=dummy(1:stau);
end
ttau=ttau{1,1}(1:stau);
avarav=mean(avarmat,1);
figure
hold on
% loglog(ttau*w_m/pi,sqrt(avarmat),'bl:','LineWidth',.1,'HandleVisibility','off')
% loglog(ttau*w_m/pi,sqrt(avarav),'b-','LineWidth',2)
loglog(ttau,sqrt(avarmat),'bl:','LineWidth',.1,'HandleVisibility','off')
loglog(ttau,sqrt(avarav),'b-','LineWidth',2)

xlabel('Averaging Time (\tau\omega_m/\pi)')
ylabel('\sigma')
title('Allan Deviation')
grid on
%%%%%%
%According to Mitchison's paper, the Allan Variance at long times should
%converge to avar-->mu/(NT), where N=mu^2/var. with mu and var being the
%usual mean and variance. Let's check this for our signal
%mu=mean(dtj);var=std(dtj)^2;
avar_asym=mu_./(ttau*N);
%loglog(ttau*w_m/pi,sqrt(avar_asym),'r--','LineWidth',2)
loglog(ttau,sqrt(avar_asym),'r--','LineWidth',2)
legend('$\sigma_A$','$\sigma_{\rm asym}$','interpreter','latex')
