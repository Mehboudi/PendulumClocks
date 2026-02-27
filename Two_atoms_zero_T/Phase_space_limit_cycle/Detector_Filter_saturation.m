%%% This code filters the observed data to take into account
% the detectors death time. We model this as an exponential decay
% that respects causality. As a result, two shortly observed
% clicks will not be counted twices; as the detector is blind.
%dec_r=w_m*(1+log10(iM));%
%dec_r=w_m*(1+(iM-1)/(-1.6*log10(n_c)));%
dec_r=w_m;
ti=[0,tvec_dN1];
n=10;
I2=[];
tt=[];
for iti=1:length(ti)-1
    tti=linspace(ti(1,iti),ti(1,iti+1),n);
    tti(1,end)=ti(1,iti+1)-1e-3*(ti(1,iti+1)-ti(1,iti));%just to avoid repeating it again in the next iteration
    dummy=exp(-dec_r*(tti-ti(iti)));
    %if iti>iM
    for jM=1:iM-1
        if iti-jM>0
            dummy=dummy+exp(-dec_r*(tti-ti(iti-jM)));%residue from the previous tcik
        end
    end
    %end
    I2=[I2,dummy];
    tt=[tt,tti];
end
%
%%%%%%%%%%%%%
%Put a Cap (threshold) for the current, above which we consider as a click
%   This code will consider the first time I2 goes above threshold
%   a click. It will not count clicks until I2 is below threshold.
threshholdup=max(.3*sqrt(iM+1),.99);%You can change this to match the data better
% The above two values work great for iM=1:6
threshholddown=.25*sqrt(iM+1);
if threshholddown <= .5
    threshholddown=iM/(iM+3);
end
index=1;
tvec_dN1_I2=[];
for itt=1:length(tvec_dN1)
    if and(I2(1,(itt-1)*n+1)>threshholdup,index==1)
        tvec_dN1_I2=[tvec_dN1_I2,tt(1,(itt-1)*n+1)];
        index=0;
    end
    if and(I2(1,itt*n)<threshholddown,index==0)
        index=1;
    end
end
%%
%plot the currents
if plot_filter==1 && i1==1
    figure(30)
    l = 0;
    u = 5;
    x = tt*w_m/pi;%this is our x axis (normalised time)
    % Logical indexing to find elements within the range (l, u)
    indices = (x > l) & (x < u);% for the times and the current
    %
    tvec_dN1_limit=tvec_dN1*w_m/pi;
    tvec_dN1_limit=tvec_dN1_limit(tvec_dN1_limit<u);
    tvec_dN1_limit=tvec_dN1_limit(tvec_dN1_limit>l);
    %
    tvec_dN1_I2_limit=tvec_dN1_I2*w_m/pi;
    tvec_dN1_I2_limit=tvec_dN1_I2_limit(tvec_dN1_I2_limit<u);
    tvec_dN1_I2_limit=tvec_dN1_I2_limit(tvec_dN1_I2_limit>l);
    %
    plot(x(indices),I2(indices),'LineWidth',2,'Color','b')
    current_ylim = ylim;
    hold on
    xline(tvec_dN1_limit,'LineWidth',.01,'LineStyle','--','color',[0,0,0],'HandleVisibility','off')
    xline(tvec_dN1_I2_limit,'LineWidth',2,'color',[0,0,0],'HandleVisibility','off')
    xlim([l,u])
    ylim([0,1.1])
    fontsize(20,"points")
    set(gca,'linewidth',1)
    if iM <2
        ylabel('$I(t)$','FontSize', 20,'Interpreter','latex');
    end
    xlabel('$\Omega_{\rm m} t/\pi$','Interpreter','latex','FontSize', 20);
    %yline(threshholdup,'b--','LineWidth',2)
    yline(threshholddown,'r--','LineWidth',2)
    %tname=(['$M = $',num2str(M_vec(1,iM))]);
    %title(tname,'Interpreter','latex')
    % Add a box around the plot
    box on;

    % Set the thickness of the box to 1
    ax = gca; % Get the current axes
    ax.LineWidth = 1; % Set the line width of the box
    %saveas(gcf,[pwd '/Data/Pics/Det_current_ticks_filtered.png'])
    %saveas(gcf,[pwd '/Data/Pics/Det_current_ticks_filtered.fig'])
end
%%
%%%Take every second tic as a tic (tic-toc)
%tvec_dN1_I2 = tvec_dN1_I2(3:2:end);
%tvec_dN1_I2 = tic_toc(tvec_dN1_I2);
%%
function tt=over_sample(ti,n)
%Take a vector ti, and produce a number between any two
% consecutive numbers in it
% n is the oversampling size
sti=length(ti);
tt=[];
for iti=1:sti-1
    tt_dumy=linspace(ti(1,iti),ti(1,iti+1),n);
    %tt_dumy=tt_dumy(2:end);
    tt=[tt,tt_dumy];
end
tt=sort(unique(tt));
end

function tictoc = tic_toc(tic)
% Calculate the length of the output vector u
n = floor(length(tic) / 2);

% Initialize the output vector u
tictoc = zeros(1, n);

% Sum successive elements of V
for i = 1:n
    tictoc(i) = tic(2*i - 1) + tic(2*i);
end
end