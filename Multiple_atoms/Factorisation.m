%% important: set ur=1 if you want to unravel (conditional evolution). set ur=0 if not.
%
w_m=1;%Don't touch! If you want  adifferent value, chnage it at the end of this section.
% Else everything else will be rescaled too!
f=20*w_m;%maybe reduce f to have the jumps 3-->2 more regularly?
g=30*w_m;
%%%To have local ME valid, we need interaction terms such as f and g much
%%%smaller than the other gaps
epsilon_1=0*max(f,g);
epsilon_2=4*max(f,g);
epsilon_3=8*max(f,g);
Delta=0*max(f,g);
w_cav=w_hot-w_cold+Delta;
k=10*w_m;
g_h=1*w_m;
g_c=100*w_m;
%n_h=10;This is decided from the other code
%n_c=0;%This is decided from the other code
n_opt=n_c;
n_m=n_c;
g_m=0.01*w_m;
%%%%%%%
%%AFTER SETTING ALL INITIAL PARAMETERS, I CHANGE w_m to what I want
r_alpha=.9;
w_m=r_alpha*w_m;
%%%%%%%
dt=1e-5;
dt=dt/M;%Just in case!
Dt=dt;
%%
% The initial conditions
% We first run the code with ur=0, i.e., no unraveling, such that we reach
% a steady state. Then, we use that as initial point for other [stable] simulations.
if ur==0
    % Let's take everything to be uncorrelated
    [p1,p2,p3]=deal(1,0,0);  %The atom in the GS
    na=0;                   %<a'a> Cavity is empty
    re_ad_s12=0;          %Re<a'sigma_12>
    im_ad_s12=0;           %Im<a'sigma_12>
    na_p3=0;             %<a'a p_3>
    x_m=0;                 %<b+b'>
    p_m=0*1i*sqrt(10*w_m);   %<b-b'>, This should be always imaginary, careful!!
else
    myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m"};
    load([sub_folder_name,'/in_cond_M',num2str(iM)],myVars{:})
end
%%
%%% The jump probability to the cold bath
% This is p_k=dt*tr[L_k*rho*L_k'], with L_k = sqrt(g_c*(n_c+1))sigma_23
p_j2c=@(dt,p3) dt*g_c*(n_c+1)*p3*M;
if ur==0
    %   This will be run only once, in order to reach a steady state
    tmax=1e2/w_m*pi;
else
    %   Then this will be repeated many times to get a good statistics!
    tmax=1e4/w_m*pi;
end
stvec=floor(tmax/dt);
stlimit=floor(10*pi/(w_m*dt));%This is only a fraction of the trajectory, where we believe the limit cycle already is achieved.
skipby=1e4*M;%counter for the x_m_vec, and p_m_vec below
x_m_vec=zeros(1,floor(stvec/skipby));
p_m_vec=zeros(1,floor(stvec/skipby));
ivec=1;%counter for the above vectors
jump_times=[];
t_p1_old=p1;
t_p2_old=p2;
t_p3_old=p3;
t_na_old=na;
t_re_ad_s12_old=re_ad_s12;
t_im_ad_s12_old=im_ad_s12;
t_x_m_old=x_m;
t_p_m_old=p_m;
%%% THE THERMODYNAMICS QUANTITIES
J_h=0;%The entropy production is J_h(1/T_h-1/T_c)
Q_h=0;%This takes the bare non-interacting Hamiltonian
Q_m=0;
Q_h_f=0;%This takes full Hamiltonian
Q_cold=0;
Q_cav=0;
%%%NOW to calcuate Q_h, note its defined as dt*w_h*g_h*(p3*(1+n_h) - p1*n_h);
% (this originally comes from J_h=-Tr[H_0 L_h \rho]). Then we can
%%%integrate it over time, by adding it up. At the end of the day, we can
%%%normalise it by total time to get J_h.
% The easiest way to simulate the conditional dynamics is to discretize time and in each
% time-step choose a random number between 0 and 1. If the number is
% smaller than the probability for a jump to occur, then you evolve
% the quantities with dN=1 and dt=0. If the number is larger than the
% probability of a jump to occur, then you set dN=0 (and dt~=0). In principle, you
% can do more advanced stuff by first figuring out when the next jump
% occurs (drawing from the waiting time distribution) but this might
% be a bit more involved. If you want to look into this, I would google
% Gillespie algorithm.
%ur=1;%if unravel=ur=1 we unravel, else standard Lindbladian, no jumps resolved.
% This is now decided from the other code!!
dN=0;
for itt=2:1:stvec
    if ur==1
        if rand<p_j2c(dt,t_p3_old)
            jump_times=[jump_times,dt*itt];
            dN=1;
            Dt=0;
        else
            dN=0;
            Dt=dt;
        end
    end
    %
    t_na=Dt*( 2*M*f*t_im_ad_s12_old + k*(n_opt-t_na_old));
    t_na=t_na+t_na_old;
    %---
    %---
    t_p1 =Dt*( 2*f*t_im_ad_s12_old +g_h*(n_h+1)*t_p3_old ...
        -g_h*n_h*t_p1_old  ) - t_p1_old/M *(dN-Dt*M*g_c*(n_c+1)*t_p3_old )*ur;
    t_p1 =t_p1 + t_p1_old ;
    t_p2 =Dt*( -2*f*t_im_ad_s12_old +g_c*(n_c+1)*t_p3_old ...
        -g_c*n_c*t_p2_old  ) + (1-t_p2_old )/M*( dN-Dt*M*g_c*(n_c+1)*t_p3_old  )*ur;
    t_p2 =t_p2 + t_p2_old ;
    t_p3 =Dt*( -g_h*(n_h+1)*t_p3_old -g_c*(n_c+1)*t_p3_old +...
        g_h*n_h*t_p1_old  + g_c*n_c*t_p2_old  ) - ...
        t_p3_old/M *( dN-Dt*M*g_c*(n_c+1)*t_p3_old  )*ur;
    t_p3 =t_p3 +t_p3_old ;
    %t_p3 =1-t_p2 -t_p1 ;%Alternatively, Probabilities should sum up to one
    %--
    %--
    [t_p1 ,t_p2 ,t_p3, t_na ]=...
        pop_regulate(t_p1 ,t_p2 ,t_p3, t_na, dt);
    % violation=population_test(t_p1 ,t_p2 ,t_p3 );
    % if violation==1
    %     [p1v,p2v,p3v,nav,re_ad_s12v,im_ad_s12v,x_mv,p_mv]=deal(t_p1(1,itt-1:itt)...
    %         ,t_p2(1,itt-1:itt),t_p3(1,itt-1:itt),t_na(1,itt-1:itt)...
    %         ,t_re_ad_s12(1,itt-1:itt),t_im_ad_s12(1,itt-1:itt),t_x_m(1,itt-1:itt)...
    %         ,t_p_m(1,itt-1:itt))
    %     return
    % end
    % %--
    %--
    t_re_ad_s12 =Dt*( -Delta*t_im_ad_s12_old -.5*(k+g_h*n_h+g_c*n_c)*...
        t_re_ad_s12_old -g*t_im_ad_s12_old *t_x_m_old ) - ...
        1/M*(dN-Dt*M*g_c*(n_c+1)*t_p3_old )*ur*t_re_ad_s12_old ;
    t_re_ad_s12 =t_re_ad_s12 +t_re_ad_s12_old ;
    %--
    t_im_ad_s12 =Dt*( Delta*t_re_ad_s12_old -.5*(k+g_h*n_h+g_c*n_c)*...
        t_im_ad_s12_old +g*t_re_ad_s12_old *t_x_m_old  + ...
        f*(t_p2_old +t_na_old *(t_p2_old -t_p1_old )) ) - ...
        1/M*(dN-Dt*M*g_c*(n_c+1)*t_p3_old )*ur*t_im_ad_s12_old ;
    t_im_ad_s12 =t_im_ad_s12 +t_im_ad_s12_old ;
    %---
    %---
    t_x_m =Dt*(-1i*w_m*t_p_m_old  - g_m/2*t_x_m_old );
    t_x_m =t_x_m +t_x_m_old ;
    %--
    t_p_m =Dt*(-1i*w_m*t_x_m_old -2i*g*t_na_old -g_m/2*t_p_m_old );
    t_p_m =t_p_m +t_p_m_old ;
    %--
    %%%This will be used for the next round
    t_p1_old=t_p1;
    t_p2_old=t_p2;
    t_p3_old=t_p3;
    t_na_old=t_na;
    t_re_ad_s12_old=t_re_ad_s12;
    t_im_ad_s12_old=t_im_ad_s12;
    t_x_m_old=t_x_m;
    t_p_m_old=t_p_m;
    if ur==0
        %%%Do we actually have a limit cycle?
        if itt>stvec-stlimit
            x_m_vec(1,itt+stlimit-stvec)=t_x_m;
            p_m_vec(1,itt+stlimit-stvec)=t_p_m;
        end
        %%%un_conditional_currents
        if itt>stvec-stlimit
            Q_h=Q_h+dt*w_hot*g_h*(t_p3*(1+n_h) - t_p1*n_h);
            Q_m=Q_m+dt*1i*w_m*g*t_na*t_p_m;
            Q_cold=Q_cold+dt*w_cold*g_c*((n_c+1)*t_p3 - n_c*t_p2);
            Q_cav=Q_cav+dt*k*w_cav*(n_opt-t_na);
        end
    else
        %%%THE heat current and the total heat For the conditional case
        %%%probbaly these will have no later use, but just in case!
        Q_h=Q_h+dt*w_hot*g_h*(t_p3*(1+n_h) - t_p1*n_h);
        Q_m=Q_m+dt*1i*w_m*g*t_na*t_p_m;
        %
        %%%The other currents
        Q_cold=Q_cold+dt*w_cold*g_c*((n_c+1)*t_p3 - n_c*t_p2);
        Q_cav=Q_cav+dt*w_cav*k*(n_opt-t_na);
        if mod(itt,skipby)==1
            x_m_vec(1,ivec)=t_x_m;
            p_m_vec(1,ivec)=t_p_m;
            ivec=ivec+1;
        end
    end
end
%%% Normalise the heat, to get the heat current
if ur==0
    J_h=Q_h/(dt*stlimit);
    J_m=Q_m/(dt*stlimit);
    J_cold=Q_cold/(dt*stlimit);
    J_cav=Q_cav/(dt*stlimit);
else
    J_h=Q_h/(dt*stvec);
    J_m=Q_m/(dt*stvec);
    J_cold=Q_cold/(dt*stvec);
    J_cav=Q_cav/(dt*stvec);
end
%%% RESET the initial conditions to the steady ones for next rounds (unravelling)
if ur==0
    p1=t_p1;
    p2=t_p2;
    p3=t_p3;
    na=t_na;
    re_ad_s12=t_re_ad_s12;
    im_ad_s12=t_im_ad_s12;
    x_m=t_x_m;
    p_m=t_p_m;
end
%%
function     [tp1,tp2,tp3,tna]=pop_regulate(t_p1,t_p2,t_p3,t_na,dt)
tp1=t_p1;tp2=t_p2;tp3=t_p3;tna=t_na;
if (t_p1<0) %&& (t_p1>-100*dt)
    tp1=0;
end
if (t_p2<0) %&& (t_p2>-100*dt)
    tp2=0;
end
if (t_p3<0) %&& (t_p3>-100*dt)
    tp3=0;
end
if (t_na<0) %&& (t_na>-100*dt)
    tna=0;
end
s=sum([tp1,tp2,tp3]);
[tp1,tp2,tp3]=deal(tp1/s,tp2/s,tp3/s);
end