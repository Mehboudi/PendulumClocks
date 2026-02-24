n=1000;%Downsample to 1/nth size
t_x_m_down=downsample(t_x_m,n);
t_p_m_down=downsample(t_p_m,n);
t_dN_down=downsample(t_dN,n);
t_p1_down=downsample(t_p1,n);
t_p2_down=downsample(t_p2,n);
t_p3_down=downsample(t_p3,n);
t_na_down=downsample(t_na,n);
t_re_ad_s12_down=downsample(t_re_ad_s12,n);
t_im_ad_s12_down=downsample(t_im_ad_s12,n);
tvec_down=downsample(tvec,n);
%
%-----
%The jumps cannot be downsamples, since one misses them easily!
dN1 = find(t_dN ~= 0);
tvec_dN1= tvec(dN1);
% %% 
% %Some of the plots
% figure
% plot(tvec_down,[t_p1_down;t_p2_down;t_p3_down])
% hold on
% xline(tvec_dN1,'-k','FontSize',2,'HandleVisibility','off')
% %----
% figure
% hold on
% plot(tvec_down,[t_x_m_down;1i*t_p_m_down])
% xline(tvec_dN1,'-k','FontSize',2,'HandleVisibility','off')
% %----
% figure
% plot(t_x_m_down,1i*t_p_m_down)
% %-----------------
% %   Tick stats
% %-----------------
% dtjump=[diff([0,tvec_dN1])];
% res=mean(dtjump)
% std_=std(dtjump)
% histogram(dtjump,30)
% %Secondly, if this was to be simulated by a poisson process, with rate
% %<p_3> (time average), how right will we be?
% bar_p3=mean(t_p3_down)
% res_theory=1/(g_c*bar_p3)
%% 
%   Clear old data before saving to reduce file size
clear t_dN t_na t_p_m t_x_m tvec t_im_ad_s12 t_re_ad_s12 ...
     t_p1 t_p2 t_p3
fname=(['new_data_down_resloved_',num2str(i1)]);
save(fname)