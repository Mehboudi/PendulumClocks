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
xt=tvec_dN1';
maxNumM = 100;
L = length(xt);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.

t0=1;
tau = m*t0;

avard = zeros(numel(m), 1);
for i_m = 1:numel(m)
    mi = m(i_m);
    avard(i_m,:) = sum( ...
        (xt(1:L-2*mi) - 2*xt(1+mi:L-mi) + xt(1+2*mi:L)).^2, 1);
end
allvar = avard ./ (2*tau.^2 .* (L - 2*m));allvar=transpose(allvar);
tauvec=tau;