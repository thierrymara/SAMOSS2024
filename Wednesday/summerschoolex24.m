0; % mark as script for matlab
%% some helper functions (octave needs them before, matlab after use)
function [Si,Vi]=cosi(x,y,M,gfx)
% cosine sensitivity for first order effects
try
    pkg load signal
end
[xr,index]=sort(x);
yr=y(index); % Reorder output
allcoeff=dct(yr); % Compute transformation
% Unconditional variance
V = sum(allcoeff(2:end,:).^2);
% Conditional variance with M resonating harmonics
Vi= sum(allcoeff(1+(1:M),:).^2); Si= Vi./V;

% graphics
allcoeff(M+2:end,:)=0;
yhat=idct(allcoeff);
[n,k]=size(xr);
for i=1:k
    if(k>1),subplot(round(sqrt(k)),ceil(sqrt(k)),i);end
    plot(xr(:,i),yr(:,i),'b.',xr(:,i),yhat(:,i),'r-');
    title('gfx');xlabel(['x_{' num2str(i) '}']);
end
end
%%
function d=deltamim(x,y,M)
% pdf-based mim
Cutoff=1.5; % <2, Kolmogorov-Smirnov distribution
Kernel=@(x)3/(4*sqrt(5))*max(1-(x.^2/5),0);[n,k]=size(x);
iqry=median(abs(median(y)-y));
% bandwidth estimate (rule of thumb)
stdy=min(std(y),iqry/(2*0.675)); h=stdy*((4/(3*n))^(1/5));
z=linspace(min(y),max(y),100); % quadrature points
W=Kernel(bsxfun(@minus,z,y)/h)/h; densy=mean(W); %KDE
[xr,indxx]=sort(x);
for i=1:k;   xr(indxx(:,i),i)=1:n; end % ranks
for j=1:M
   indx=((j-1)*n/M<xr) & (xr <= j*n/M);   nm(:,j)=sum(indx);
   for i=1:k
   densc=mean(W(indx(:,i),:)); % conditional density
   Sm(i,j)=trapz(z,max(densy-densc,0)); %only positive part
   end
end
Sm(Sm<Cutoff.*sqrt(1/n+1./nm))=0; d=sum(Sm.*nm,2)'/n;
end
%%
function lst=copulasi(x,y)
[n,d]=size(x); uv=empcdf([x,y]); m=100; % Bernstein polynomial order
u=linspace(1/(m+1),1,m); N=max(2*m,20); t=linspace(1/(2*N),1-1/(2*N),N)';%Grids
[F,F0]=bern(t,m);  Z=zeros(m+1,m+1);
for k=1:d
 for i=1:m,for j=1:m, Z(j+1,i+1)=sum( uv(:,k)<u(i) & uv(:,end)<u(j))/n;end;end
    L=F*Z*F'-t*t'; Q=F0*Z*F0'; R=F*Z*F0'; %c. minus indep., c.density, cond.c.
    S=R-t*ones(1,N); % cdf difference
    discrepancy(k)=4*max(max(abs(L))); rho(k)=12*mean(mean(L));
    eta(k)=12*mean( mean(R).^2)-3;  delta(k)=mean(mean(abs(Q-1)))/2;
    smirnov(k)=mean(max(abs(S)));   kuiper(k)=mean(max(S)-min(S));
    gini(k)=6*mean(mean(S.^2));     hellinger(k)=1-mean(mean(sqrt(max(Q,0))));
    % gfx
    subplot(round(sqrt(d)),ceil(sqrt(d)),k);imagesc(Q);% copula density
end
% return a list
lst=struct('discrepancy',discrepancy,'rho',rho,'eta',eta,'delta',delta,...
    'smirnov',smirnov,'kuiper',kuiper,'gini',gini,'hellinger',hellinger);
end
%
function [b,B]=bern(t,n)
%BERN Bernstein base polynomials
 O=zeros(length(t),1); b=[1-t,t];
 for k=2:n
    T=repmat(t,1,k+1);
    B=b; b=(1-T).*[b,O]+T.*[O,b];
 end
 B=n*([O,B]-[B,O]); % derivative
end
%
function [ xc ] = empcdf( x )
%EMPCDF Computes the empirical cumulative distribution function
%   C=EMPCDF(X) assigns ranks (including ties)

[n,k]=size(x);
xc=zeros(n,k);
for j=1:k
 [xs,indj]=sort(x(:,j));
 xr=(1:n)';						 % ranks
 tie_loc=[find(diff(xs)==0);n+2]; % stopper
 tie_next=diff(tie_loc);
 maxt=numel(tie_loc);
 i=1;while(i < maxt)
	run=tie_loc(i);len=1;
	while(tie_next(i)==1), i=i+1;len=len+1; end
	xr( run : run + len) = run+len/2;
	i=i+1;
 end
 xc(indj,j)=(xr)/n;
end
end
%%
%
n=1024;

% constant conditional mean
%x=sobolpoints(n,3);
x=rand(n,3);
m=5;
y1=x(:,1).*(x(:,2).^m-(1-m/(m+1)))+.1*x(:,3);
subplot(1,4,1)
cosi(x(:,1),y1,4,'Constant conditional mean')

y2=x(:,1)-2*x(:,1).*(x(:,2)>.6)+sqrt(x(:,3));

subplot(1,4,2)
cosi(x(:,1),y2,4,'Bifurcation')

subplot(1,4,3)
y3=exp(-10*exp(y2));
cosi(x(:,1),y3,4,'Orders of magnitude')

try
  pkg load statistics
end
subplot(1,4,4)
y4=((x(:,2)+x(:,3)).*betainv(x(:,1),x(:,2),x(:,3))-x(:,2))./sqrt(x(:,2).*x(:,3)./...
(x(:,2)+x(:,3)+1));

cosi(x(:,2),y4,4,'Conditional Skewness')
%%
deltamim(x,y1,8)
%%
copulasi(x,y1)

