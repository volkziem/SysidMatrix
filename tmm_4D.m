% tmm_4D.m, V. Ziemann, 250310
clear all; close all;

Niter=10000;  % Number of iterations
sigm=1e-2;    % BPM noise
alpha=1;      % forgetting factor
dither=1;     % rms amplitude of the orbit variation
%rng(1234);   % start with a definite seed

L1=2; L2=1; L3=1; L4=3;  % Length of drift spaces
DD=@(L)[1,L;0,1];        % 2x2 drift matrix
D1=[DD(L1),zeros(2);zeros(2),DD(L1)];  % 4x4 matrices
D2=[DD(L2),zeros(2);zeros(2),DD(L2)];
D3=[DD(L3),zeros(2);zeros(2),DD(L3)];
D4=[DD(L4),zeros(2);zeros(2),DD(L4)];

DUT=randn(4);  % Random matrix to be recovered
tmp=det(DUT);
if tmp < 0
  DUT=[DUT(:,2),DUT(:,1),DUT(:,3),DUT(:,4)];
  DUT=DUT/sqrt(det(DUT));
else
  DUT=DUT/sqrt(det(DUT));
end
DUT
% or use any other matrix
ndim=4;                % dimension of phase space
nfit=ndim*ndim;        % number of fitted matrix elements
z0=zeros(1,ndim);
xx=zeros(ndim,1); zz=zeros(ndim,1);
P=eye(nfit);           % initial value of Pt
qhat=zeros(nfit,1);    % initial parameter extimate
data=zeros(Niter,4);   % storage for later plotting

for iter=1:Niter            % loop over iterations
  x0=dither*randn(ndim,1);  % different starting orbits, dither 
  % calculate the trajectory at BPM
  x1=x0;
  x2=D1*x1;
  x=D2*x2;
  y=DUT*x;   % the "unknown" matrix
  x3=D3*y;
  x4=D4*x3;
% add noise to BPM measurements
  x1=x1+sigm*randn(ndim,1);
  x2=x2+sigm*randn(ndim,1);
  x3=x3+sigm*randn(ndim,1);
  x4=x4+sigm*randn(ndim,1);
% reconstruct x from x1 and x2
  xx(1)=x2(1)+L2*(x2(1)-x1(1))/L1;
  xx(2)=(x2(1)-x1(1))/L1;
  xx(3)=x2(3)+L2*(x2(3)-x1(3))/L1;
  xx(4)=(x2(3)-x1(3))/L1;
  %comp=[x,xx]  % just verifying
  % reconstruct y from x3 and x4
  zz(1)=x3(1)-L3*(x4(1)-x3(1))/L4;
  zz(2)=(x4(1)-x3(1))/L4;
  zz(3)=x3(3)-L3*(x4(3)-x3(3))/L4;
  zz(4)=(x4(3)-x3(3))/L4;
  %comp=[y,zz]  % just verifying
  %....................System Identification                                      
  G=[xx',z0,z0,z0;z0,xx',z0,z0;z0,z0,xx',z0;z0,z0,z0,xx'];
  tmp2=eye(nfit)-P*G'*inv(alpha*eye(ndim)+G*P*G')*G;    
  Pnew=tmp2*P/alpha;                              
  qhat=tmp2*(qhat+P*G'*zz/alpha);
  P=Pnew;
  % display
  dDUT=DUT-[qhat(1:4)';qhat(5:8)';qhat(9:12)';qhat(13:16)'];
  data(iter,1)=trace(P'*P)/nfit^2;
  data(iter,2)=trace(dDUT'*dDUT)/ndim^2;
  data(iter,3)=P(1,1);
end
DUT_T=[qhat(1:4)';qhat(5:8)';qhat(9:12)';qhat(13:16)']
mm=1:Niter;
figure("Name","|P_T| and estimation error |dR_T|")
loglog(mm,sqrt(data(:,1)),'k',mm,sqrt(data(:,2)),'r','Linewidth',2)
xlabel('Iterations'); ylabel('|P_T|, |dR_T|')
legend('|P_T|','|dR_T|')
set(gca,'FontSize',14)
%figure("Name","Estimation error and error bar")
%loglog(mm,sqrt(data(:,2)),'k',mm,sqrt(data(:,3)),'r','Linewidth',2)