% tmm0.m, V. Ziemann, 250228
clear all; close all;

Niter=10000;
sigm=1e-2;
alpha=1;
dither=1;
%rng(1234);   % start with a definite seed

L1=2; L2=1; L3=1; L4=3;
DD=@(L)[1,L;0,1];
D1=DD(L1);
D2=DD(L2);
D3=DD(L3);
D4=DD(L4);

DUT=randn(2);
tmp=det(DUT);
if tmp < 0
  DUT=[DUT(:,2),DUT(:,1)];
  DUT=DUT/sqrt(det(DUT));
else
  DUT=DUT/sqrt(det(DUT));
end
DUT
% or use any other matrix

xx=zeros(2,1); zz=zeros(2,1);
P=eye(4);             % initial value of Pt
qhat=zeros(4,1);      % initial parameter extimate
data=zeros(Niter,4);  % storage for later plotting

for iter=1:Niter
  x0=dither*randn(2,1); % different starting orbits, dither 
  x1=x0;
  x2=D1*x1;
  x=D2*x2;
  y=DUT*x;
  x3=D3*y;
  x4=D4*x3;
% add noise to BPM measurements
  x1=x1+sigm*randn(2,1);
  x2=x2+sigm*randn(2,1);
  x3=x3+sigm*randn(2,1);
  x4=x4+sigm*randn(2,1);
% reconstruct x from x1 and x2
  xx(1)=x2(1)+L2*(x2(1)-x1(1))/L1;
  xx(2)=(x2(1)-x1(1))/L1;
  %comp=[x,xx]  % just verifying
% reconstruct y from x3 and x4
  zz(1)=x3(1)-L3*(x4(1)-x3(1))/L4;
  zz(2)=(x4(1)-x3(1))/L4;
  %comp=[y,yy]  % just verifying
  %....................System Identification                                      
  G=[xx(1),xx(2),0,0;0,0,xx(1),xx(2)];
  tmp2=eye(4)-P*G'*inv(alpha*eye(2)+G*P*G')*G;    
  Pnew=tmp2*P/alpha;                              
  qhat=tmp2*(qhat+P*G'*zz/alpha);
  P=Pnew;
  % display
  dDUT=DUT-[qhat(1), qhat(2);qhat(3),qhat(4)];
  data(iter,1)=trace(P'*P)/16;
  data(iter,2)=trace(dDUT'*dDUT)/4;
  data(iter,3)=P(1,1);
end
DUT_T=[qhat(1), qhat(2);qhat(3),qhat(4)]
mm=1:Niter;
figure("Name","|P_T| and estimation error |dR_T|")
loglog(mm,sqrt(data(:,1)),'k',mm,sqrt(data(:,2)),'r','Linewidth',2)
xlabel('Iterations'); ylabel('|P_T|, |dR_T|')
legend('|P_T|','|dR_T|')
set(gca,'FontSize',14)
