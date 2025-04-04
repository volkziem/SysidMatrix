% tmm_2D_nonlin.m, V. Ziemann, 250317
clear all;% close all;

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

ndim=2;
nfit=10;
xx=zeros(ndim,1); zz=zeros(ndim,1);
P=eye(nfit);             % initial value of Pt
qhat=zeros(nfit,1);      % initial parameter extimate
data=zeros(Niter,4);     % storage for later plotting
offs=0.0*randn(4,1);     % BPM offset, non-zero value corrupts
sext=0.1;                % sextupolar kick

for iter=1:Niter
  x0=dither*randn(2,1); % different starting orbits, dither 
  x1=x0;
  x2=D1*x1;
  x=D2*x2;
  y=DUT*x+[0;sext*x(1).^2];  % added nonlinearity
  x3=D3*y;
  x4=D4*x3;
% add noise to BPM measurements
  x1=x1+sigm*randn(2,1)+[offs(1);0];
  x2=x2+sigm*randn(2,1)+[offs(2);0];
  x3=x3+sigm*randn(2,1)+[offs(3);0];
  x4=x4+sigm*randn(2,1)+[offs(4);0];
% reconstruct x from x1 and x2
  xx(1)=x2(1)+L2*(x2(1)-x1(1))/L1;
  xx(2)=(x2(1)-x1(1))/L1;
% reconstruct y from x3 and x4
  zz(1)=x3(1)-L3*(x4(1)-x3(1))/L4;
  zz(2)=(x4(1)-x3(1))/L4;
  %....................System Identification                                      
  G=[xx(1),xx(2),xx(1)^2,xx(1)*xx(2),xx(2)^2,zeros(1,5);
    zeros(1,5),xx(1),xx(2),xx(1)^2,xx(1)*xx(2),xx(2)^2];
  tmp2=eye(nfit)-P*G'*inv(alpha*eye(ndim)+G*P*G')*G;    
  Pnew=tmp2*P/alpha;                              
  qhat=tmp2*(qhat+P*G'*zz/alpha);
  P=Pnew;
  % display
  dDUT=DUT-[qhat(1), qhat(2);qhat(6),qhat(7)];
  data(iter,1)=trace(P'*P)/16;
  data(iter,2)=trace(dDUT'*dDUT)/4;
  data(iter,3)=P(1,1);
  data(iter,4)=qhat(8);
  data(iter,5)=abs(qhat(8)-sext);
end
DUT_T=[qhat(1), qhat(2);qhat(6),qhat(7)]
mm=1:Niter;
figure(1)
loglog(mm,sqrt(data(:,1)),'k',mm,sqrt(data(:,2)),'r', ...
  mm,data(:,5),'b','Linewidth',2)
xlabel('Iterations'); ylabel('|P_T|, |dR_T|, |dm_T|')
legend('|P_T|','|dR_T|','|dm_T|')
ylim([1e-5,2.3])
set(gca,'FontSize',14)

figure(2)
semilogx(mm,data(:,4),'k','LineWidth',2)
ylim([0,0.22])
xlabel('Iterations'); ylabel('Sextupole strength')
set(gca,'FontSize',14)

qhat=qhat'