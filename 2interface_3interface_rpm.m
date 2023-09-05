clear
%set up program for solving solid state diffusion

%This script obtains diffusion and kinetic constants for the data in the excel file named Fe2O3_Fe3O4_3CO_15CO2_2Cycle_150_300um_NQW_900C_10032022
Sheetnumber='Sheet1';
data=xlsread('Fe2O3_FeO_5CO_5CO2_1Cycle_150_300um_NQW_850C_24032022.xlsx',Sheetnumber);

Time=data(:,2);
Sdata=data(:,4);

endtime=Time(end);
numdata=length(Time);

numdatacells=421
endtime=Time(421);
numdata=43;

T=1123;
t=Time(1:10:421);

p=[1.3660e-06*exp(-4.8267e+04./(8.3145*T))*0.7, 1.7216e-09*exp(-3.1276e+04./(8.3145*T)), 4.450100406743612e-09*1,8.995334033932562e-11*0.6];
p0=[1.3660e-06*exp(-4.8267e+04./(8.3145*T))*0.7, 1.7216e-09*exp(-3.1276e+04./(8.3145*T)), 4.450100406743612e-09*1,8.995334033932562e-11*0.6];

options = optimoptions('lsqcurvefit','StepTolerance',1e-6);
ub = [1e-5,1e-10,1e-5,1e-10];
lb = [0,0,0,0];


%[p,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@kineticsfun,p0,Time(1:10:421),Sdata(1:10:421),lb,ub,options);
%ci=nlparci(p,Rsd,'jacobian',Jmat);

results=kineticsfun(p,t)
figure(15)
plot(Time(1:10:421),Sdata(1:10:421),'o',Time(1:10:421),results,'LineWidth',2)
set(gca,'FontSize',16)
l=[legend('TGA Experiment','Random Pore Model') ylabel('Conversion $X_{2}$') xlabel('Time /s')]
set(l,'Interpreter','latex')
legend('TGA Experiment','Random Pore Model')

function S=kineticsfun(p,t)

syms r1(t) r2(t) r3(t) r3dash(t) Vol1(t) Vol2(t) Vol3(t) Ci1(t) Ci2(t) Ci3(t) Ci1H2O(t) Ci2H2O(t) Ci3H2O(t) D3 ks3 V1 V2 V3 Cis CisH2O structuralparameter Z1 Z2 Z3 B1 B2 B3 porosity S01 r0 x K1 K2 K3 numdatacells

%syms r1(t) r2(t) r3(t) r3dash(t) Vol1(t) Vol2(t) Vol3(t) Ci1(t) Ci2(t) Ci3(t) D1 D2 D3 ks2 ks3 V1 V2 V3 Cis structuralparameter Z1 Z2 Z3 B1 B2 B3 porosity S01 r0 x K1 K2 K3
assume(r1(t)-r2(t),'positive');
assume(r2(t)-r3(t),'positive');
assume(r3(t),'positive');
assume(Ci1(t),'positive');
assume(Ci2(t),'positive');
assume(Ci3(t),'positive');

Sheetnumber='Sheet1';
data=xlsread('Fe2O3_FeO_5CO_5CO2_1Cycle_150_300um_NQW_850C_24032022.xlsx',Sheetnumber);

Time=data(:,2);
Sdata=data(:,4);

Time=Time(1:10:421);
Sdata=data(1:10:421);


endtime=Time(end);
numdata=length(Time);

endtime=Time(end);
numdata=43;

Vol0=porosity;

ks1=p(1);
D1=p(2);
ks2=p(3);
D2=p(4);

Vol1dash=Vol1(t)-Z1*(Vol1(t)-porosity);
Vol2dash=Vol2(t)-Z2*(Vol2(t)-Vol1dash);
Vol3dash=Vol3(t)-Z3*(Vol3(t)-Vol2dash);

r1(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol1(t))/(1-porosity))))-1)+r0;
r2(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol2(t))/(1-porosity))))-1)+r0;
r3(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol3(t))/(1-porosity))))-1)+r0;
r3dash(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity))))-1)+r0;


%{
% This set of equations for all positive chemical driving forces

eqn1=(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*(Cis-Ci3(t))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))+ks3*(Ci3(t)-Ci3H2O(t)/K3)/V3;
eqn2=(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))+((4*x-3)/x)*ks2*(Ci2(t)-Ci2H2O(t)/K2)/V2;
eqn3=(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))==(ks1*(Ci1(t)-Ci1H2O(t)/K1))/(3*V1);

eqn4=-(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*((CisH2O-Ci3H2O(t)))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3H2O(t)-Ci2H2O(t))+(ks3*(Ci3(t)-Ci3H2O(t)/K3)/V3);
eqn5=-(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*((Ci3H2O(t)-Ci2H2O(t)))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2H2O(t)-Ci1H2O(t))+((4*x-3)/x)*(ks2*(Ci2(t)-Ci2H2O(t)/K2))/V2;
eqn6=-(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*((Ci2H2O(t)-Ci1H2O(t)))==(ks1*(Ci1(t)-Ci1H2O(t)/K1)/(3*V1);
%}


% For no second reaction

%This set of equations for (Cis-CisH2O/K2)<0
%{
eqn1=(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*(Cis-Ci3(t))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))+0;
eqn2=(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))+0;
eqn3=(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))==ks1*(Ci1(t)-Ci1H2O(t)/K1)/(3*V1);

eqn4=-(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*((CisH2O-Ci3H2O(t)))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3H2O(t)-Ci2H2O(t))+0;
eqn5=-(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*((Ci3H2O(t)-Ci2H2O(t)))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2H2O(t)-Ci1H2O(t))+0;
eqn6=-(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*((Ci2H2O(t)-Ci1H2O(t)))==(ks1*(Ci1(t)-Ci1H2O(t)/K1))/(3*V1);
%}

%{

% For no third reaction
%This set of equations for (Cis-CisH2O/K3)<0

    
    eqn1=(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*(Cis-Ci3(t))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))+0;
eqn2=(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))+((4*x-3)/x)*ks2*(Ci2(t)-Ci2H2O(t)/K2)/V2;
eqn3=(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))==ks1*(Ci1(t)-Ci1H2O(t)/K1)/(3*V1);

eqn4=-(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*((CisH2O-Ci3H2O(t)))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3H2O(t)-Ci2H2O(t))+0;
eqn5=-(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*((Ci3H2O(t)-Ci2H2O(t)))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2H2O(t)-Ci1H2O(t))+((4*x-3)/x)*(ks2*(Ci2(t)-Ci2H2O(t)/K2))/V2;
eqn6=-(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*((Ci2H2O(t)-Ci1H2O(t)))==(ks1*(Ci1(t)-Ci1H2O(t)/K1))/(3*V1);

%}

 %eqn1=(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*(Cis-Ci3(t))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))+0;
eqn2=(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Cis-Ci2(t))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))+((4*x-3)/x)*ks2*(Ci2(t)-Ci2H2O(t)/K2)/V2;
eqn3=(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))==ks1*(Ci1(t)-Ci1H2O(t)/K1)/(3*V1);

%eqn4=-(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*((CisH2O-Ci3H2O(t)))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3H2O(t)-Ci2H2O(t))+0;
eqn5=-(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*((CisH2O-Ci2H2O(t)))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2H2O(t)-Ci1H2O(t))+((4*x-3)/x)*(ks2*(Ci2(t)-Ci2H2O(t)/K2))/V2;
eqn6=-(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*((Ci2H2O(t)-Ci1H2O(t)))==(ks1*(Ci1(t)-Ci1H2O(t)/K1))/(3*V1);
%}

[solCi1(t),solCi2(t),solCi1H2O(t),solCi2H2O(t)]=solve([eqn2,eqn3,eqn5,eqn6],[Ci1(t),Ci2(t),Ci1H2O(t),Ci2H2O(t)]);

Ci1(t)=solCi1(t);
Ci2(t)=solCi2(t);
%Ci3(t)=solCi3(t);
Ci1H2O(t)=solCi1H2O(t);
Ci2H2O(t)=solCi2H2O(t);
%Ci3H2O(t)=solCi3H2O(t);





S1=S01*(1-Vol1(t))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol1(t))/(1-porosity)));
S2=S01*(1-Vol2(t))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol2(t))/(1-porosity)));
S3=S01*(1-Vol3(t))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol3(t))/(1-porosity)));
S3dash=S01*(1-Vol3dash)/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity)));

%{ 
% For all three reactions
eqn7=diff(Vol1(t),t)==(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
dVol1dt=(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
eqn8=diff(Vol2(t),t)==(1-Z1)*dVol1dt+ks2*(Ci2(t)-Ci2H2O(t)/K2)*S2;
dVol2dt=(1-Z1)*dVol1dt+ks2*(Ci2(t)-Ci2H2O(t)/K2)*S2;
eqn9=diff(Vol3(t),t)==(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+ks3*(Ci3(t)-Ci3H2O(t)/K3)*S3;
dVol3dt=(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+ks3*(Ci3(t)-Ci3H2O(t)/K3)*S3;
%}


% For no second reaction

%{
eqn7=diff(Vol1(t),t)==(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
dVol1dt=(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
eqn8=diff(Vol2(t),t)==(1-Z1)*dVol1dt+0;
dVol2dt=(1-Z1)*dVol1dt+0;
eqn9=diff(Vol3(t),t)==(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
dVol3dt=(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
%}


% For no third reaction

eqn7=diff(Vol1(t),t)==(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
dVol1dt=(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
eqn8=diff(Vol2(t),t)==(1-Z1)*dVol1dt+ks2*(Ci2(t)-Ci2H2O(t)/K2)*S2;
dVol2dt=(1-Z1)*dVol1dt+ks2*(Ci2(t)-Ci2H2O(t)/K2)*S2;
eqn9=diff(Vol3(t),t)==(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
dVol3dt=(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
dVol3dash=dVol3dt*(1-Z3)+dVol2dt*(1-Z2)*Z3+dVol1dt*Z2*Z3*(1-Z1);

eqs=[eqn7,eqn8,eqn9];
vars=[Vol1(t) Vol2(t) Vol3(t)];

[M,F]=massMatrixForm(eqs,vars);

M=odeFunction(M,vars);
F=odeFunction(F,vars,D3,ks3,V1,V2,V3,Cis,CisH2O,structuralparameter,Z1,Z2,Z3,porosity,S01,r0,x,K1,K2,K3);

T=1123;


%ks3=1*(1.1*10^-6).*exp(-64000./(8.3145.*T)); %Liu et al, 2014
%D3=1*5*10^(-7)*exp(-90*10^3/(8.3145*T)); %from Liu 2014, page 160 !!check this!!!

%new data

K1=7.583567432175545e+04;
K2=1.704654782419324;
K3=0.491055794220954;

%D1=1.7216e-09*exp(-3.1276e+04./(8.3145*T))
%ks1=1.3660e-06*exp(-4.8267e+04./(8.3145*T))


Ci1(t)=subs(Ci1(t));
Ci2(t)=subs(Ci2(t));
%Ci3(t)=subs(Ci3(t));
Ci1H2O(t)=subs(Ci1H2O(t));
Ci2H2O(t)=subs(Ci2H2O(t));
%Ci3H2O(t)=subs(Ci3H2O(t));

V1=0.159692/5250; %Paul fennell paper
V2=0.231533/5100;
V3=0.068885/5600;
V4=0.055485/7870;


S01=2*4041240;
S01=1.33*10^7;
S01=7.120109100000000e+06;

x=0.947;

Z1=2.13/2.14; %sourced from paper by fennell hayhurst dennis 2011. Maybe not interpreted correctly here
Z2=1.81/2.13;
Z3=1/1.81;

Z1=V2/V1;
Z2=V3/V2;
Z3=V4/V3;


porosity=0.22
structuralparameter=4.04;

r0=sqrt(porosity/(pi()*(structuralparameter*S01^2/(4*pi()*(1-porosity)))));


Cis=(0.05*10^5)/(8.3145*T);
CisH2O=(0.05*10^5)/(8.3145*T);

Cis-CisH2O/K1;
Cis-CisH2O/K2;
Cis-CisH2O/K3;

Ci1(t)=subs(Ci1(t));
Ci2(t)=subs(Ci2(t));
%Ci3(t)=subs(Ci3(t));

Ci1H2O(t)=subs(Ci1H2O(t));
Ci2H2O(t)=subs(Ci2H2O(t));
%Ci3H2O(t)=subs(Ci3H2O(t));

F=@(t,Y) F(t,Y,D3,ks3,V1,V2,V3,Cis,CisH2O,structuralparameter,Z1,Z2,Z3,porosity,S01,r0,x,K1,K2,K3);

y0=[porosity;porosity;porosity];
opt=odeset('mass',M);
[t Y]=ode45(F,linspace(0,endtime,numdata),y0,opt);

X1=(Y(:,1)-porosity)/(1-porosity);
X2=(Y(:,2)-(1-Z1)*Y(:,1)-Z1*porosity)/(Z1*(1-porosity));
X3=(Y(:,3)-(1-Z2)*Y(:,2)-Z2*(1-Z1)*Y(:,1)-Z1*Z2*porosity)/(Z1*Z2*(1-porosity));

%The volume enclosed by the wustite layer interface is less than (1-Z2) times the
%volume enclosed by the 
%magnetite layer inteface

figure(1)
plot(t,[X1,X2,X3],'LineWidth',2)
legend('X_1','X_2','X_3')
xlabel('Time /s','FontSize',16)
ylabel('Conversion','FontSize',16)
xlim([0 200])
set(gca,'FontSize',16)
xlabel('Time /s')
ylabel('Conversion X_j')

figure(6)
plot(t,Y(:,1),t,Y(:,2),t,Y(:,3),'LineWidth',2)
legend('V_1','V_2','V_3')
xlim([0 200])
set(gca,'FontSize',16)
xlabel('Time /s')
ylabel('V_i (dimensionless)')

%calculating the derivatives

Vol1dash=Y(:,1)-Z1*(Y(:,1)-porosity);
Vol2dash=Y(:,2)-Z2*(Y(:,2)-Vol1dash);
Vol3dash=Y(:,3)-Z3*(Y(:,3)-Vol2dash);

S1=S01.*(1-Y(:,1))./(1-porosity).*sqrt(1-structuralparameter.*log((1-Y(:,1))./(1-porosity)));
S2=S01.*(1-Y(:,2))./(1-porosity).*sqrt(1-structuralparameter.*log((1-Y(:,2))./(1-porosity)));
S3=S01.*(1-Y(:,3))./(1-porosity).*sqrt(1-structuralparameter.*log((1-Y(:,3))./(1-porosity)));
S3dash=S01.*(1-Vol3dash)./(1-porosity).*sqrt(1-structuralparameter.*log((1-Vol3dash)./(1-porosity)));

figure(8)
plot(t,S1,t,S2,t,S3,'LineWidth',2)
legend('S1','S2','S3')
xlim([0 200])
set(gca,'FontSize',16)
xlabel('Time /s')
ylabel('Surface area /m^-^1')

conversion=0.1111*X1+0.1889*X2+0.7*X3;
S=real(100*(0.1111*X1+0.1889*X2)/0.3);

end