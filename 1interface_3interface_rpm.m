clear
%\set up program for solving solid state diffusion

%This script obtains diffusion and kinetic constants for the data in the excel file named Fe2O3_Fe3O4_3CO_15CO2_2Cycle_150_300um_NQW_900C_10032022
Sheetnumber='Sheet1';
data=xlsread('Fe2O3_Fe3O4_3CO_15CO2_2Cycle_150_300um_NQW_800C_10032022',Sheetnumber);

Time=data(:,2); %time data
Sdata=data(:,4); %conversion data

endtime=Time(201);
numdata=201;

T=1073;
t=Time;

%ake notes for Jose on fitting process. Start by visually estimating good
%parameters adjusting p (with lsqcurvefit and ci lines of code blocked out).
%Then intput these estimates into p0 to find optimal values of
%parameters with lsqcurve fit function.
p=[3.853594175840388e-09*1.45,5.000000000000000e-11*21]; %resultant kinetic constants ks1 in column 1 and diffusion coefficient D1 in column 2
p0=[3.853594175840388e-09*1.45,5.1112342356e-11*20]; %guesses kinetic constants ks1 in column 1 and diffusion coefficient D1 in column 2

options = optimoptions('lsqcurvefit','StepTolerance',1e-6);

ub = [1e-5,1e-10];
lb = [0,0];

[p,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@kineticsfun,p0,Time(1:201),Sdata(1:201),lb,ub,options);
ci=nlparci(p,Rsd,'jacobian',Jmat);

results=kineticsfun(p,t);

figure(15)
plot(Time(1:201),Sdata(1:201),'o',Time(1:201),results,'LineWidth',2)
set(gca,'FontSize',16)
l=[legend('TGA Experiment','Random Pore Model') ylabel('Conversion $X_{1}$') xlabel('Time /s')]
xlim([0 100])
set(l,'Interpreter','latex')
legend('TGA Experiment','Random Pore Model')


dX1dt=gradient(results(1:10:201))./gradient(Time(1:10:201))
dX1dt_data=gradient(Sdata(1:10:201))./gradient(Time(1:10:201))


figure(700)
plot(results(1:10:201),dX1dt)
hold on
plot(Sdata(1:10:201),dX1dt_data,'o')
hold off


function S=kineticsfun(p,t)


syms r1(t) r2(t) r3(t) r3dash(t) Vol1(t) Vol2(t) Vol3(t) Ci1(t) Ci2(t) Ci3(t) Ci1H2O(t) Ci2H2O(t) Ci3H2O(t) D2 D3 ks2 ks3 V1 V2 V3 Cis CisH2O structuralparameter Z1 Z2 Z3 B1 B2 B3 porosity S01 r0 x K1 K2 K3

assume(r1(t)-r2(t),'positive');
assume(r2(t)-r3(t),'positive');
assume(r3(t),'positive');
assume(Ci1(t),'positive');
assume(Ci2(t),'positive');
assume(Ci3(t),'positive');

Sheetnumber='Sheet1';
data=xlsread('Fe2O3_Fe3O4_3CO_15CO2_2Cycle_150_300um_NQW_800C_10032022',Sheetnumber);

Time=data(:,2);
Sdata=data(:,4);

endtime=Time(201);
numdata=201;

Vol0=porosity;

ks1=p(1);
D1=p(2);

Vol1dash=Vol1(t)-Z1*(Vol1(t)-porosity);
Vol2dash=Vol2(t)-Z2*(Vol2(t)-Vol1dash);
Vol3dash=Vol3(t)-Z3*(Vol3(t)-Vol2dash);

r1(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol1(t))/(1-porosity))))-1)+r0;
r2(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol2(t))/(1-porosity))))-1)+r0;
r3(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol3(t))/(1-porosity))))-1)+r0;
r3dash(t)=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity))))-1)+r0;

%{
% This set of equations for all positive chemical driving forces, i.e. full
reduction from Fe2O3 to Fe

eqn1=(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*(Cis-Ci3(t))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))+ks3*(Ci3(t)-Ci3H2O(t)/K3)/V3;
eqn2=(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))+((4*x-3)/x)*ks2*(Ci2(t)-Ci2H2O(t)/K2)/V2;
eqn3=(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))==(ks1*(Ci1(t)-Ci1H2O(t)/K1))/(3*V1);

eqn4=-(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*((CisH2O-Ci3H2O(t)))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3H2O(t)-Ci2H2O(t))+(ks3*(Ci3(t)-Ci3H2O(t)/K3)/V3);
eqn5=-(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*((Ci3H2O(t)-Ci2H2O(t)))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2H2O(t)-Ci1H2O(t))+((4*x-3)/x)*(ks2*(Ci2(t)-Ci2H2O(t)/K2))/V2;
eqn6=-(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*((Ci2H2O(t)-Ci1H2O(t)))==(ks1*(Ci1(t)-Ci1H2O(t)/K1)/(3*V1);
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


% For no second or third reaction

%This set of equations for (Cis-CisH2O/K2)<0

eqn1=(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*(Cis-Ci3(t))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))+0;
eqn2=(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3(t)-Ci2(t))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))+0;
eqn3=(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2(t)-Ci1(t))==ks1*(Ci1(t)-Ci1H2O(t)/K1)/(3*V1);

eqn4=-(D3/((abs(r3(t)-r3dash(t))+r3(t)-r3dash(t))/2))*((CisH2O-Ci3H2O(t)))==(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*(Ci3H2O(t)-Ci2H2O(t))+0;
eqn5=-(D2/((abs(r2(t)-r3(t))+r2(t)-r3(t))/2))*((Ci3H2O(t)-Ci2H2O(t)))==(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*(Ci2H2O(t)-Ci1H2O(t))+0;
eqn6=-(D1/((abs(r1(t)-r2(t))+r1(t)-r2(t))/2))*((Ci2H2O(t)-Ci1H2O(t)))==(ks1*(Ci1(t)-Ci1H2O(t)/K1))/(3*V1);

[solCi1(t),solCi2(t),solCi3(t),solCi1H2O(t),solCi2H2O(t),solCi3H2O(t)]=solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[Ci1(t),Ci2(t),Ci3(t),Ci1H2O(t),Ci2H2O(t),Ci3H2O(t)]);

Ci1(t)=solCi1(t);
Ci2(t)=solCi2(t);
Ci3(t)=solCi3(t);
Ci1H2O(t)=solCi1H2O(t);
Ci2H2O(t)=solCi2H2O(t);
Ci3H2O(t)=solCi3H2O(t);

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

% For no second or third reaction

eqn7=diff(Vol1(t),t)==(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
dVol1dt=(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
eqn8=diff(Vol2(t),t)==(1-Z1)*dVol1dt+0;
dVol2dt=(1-Z1)*dVol1dt+0;
eqn9=diff(Vol3(t),t)==(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
dVol3dt=(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;

%{ 
% For no third reaction
eqn7=diff(Vol1(t),t)==(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
dVol1dt=(ks1*(Ci1(t)-Ci1H2O(t)/K1))*S1;
eqn8=diff(Vol2(t),t)==(1-Z1)*dVol1dt+ks2*(Ci2(t)-Ci2H2O(t)/K2)*S2;
dVol2dt=(1-Z1)*dVol1dt+ks2*(Ci2(t)-Ci2H2O(t)/K2)*S2;
eqn9=diff(Vol3(t),t)==(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
dVol3dt=(1-Z2)*dVol2dt+Z2*(1-Z1)*dVol1dt+0;
%}

dVol3dash=dVol3dt*(1-Z3)+dVol2dt*(1-Z2)*Z3+dVol1dt*Z2*Z3*(1-Z1);

eqs=[eqn7,eqn8,eqn9];
vars=[Vol1(t) Vol2(t) Vol3(t)];

[M,F]=massMatrixForm(eqs,vars);

M=odeFunction(M,vars);
F=odeFunction(F,vars,D2,D3,ks2,ks3,V1,V2,V3,Cis,CisH2O,structuralparameter,Z1,Z2,Z3,porosity,S01,r0,x,K1,K2,K3);

T=1073; %CHANGE THIS!

%these are obsolete for the first reaction
ks2=1*(0.5*1.92*10^-5).*exp(-35000./(8.3145.*T)); %assumed, since data is not available
%ks3=1*(0.2*1.92*10^-5).*exp(-35000./(8.3145.*T)); %assumed, since data is not available
ks3=1*(1.1*10^-6).*exp(-64000./(8.3145.*T)); %Liu et al, 2014

D3=1*5*10^(-7)*exp(-90*10^3/(8.3145*T)); %from Liu 2014, page 160 !!check this!!! 
D2=2.56*10^(-6)*exp(-100*10^3/(8.3145*T));
D1=2.56*10^(-6)*exp(-100*10^3/(8.3145*T)); %From New mechanism paper

%Change these for every different temperature!

K1=3.839554194322706e+05; %data from equilibrium plots file, make sure set at the temperature of the experiment. Fe2O3 to Fe3O4, with CO
K2=0.688550105874293; %Fe3O4 to FeO, with CO
K3=1.144396086018002; %FeO to Fe with CO

Ci1(t)=subs(Ci1(t));
Ci2(t)=subs(Ci2(t));
Ci3(t)=subs(Ci3(t));
Ci1H2O(t)=subs(Ci1H2O(t));
Ci2H2O(t)=subs(Ci2H2O(t));
Ci3H2O(t)=subs(Ci3H2O(t));

%These might be slightly different for Tungsten/Iron system. Get Fennell
%paper link for Jose. 

V1=0.159692/5250; %Molar volume. Data from Paul fennell paper. Molar mass/Density ((mol/kg)/(kg/m^3))
V2=0.231533/5100;
V3=0.068885/5600;
V4=0.055485/7870;

S01=1.33*10^7;
S01=8.88*10^6;
S01=7.120109100000000e+06; %specific surface area (m^2/m^3)

x=0.947; %wustite non-stoichiometry

Z1=V2/V1;
Z2=V3/V2;
Z3=V4/V3;


porosity=0.22;
structuralparameter=4.04;
%structuralparameter=4;

r0=sqrt(porosity/(pi()*(structuralparameter*S01^2/(4*pi()*(1-porosity))))); %initial radius of the pores as initial condition

Cis=(0.03*10^5)/(8.3145*T); %gas concentration of fuel gas (i.e. CO)
CisH2O=(0.15*10^5)/(8.3145*T); %%gas concentration of product gas (i.e. CO2)

Cis-CisH2O/K1; %might be obsolete
Cis-CisH2O/K2;
Cis-CisH2O/K3;

Ci1(t)=subs(Ci1(t)); %syntax to get equations into right type of variable for the solver
Ci2(t)=subs(Ci2(t));
Ci3(t)=subs(Ci3(t));

Ci1H2O(t)=subs(Ci1H2O(t));
Ci2H2O(t)=subs(Ci2H2O(t));
Ci3H2O(t)=subs(Ci3H2O(t));

F=@(t,Y) F(t,Y,D2,D3,ks2,ks3,V1,V2,V3,Cis,CisH2O,structuralparameter,Z1,Z2,Z3,porosity,S01,r0,x,K1,K2,K3); %Set of DEs to be solved 

y0=[porosity;porosity;porosity]; %initial condition enclosed volume
opt=odeset('mass',M); %syntax/options for solver
[t Y]=ode45(F,linspace(0,endtime,numdata),y0,opt); %ODE solver

X1=(Y(:,1)-porosity)/(1-porosity);
X2=(Y(:,2)-(1-Z1)*Y(:,1)-Z1*porosity)/(Z1*(1-porosity));
X3=(Y(:,3)-(1-Z2)*Y(:,2)-Z2*(1-Z1)*Y(:,1)-Z1*Z2*porosity)/(Z1*Z2*(1-porosity));

conversion=0.1111*X1+0.1889*X2+0.7*X3;
S=real(100*X1);

end