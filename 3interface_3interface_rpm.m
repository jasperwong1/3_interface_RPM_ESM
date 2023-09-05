clear

Sheetnumber='Sheet1';
data=xlsread('Fe2O3_Fe_15CO_1Cycles_150_300um_QW_900C_29032022.xlsx',Sheetnumber);

T=1173;

Time=data(2:901,2);
Sdata=data(2:901,4);
t=Time;

endtime=Time(end);
numdata=length(Time);

porosity=0.22;
V1=0.159692/5250; %Paul fennell paper
V2=0.231533/5100;
V3=0.068885/5600;
V4=0.055485/7870;

Z1=2.13/2.14; %sourced from paper by fennell hayhurst dennis 2011. Maybe not interpreted correctly here
Z2=1.81/2.13;
Z3=1/1.81;

Z1=V2/V1;
Z2=V3/V2;
Z3=V4/V3;

p0=[3.820607170308743e-09*0.024,1.281721319870089e-11*10.8];
p=[1.48e-9*exp((-25.9*10^3)/(8.3145*T)), 1.281721319870089e-11*10.8]

options = optimoptions('lsqcurvefit','StepTolerance',1e-6);

ub = [1e-5,1e-9];
lb = [1e-11,1e-13];

%ub = [1,1];


%[p,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@(p,t) kineticsfun(p,t)./Sdata,p0,Time,Sdata./Sdata,lb,ub,options);
%ci=nlparci(p,Rsd,'jacobian',Jmat);
results=kineticsfun(p,t);


figure(15)
%plot(Time,Sdata,'o',t,results,'LineWidth',2)
plot(t,Sdata,'o',t,results,'LineWidth',2)
set(gca,'FontSize',16)
l=[legend('TGA Experiment','Random Pore Model') ylabel('Conversion $X$') xlabel('Time /s')]
set(l,'Interpreter','latex')
legend('TGA Experiment','Random Pore Model')
%function S=kineticsfun(t,p)


function S=kineticsfun(p,t)




[t,Vol] = ode45(@odefun,linspace(0.5004,length(t),900),[0.22; 0.22; 0.22]);



function dVoldt=odefun(t,Vol)

r1=[];
r2=[];
r3=[];
r3dash=[];
S1=[];
S2=[];
S3=[];
S3dash=[];
structuralparameter=[]
D1=[];
D2=[];
D3=[];
ks1=[];
ks2=[];
ks3=[];


syms D1 D2 D3 ks1 ks2 ks3 r1 r2 r3 r3dash Ci1 Ci2 Ci3 Ci1H2O Ci2H2O Ci3H2O Cis CisH2O K1 K2 K3 V1 V2 V3 S1 S2 S3 S3dash structuralparameter

%syms D1 D2 D3  Ci1 Ci2 Ci3 Ci1H2O Ci2H2O Ci3H2O Cis CisH2O K1 K2 K3 V1 V2 V3

%{
ks1=p(1);
D1=p(2);
ks2=p(3);
D2=p(4);
ks3=p(5);
D3=p(6);
%}


%p=[1.042429154960800e-08*0.8,4.759700427821800e-11,4.904349059047330e-09*1.1,2.193367297224893e-11*0.5, 2.061772645834602e-09*0.25,3.531545792863986e-12*0.2];

%[4.09583377831375e-09,5.02137397654932e-09;5.09174900898165e-11,5.16263232629492e-11;3.69792628151003e-09,3.94328805910746e-09;6.99285810797068e-12,1.86415682894311e-11]
%p=[4.558603877431536e-09,5.127190667638285e-11,3.820607170308743e-09,1.281721319870089e-11,3.820607170308743e-09*0.02,1.281721319870089e-11*0.8];



T=1173;
Cis=(0.15*10^5)/(8.3145*T);
CisH2O=(0*10^5)/(8.3145*T);

ks1=1.3660e-06*exp(-4.8267e+04./(8.3145*T));
D1=1.7216e-09*exp(-3.1276e+04./(8.3145*T));
ks2=5.751709896875172e-09*0.4;
D2=5.540082985371705e-11*0.4;
%ks2=4.020243100681743e-09;
%D2=2.6950e-11;

K1=7.583567432175545e+04;
K2=1.704654782419324;
K3=0.491055794220954;

K1=7.583567432175545e+04;
K2=1.704654782419324;
K3=0.491055794220954;


x=0.947;

V1=0.159692/5250; %Paul fennell paper
V2=0.231533/5100;
V3=0.068885/5600;
V4=0.055485/7870;

Z1=2.13/2.14; %sourced from paper by fennell hayhurst dennis 2011. Maybe not interpreted correctly here
Z2=1.81/2.13;
Z3=1/1.81;

Z1=V2/V1;
Z2=V3/V2;
Z3=V4/V3;

%{
ks1=p(1);
D1=p(2);
ks2=p(3);
D2=p(4);
ks3=p(5);
D3=p(6);
%}

structuralparameter=4.04;
porosity=0.22;
S01=7.120109100000000e+06;
r0=sqrt(porosity/(pi()*(structuralparameter*S01^2/(4*pi()*(1-porosity)))));

Ci1H2O=Cis+CisH2O-Ci1;
Ci2H2O=Cis+CisH2O-Ci2;
Ci3H2O=Cis+CisH2O-Ci3;

%r1=2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(1))/(1-porosity))))-1)+r0

eqn1=(D3/((abs(r3-r3dash)+r3-r3dash)/2))*(Cis-Ci3)==(D2/((abs(r2-r3)+r2-r3)/2))*(Ci3-Ci2)+ks3*(Ci3-Ci3H2O/K3)/V3;
eqn2=(D2/((abs(r2-r3)+r2-r3)/2))*(Ci3-Ci2)==(D1/((abs(r1-r2)+r1-r2)/2))*(Ci2-Ci1)+((4*x-3)/x)*ks2*(Ci2-Ci2H2O/K2)/V2;
eqn3=(D1/((abs(r1-r2)+r1-r2)/2))*(Ci2-Ci1)==(ks1*(Ci1-Ci1H2O/K1))/(3*V1);



[solCi1,solCi2,solCi3]=solve([eqn1,eqn2,eqn3],[Ci1,Ci2,Ci3]);



Ci1=solCi1;
Ci2=solCi2;
Ci3=solCi3;

Ci1H2O=Cis+CisH2O-Ci1;
Ci2H2O=Cis+CisH2O-Ci2;
Ci3H2O=Cis+CisH2O-Ci3;


dVoldt=zeros(3,1);
%Vol=zeros(3,1)



Vol1dash=Vol(1)-Z1*(Vol(1)-porosity);
Vol2dash=Vol(2)-Z2*(Vol(2)-Vol1dash);
Vol3dash=Vol(3)-Z3*(Vol(3)-Vol2dash);

%{
rvector=[2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(1))/(1-porosity))))-1)+r0; 2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(2))/(1-porosity))))-1)+r0; 2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(3))/(1-porosity))))-1)+r0; 2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity))))-1)+r0]

Ci1=subs(Ci1,[r1;r2;r3;r3dash],rvector)
Ci2=subs(Ci2,[r1;r2;r3;r3dash],rvector);
Ci3=subs(Ci3,[r1;r2;r3;r3dash],rvector);
%}

%{
S1=S01*(1-Vol(1))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol(1))/(1-porosity)));
S2=S01*(1-Vol(2))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol(2))/(1-porosity)));
S3=S01*(1-Vol(3))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol(3))/(1-porosity)));
S3dash=S01*(1-Vol3dash)/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity)));
%}
rvector=[2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(1))/(1-porosity))))-1)+r0;
    2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(2))/(1-porosity))))-1)+r0; 
    2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol(3))/(1-porosity))))-1)+r0; 
    2/(S01*structuralparameter/(1-porosity))*((sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity))))-1)+r0;
    S01*(1-Vol(1))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol(1))/(1-porosity)));
    S01*(1-Vol(2))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol(2))/(1-porosity)));
    S01*(1-Vol(3))/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol(3))/(1-porosity)));
    S01*(1-Vol3dash)/(1-porosity)*sqrt(1-structuralparameter*log((1-Vol3dash)/(1-porosity)));
    p(1);
    p(2)];

%{
dVoldt(1)=(ks1*(Ci1-Ci1H2O/K1))*S1
dVoldt(2)=(1-Z1)*dVoldt(1)+ks2*(Ci2-Ci2H2O/K2)*S2
dVoldt(3)=(1-Z2)*dVoldt(2)+Z2*(1-Z1)*dVoldt(1)+ks3*(Ci3-Ci3H2O/K3)*S3
%}


    
dVoldt(1)=subs((ks1*(Ci1-Ci1H2O/K1))*S1,[r1;r2;r3;r3dash;S1;S2;S3;S3dash;ks3;D3],rvector);
dVoldt(2)=subs((1-Z1)*dVoldt(1)+ks2*(Ci2-Ci2H2O/K2)*S2,[r1;r2;r3;r3dash;S1;S2;S3;S3dash;ks3;D3],rvector);
dVoldt(3)=subs((1-Z2)*dVoldt(2)+Z2*(1-Z1)*dVoldt(1)+ks3*(Ci3-Ci3H2O/K3)*S3,[r1;r2;r3;r3dash;S1;S2;S3;S3dash;ks3;D3],rvector);

    %}
    
%{
dVoldt(2)=(1-Z1)*dVoldt(1)+ks2*(Ci2-Ci2H2O/K2)*S2
dVoldt(3)=(1-Z2)*dVoldt(2)+Z2*(1-Z1)*dVoldt(1)+ks3*(Ci3-Ci3H2O/K3)*S3

toc

%}
end

porosity=0.22;

V1=0.159692/5250; %Paul fennell paper
V2=0.231533/5100;
V3=0.068885/5600;
V4=0.055485/7870;

Z1=2.13/2.14; %sourced from paper by fennell hayhurst dennis 2011. Maybe not interpreted correctly here
Z2=1.81/2.13;
Z3=1/1.81;

Z1=V2/V1;
Z2=V3/V2;
Z3=V4/V3;

X1=(Vol(:,1)-porosity)/(1-porosity);
X2=(Vol(:,2)-(1-Z1)*Vol(:,1)-Z1*porosity)/(Z1*(1-porosity));
X3=(Vol(:,3)-(1-Z2)*Vol(:,2)-Z2*(1-Z1)*Vol(:,1)-Z1*Z2*porosity)/(Z1*Z2*(1-porosity));

figure(1)
plot(t,[X1,X2,X3],'LineWidth',2)
legend('X_1','X_2','X_3')
xlabel('Time /s','FontSize',16)
ylabel('Conversion','FontSize',16)
xlim([0 40])
set(gca,'FontSize',16)
xlabel('Time /s')
ylabel('Conversion X_j')

%conversion=0.1111*X1+0.1889*X2+0.7*X3;
S=real(0.1111*X1+0.1889*X2+0.7*X3)*100;


Vol1dash=Vol(:,1)-Z1*(Vol(:,1)-porosity);
Vol2dash=Vol(:,2)-Z2*(Vol(:,2)-Vol1dash);
Vol3dash=Vol(:,3)-Z3*(Vol(:,3)-Vol2dash);

structuralparameter=4.04;
S01=7.120109100000000e+06;

S1=S01.*(1-Vol(:,1))./(1-porosity).*sqrt(1-structuralparameter.*log((1-Vol(:,1))./(1-porosity)));
S2=S01.*(1-Vol(:,2))./(1-porosity).*sqrt(1-structuralparameter.*log((1-Vol(:,2))./(1-porosity)));
S3=S01.*(1-Vol(:,3))./(1-porosity).*sqrt(1-structuralparameter.*log((1-Vol(:,3))./(1-porosity)));
S3dash=S01.*(1-Vol3dash)./(1-porosity).*sqrt(1-structuralparameter.*log((1-Vol3dash)./(1-porosity)));


figure(8)
yyaxis left
plot(t,S1,t,S2,t,S3,'LineWidth',2)
l=[ylabel('Surface area $S_{i}$ /m$^{-1}$') xlabel('Time /s')]
ylim([0 8*10^6])
xlim([0 150])
set(l,'Interpreter','latex')
yyaxis right
plot(t,Vol(:,1),t,Vol(:,2),t,Vol(:,3),'LineWidth',2)
ylim([0 1])
xlim([0 150])
l=[legend('$S_{1}$','$S_{2}$','$S_{3}$','$V_{1}$','$V_{2}$','$V_{3}$'),ylabel('Enclosed volume $V_{i}$ (-)') xlabel('Time /s')]
set(gca,'FontSize',16)
set(l,'Interpreter','latex')


end
