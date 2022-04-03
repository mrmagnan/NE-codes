%% LAB 4
% NE 400
% Condensation
% Declan Miller
clc;clear;close all;

%% Data

% conversions
cm2m=.01;
mm2m=.001;
g2kg=0.001;

% Dropwise Data
ChmbPd=[0	0	0	0	0]; % psat/kN/m^2
SatTd=[99.97	99.97	99.97	99.97	99.97]; % T5; [C]
SurfTd=[79.35	86.12	92.67	97.94	103.9]; % T6; [C]
H20outTd=[36.6	38.6	38.8	40.1	42.5]; % T7; [C]
H20inTd=[30.1	30.5	29.6	29.3	28.8]; % T8; [C]
DripRtd=[330	340	360	340	320]/60; % Drops/min-->Drops/sec
CoolFloRtd=[35	30	25	20	15]*g2kg; % g/s-->kg/s

DropMat=[ChmbPd; SatTd; SurfTd; H20outTd; H20inTd; DripRtd; CoolFloRtd];

% Filmwise Data
ChmbPf=[0	0	0	0	0];
SatTf=[99.97	99.97	99.97	99.97	99.97];
SurfTf=[60.33	62.75	72.26	76.28	84.41];
H20outTf=[28.7	28.9	29.4	31.6	34.5];
H20inTf=[26.1	26.5	26.7	27.3	27.9];
DripRtf=[130	120	150	140	190]/60;
CoolFloRtf=[40	35	30	25	19]*g2kg;

FilmMat=[ChmbPf; SatTf; SurfTf; H20outTf; H20inTf; DripRtf; CoolFloRtf];

% Certain other values

DropDm=2*mm2m; %mm-->m
DropR=DropDm/2;
hf=419*1000; % kJ/kg-->J/kg
hg=2676*1000; %kJ/kg-->J/kg
hfg=hg-hf;

%Errors

ErrDripRt=10/60;
ErrT=0.1;
ErrCoolFlo=1*g2kg;
ErrChamP=0;
ErrDropDm=0.1*mm2m;
Errhf=0;
Errhg=0;

% Experimental constants
CondL=90*mm2m; % mm-->m
CondThk=0.71*mm2m; % mm-->m
CondD=12.7*mm2m; % mm-->m
CondSA=37*(cm2m)^2; % cm^2-->m^2
intVsteam=1840*(cm2m)^3; % cm^3-->m^3
intDchmb=76*mm2m; % mm-->m
kCond=360; % W/(m*K);
H20capNorm=500*(cm2m)^2; % cm^3-->m^3
SAheatelement=144*(cm2m)^2; % cm^2-->m^2
HeatLssChmb=2.5; % W/K
Patm=101.3; % kN/m^2
cp=4180; % J/(kg*K)
rho=1000; % kg/m^3

%% Dropwise Calculations

% Heat transfer rate in condensor
for i=1:5;
Qd(i)=CoolFloRtd(i)*cp*(H20outTd(i)-H20inTd(i));
% Heat flux
HtFlxd(i)=Qd(i)/CondSA;
% T drop across shell
dTshelld(i)=HtFlxd(i)*CondThk/kCond;
% Corrected steam to surface T diff
dTcrctd(i)=SatTd(i)-SurfTd(i);
% surface heat transfer coeff
hsurfd(i)=HtFlxd(i)/dTcrctd(i);
% Drop heat transfer rate
Qdripd(i)=DripRtd(i)*rho*(4/3*pi*DropR^3)*hfg;
end


% Errors
for i=1:5;
ErrQd(i)=sqrt(ErrCoolFlo^2*(cp*(H20outTd(i)-H20inTd(i)))^2+2*ErrT^2*(CoolFloRtd(i)*cp)^2);
ErrHtFlxd(i)=ErrQd(i)/CondSA;
ErrdTshelld(i)=ErrHtFlxd(i)*CondThk/kCond;
ErrdTcrctd(i)=sqrt(2*ErrT^2);
Errhsurfd(i)=sqrt(ErrHtFlxd(i)^2*(1/dTcrctd(i))^2+ErrdTcrctd(i)^2*(HtFlxd(i)/dTcrctd(i)^2)^2);
ErrQdripd(i)=sqrt(ErrDripRt^2*(rho*(4/3*pi*DropR^3)*hfg)^2+(ErrDropDm/2)^2*(3*DripRtd(i)*rho*(4/3*pi*DropR^2)*hfg)^2);
end

% Table
Runs=[1; 2; 3; 4; 5];
DropResults=table(Runs, Qd(:),ErrQd(:), HtFlxd(:),ErrHtFlxd(:), dTshelld(:),ErrdTshelld(:), dTcrctd(:),ErrdTcrctd(:),Runs, hsurfd(:),Errhsurfd(:), Qdripd(:),ErrQdripd(:));
DropResults.Properties.VariableNames={'Test No.' 'Q-dot [J/s]' 'Error Q-dot' 'q" [J/(s*m^2)]' 'Error q"' 'dT shell [C]' 'Error dT shell' 'dT corrected [C]' 'Error dT corrected' 'Test No. ' 'h_surface [W/(m^2*C]' 'Error h_surface' 'Q-dot, drip rate [J/s]' 'Error Q-dot, drip rate'}
%% Filmwise calculations

for i=1:5;
Qf(i)=CoolFloRtf(i)*cp*(H20outTf(i)-H20inTf(i));
% Heat flux
HtFlxf(i)=Qf(i)/CondSA;
% T drop across shell
dTshellf(i)=HtFlxf(i)*CondThk/kCond;
% Corrected steam to surface T diff
dTcrctf(i)=SatTf(i)-SurfTf(i);
% surface heat transfer coeff
hsurff(i)=HtFlxf(i)/dTcrctf(i);
% surface heat transfer coeff Nusselt calc
mu=0.0002822; % N*s/m^2, dynamic viscosity at T_sat
kSatH20=0.6634; % Saturated H20 Thermal Conductivity ((W/(m*K))
hsurfNu(i)=(2/3)*sqrt(2)*(kSatH20^3*rho^2*hfg*9.81/(CondL*mu*(SatTf(i)-SurfTf(i))))^(1/4);
% Drop heat transfer rate
Qdripf(i)=DripRtf(i)*rho*(4/3*pi*DropR^3)*hfg;
end

% Errors

for i=1:5;
ErrQf(i)=sqrt(ErrCoolFlo^2*(cp*(H20outTf(i)-H20inTf(i)))^2+2*ErrT^2*(CoolFloRtf(i)*cp)^2);
ErrHtFlxf(i)=ErrQf(i)/CondSA;
ErrdTshellf(i)=ErrHtFlxf(i)*CondThk/kCond;
ErrdTcrctf(i)=sqrt(2*ErrT^2);
Errhsurff(i)=sqrt(ErrHtFlxf(i)^2*(1/dTcrctf(i))^2+ErrdTcrctf(i)^2*(HtFlxf(i)/dTcrctf(i)^2)^2);
ErrQdripf(i)=sqrt(ErrDripRt^2*(rho*(4/3*pi*DropR^3)*hfg)^2+(ErrDropDm/2)^2*(3*DripRtf(i)*rho*(4/3*pi*DropR^2)*hfg)^2);
end

% Table

FilmResults=table(Runs, Qf(:),ErrQf(:), HtFlxf(:),ErrHtFlxf(:), dTshellf(:),ErrdTshellf(:), dTcrctf(:),ErrdTcrctf(:),Runs, hsurfNu(:), hsurff(:),Errhsurff(:), Qdripf(:),ErrQdripf(:));
FilmResults.Properties.VariableNames={'Test No.' 'Q-dot [J/s]' 'Error Q-dot' 'q" [J/(s*m^2)]' 'Error q"' 'dT shell [C]' 'Error dT shell' 'dT corrected [C]' 'Error dT corrected' 'Test No. ' 'h_surface, from Nusselt eqn [W/(m^2*C]' 'h_surface [W/(m^2*C]' 'Error h_surface'  'Q-dot, drip rate [J/s]' 'Error Q-dot, drip rate'}

%% Plotting/Data Presentation

% Heat flux versus Delta T (corrected)

plot(dTcrctd(:),HtFlxd(:),'-o');
hold on
plot(dTcrctf(:),HtFlxf(:),'-o');
xlabel('Delta T, corrected [C]');
ylabel('Heat Flux [J/(m^2*s)]');
title('Heat flux versus Delta T (corrected)');
legend('Dropwise','Filmwise');
axis([0 max(dTcrctf(:)) 0 max(HtFlxd(:))]);

% Heat Transfer coefficient versus Delta T (corrected)

figure

plot(dTcrctd(:),hsurfd(:),'-o');
% Leftmost point excluded due to axis limits, note its value in body of text.
hold on
plot(dTcrctf(:),hsurff(:),'-o');
xlabel('Delta T, corrected [C]');
ylabel('Heat transfer coefficient, surface [W/(m^2*K)]');
title('Heat Transfer coefficient versus Delta T (corrected)');
legend('Dropwise','Filmwise');
axis([0 max(dTcrctf(:)) 0 max(hsurfd(:))]);

% Filmwise, comparing h_surf from cond to that from Nusselt
% Heat Transfer coefficient versus Delta T (corrected)

figure

plot(dTcrctf(:),hsurff(:),'-o');
hold on
plot(dTcrctf(:),hsurfNu(:),'-o');
xlabel('Delta T, corrected [C]');
ylabel('Heat transfer coefficient, surface [W/(m^2*K)]');
title('Heat Transfer coefficient versus Delta T (corrected); Filmwise');
legend('From Condensation Experiment','From Nusselt Equation');

axis([min(dTcrctf(:)) max(dTcrctf(:)) 0 max(hsurfNu(:))]);
