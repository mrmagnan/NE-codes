% NE 400, Lab 3, Heat Convection
% March 2022
% Code by Declan Miller

% NOTE: Requires data sheet (excel) 'Lab3TableOneValues.xlsx' in order to run.
% This data sheet simply corresponds to the values given for experiment D in the given data sheet.
% However, the sheet has been modified, in order to make it more easily readable in matlab.

clc; clear;
Tab1=table2array(readtable('Lab3TableOneValues.xlsx'));

% Experiment D data
T1=70;
T2=Tab1(1,:);
T3=Tab1(2,:);
T4=Tab1(3,:);
T5=Tab1(4,:);
T6=Tab1(5,:);
T7=Tab1(6,:);
T8=Tab1(7,:);
T9=Tab1(8,:);
T10=Tab1(9,:);
Mflowhot=Tab1(10,:)/60; % L/min = kg/min --> kg/s, hot
Mflowcool=30/1000; % g/s --> kg/s, coolant
Hotout=T2;
Hotin=T3;

% Experiment C data
HotoutC=[33.9	34.3	33.9	34.9	33.5	35.8	35.6	37.1	35.2	37.4];
HotinC=[36.4	36.1	36.5	37.5	38.1	39.7	40.9	42.4	42.7	44.2];
T7C=[29.9	18.6	30.3	20.5	30.5	17	31.9	17.2	31.5	18.5];
T10C=[19.4	30.1	18.8	30.2	16.4	30.5	17.9	31.2	17.2	31.3];
HotflowC=[8	8	6	6	4	4	3	3	2	2]./60;
ColdflowC=[0.02 0.02 0.02 0.02 0.016 0.016 0.016 0.016 0.016 0.016];

Cp=4190; % specific heat of water [J/(kg*C)]

k=0.651; % thermal conductivity of copper in W/m

% Getting delta T's
deltaT=[];
Qi=[];
Qout=[];
G=[];
mu= 4.67*10^-4; % Dynamic viscosity of water in Pa*s

di=0.0079;
Ainner=di*3.1415*.870; % inner area
AXS=3.1415*(di/2)^2; % Cross Sectional Area

Pr=Cp*mu/k; % Prandl Number
j=0; % an index
z=0; % an index
for i=1:12
    % Experiment D
    deltaT(i)=((T3(i)-T7(i))-(T2(i)-T10(i)))/log((T3(i)-T7(i))/(T2(i)-T10(i)));
    Qi(i)=Mflowhot(i)*Cp*(T3(i)-T2(i));
    Qout(i)=Mflowcool*Cp*abs(T7(i)-T10(i));
    % heat transfer
    h(i)=Qi(i)/(Ainner*deltaT(i));
    G(i)=Mflowhot(i)/AXS;
    Nu(i)=h(i)*di/k;
    Re(i)=G(i)*di/mu;
    
    if rem(i,2) ==0
    hparD(j)=h(i);
    j=j+1;
    Nupar(j)=Nu(i);
    Repar(j)=Re(i);
    LMTDpar(j)=deltaT(i);
    Qipar(j)=Qi(i);
    Qoutpar(j)=Qout(i);
    else        
    z=z+1;
    hcountD(z)=h(i);
    Nu_countflow_D(z)=Nu(i);
    Re_countflow_D(z)=Re(i);
    LMTDcount(z)=deltaT(i);
    Qicount(z)=Qi(z);
    Qoutcount(z)=Qout(z);
    end
    
    % Errors:
    
    % Error h
    dhdqi(i)=1/(Ainner*deltaT(i));
    SigQi=0.25/60;
    dhdLMTD(i)=-Qi(i)/(Ainner*deltaT(i)^2);
    SigLMTD=0.1;
    error_h(i)=sqrt(dhdqi(i)^2*SigQi^2+dhdLMTD(i)^2*SigLMTD^2);
    errorNu(i)=error_h(i)*di/k;
    errorQi(i)=SigQi*Cp*(T3(i)-T2(i));
    
    SigG(i)=SigQi/AXS;
    errorRe(i)=SigG(i)*di/mu;
end

% Averaging some values yo
havgDcount=sum(hcountD)/length(hcountD);
havgDcounterrorD=sum(error_h)/length(error_h);
havgDpar=sum(h)/length(h);
havgDpar=sum(error_h)/length(error_h);


j=0;
z=0;
for i=1:10
    % Experiment C
    deltaTC(i)=((HotinC(i)-T7C(i))-(HotoutC(i)-T10C(i)))/log((HotinC(i)-T7C(i))/(HotoutC(i)-T10C(i)));
    QiC(i)=HotflowC(i)*Cp*(HotinC(i)-HotoutC(i));
    QoutC(i)=ColdflowC(i)*Cp*abs(T7C(i)-T10C(i));
    % heat transfer
    hC(i)=QiC(i)/(Ainner*deltaTC(i));
    GC(i)=HotflowC(i)/AXS;
    NuC(i)=hC(i)*di/k;
    ReC(i)=GC(i)*di/mu;
    
    if rem(i,2) ==0
    j=j+1;
    NuparC(j)=NuC(i);
    ReparC(j)=ReC(i);
    LMTDparC(j)=deltaTC(i);
    QiparC(j)=QiC(i);
    QoutparC(j)=QoutC(i);
    else        
    z=z+1;
    Nu_countflow_C(z)=NuC(i);
    Re_counterflow_C(z)=ReC(i);
    LMTDcountC(z)=deltaTC(i);
    QicountC(z)=QiC(z);
    QoutcountC(z)=QoutC(z);
    end
    
    % Errors:
    
    % Error h
    dhdqiC(i)=1/(Ainner*deltaTC(i));
    SigQiC=0.25/60;
    dhdLMTDC(i)=-QiC(i)/(Ainner*deltaTC(i)^2);
    SigLMTDC=0.1;
    error_hC(i)=sqrt(dhdqiC(i)^2*SigQiC^2+dhdLMTDC(i)^2*SigLMTDC^2);
    errorNuC(i)=error_hC(i)*di/k;
    errorQiC(i)=SigQiC*Cp*(HotinC(i)-HotoutC(i));
    
    SigGC(i)=SigQi/AXS;
    errorReC(i)=SigG(i)*di/mu;
end

%% Plotting

logReC=log(Re_counterflow_C);
logNuC=log(Nu_countflow_C);

plot(logReC,logNuC,'o-');
xlabel('ln[Re]');
ylabel('ln[Nu]');
title('ln[Nu] vs ln[Re], experiment C');

figure
logReD=log(Re_countflow_D);
logNuD=log(Nu_countflow_D);
plot(logReD,logNuD,'o-');
xlabel('ln[Re]');
ylabel('ln[Nu]');
title('ln[Nu] vs ln[Re], experiment D');


