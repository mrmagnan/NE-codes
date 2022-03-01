% NE 403 data analysis, lab 2, heat conduction
% Code by Declan Miller

close all
clc; clear;
Trad= [ 25.0,23.4,22.1,21.3,20.8,20.3;
        27.2	24.8	22.9	21.7	21.1	20.5;
        32.1	28.4	24.8	22.5	21.5	20.5;
        34.7	30.5	25.9	23.1	21.8	20.6;
       38.1	32.9	27.7	24.0	22.3	20.9];

T0=100; rh=.004; L=.003;
Rrad=[0.01 0.02 0.03 0.04 0.05];
Q= [6.5 9.5 16.5 19.9 25.3];
Trad6=Trad(1,2:end);
Trad9=Trad(2,2:end);
Trad16=Trad(3,2:end);
Trad19=Trad(4,2:end);
Trad25=Trad(5,2:end);

% Apparent results

aRad=[ 25.16 -1.912;
       27.22 -2.677;
       32.81 -4.946;
       35.99 -6.198;
       39.77 -7.583];
for i=1:5
kradcalc(i)=Q(i)/(2*pi*-aRad(i,2)*L);
end
% uncertainties in slope
sigmam=0;

% Converting the confidence intervals into standard errors
CIrad=[-1.993, -1.831;
        -2.855 -2.469;
        -5.551 -4.341;
        -6.972, -5.424;
        -8.432 -6.733];

for i=1:5
    
    for j=1:4
a2pointwise(i,j)=(Trad(i,1)-Trad(i,j+1))/log(Rrad(j)/rh);
kradpointwise(i,j)=Q(i)/(2*pi*a2pointwise(i,j)*L);

 % Error propagation for analytical values
    siga2(i)=0.1/log(Rrad(i)/rh);
    sigQ=0.5;
    sigkradpointwise(i,j)=sqrt((-Q(i)/(2*pi*L*a2pointwise(i,j)^2))^2*siga2(i)+(1/(2*pi*a2pointwise(i,j)*L))^2*0.5^2);
 % Average deviations etc
    sigkradavgptwise(i)=sum(sigkradpointwise(i,:))/4;
    end
    kradptwiseavg(i)=sum(kradpointwise(i,:))/4;
    % Error propagation for calculated values
    SE(i)=(CIrad(i,2)-CIrad(i,1))/3.92;
    sigkradcalc(i)=sqrt((-Q(i)/(2*pi*L*aRad(i,2)^2))^2*SE(i)+(1/(2*pi*aRad(i,2)*L))^2*0.5^2);
  
end

kradcalcavg=sum(kradcalc)/5;
sigkradcalcavg=sum(sigkradcalc)/5;

% Analytical temp distribution;
syms rradi
figure
for i=1:5
    hold on
Tradi(i)=Trad(i,1)-(Q(i)/(2*pi*L*kradcalc(i)))*(log(rradi/rh));
fplot(Tradi(i),[0.004 0.05]);
end
xlabel('Distance from center (mm)');
ylabel('Temperature (C)');
title('Analytical Temperature Distribution with calculated k');
legend('6.5 W','9.5 W','16.5 W','19.9 W','25.3 W');


% K vs average temperature
for i=1:5
    Tavg(i)=sum(Trad(i,:))/5; 
end
figure
plot(Tavg(:),kradcalc(:));
ylabel('Thermal conductivity');
xlabel('Average Temperature');
title('k vs T avg for radial experiment');


syms rn zn;

kradptwiseavgoverall=sum(kradptwiseavg(1,:))/5;
sigptwiseavg=sum(sigkradavgptwise(:))/5;
% Tr=T0-(1/(2*pi*k*L))*log(rn/rh);

Tax= [30.4	29.7	28.8	17.1	16.5	15.7	14.5	13.4	12.7;
42.3	41.2	39.6	20.5	19.2	17.9	16.5	14.8	13.3;
48.4	46.6	43.9	27.7	25.3	22.7	20.1	18.1	16.3;
59.9	57.9	54.7	31.5	28.7	26.1	22.5	19.8	17.4;
67.0	65.0	61.5	33.5	30.4	27.6	23.7	20.7	18.0];

Tax7=Tax(1,2:end);
Tax10=Tax(2,2:end);
Tax15=Tax(3,2:end);
Tax20=Tax(4,2:end);
Tax23=Tax(5,2:end);

Qax=[7.2 9.9 15.4 20.1 22.5];

d=0.025;
A0=pi*(d/2)^2;
zax=0.020:0.010:0.090;
z0=0.01; z1gap=0.035; z2gap=0.065;

aAx=[30.47 81.11 10.79 0.4667;
    42.6 146.7 17.53 -0.06667;
    48.63 225.6 14.38 0.3;
    60.33 268.9 20.81 0.8;
    67.7 296.7 25.33 0.8];

CIax2=[62.86 99.36;
       125.3 168;
       180 271.1;
       240.3 297.5;
       267.4 325.9];
    
CIax3=[10.21 11.37;
       16.86 18.21;
       12.94 15.82;
       19.91 21.72;
       24.41 26.26;
            ];

CIax4=[-.1656 1.099;
       -0.8071 0.6737;
       -1.278 1.878;
       -0.191 1.791;
       -0.2138 1.814;
            ];
        
for i=1:5
SEax2(i)=(CIax2(i,2)-CIax2(i,1))/3.92;    
SEax3(i)=(CIax3(i,2)-CIax3(i,1))/3.92;    
SEax4(i)=(CIax4(i,2)-CIax4(i,1))/3.92;    

% calculated
kaxsigma(i)=sqrt((1/(A0*aAx(i,2)))^2*(0.5)^2+(-Q(i)/(A0*aAx(i,2)^2))^2*SEax2(i)^2);

kaxcalc(i)=Qax(i)/(aAx(i,2)*A0);

% Hg1
Hg1(i)=Qax(i)/(aAx(i,3)*A0);
Hg1sigma(i)=sqrt((1/(A0*aAx(i,3)))^2*(0.5)^2+(-Q(i)/(A0*aAx(i,3)^2))^2*SEax3(i)^2);

% Hg2
Hg2(i)=Qax(i)/(aAx(i,4)*A0);
Hg2sigma(i)=sqrt((1/(A0*aAx(i,4)))^2*(0.5)^2+(-Q(i)/(A0*aAx(i,4)^2))^2*SEax4(i)^2);

end

kaxavg=sum(kaxcalc(:))/5;
kaxsigmavg=sum(kaxsigma(:))/5;
Hg1ax=sum(Hg1(:))/5;
Hg1sigmaax=sum(Hg1sigma(:))/5;
Hg2ax=sum(Hg2(:))/5;
Hg2sigmaax=sum(Hg2sigma(:))/5;
% Pointwise for axial functions

% Analytical temp distribution;
syms zzax
figure
for i=1:5
    hold on
Taxz(i)=Tax(i,1)-(Qax(i)/(kaxcalc(i)*A0)*(zzax-z0)+(Q(i)/(Hg1(i)*A0))*heaviside(zzax-z1gap)+(Q(i)/(Hg2(i)*A0))*heaviside(zzax-z2gap));
fplot(Taxz(i),[0.01 0.1]);
end
xlabel('Position along axis (mm)');
ylabel('Temperature (C)');
title('Analytical Temperature Distribution with calculated k''s, Hg1''s, Hg2''s');
legend('7.2 W','9.9 W','15.4 W','20.1 W','22.5 W');



for i=1:5
a2axpointwise(i)=(Tax(i,1)-Tax(i,2))/(zax(1)-z0);
sigmaa2pt(i)=0.1/(zax(1)-z0);
kpointwiseax(i)=Qax(i)/(A0*a2axpointwise(i));
% pointwise
kaxsigma(i)=sqrt((1/(A0*a2axpointwise(i))^2*(0.5)^2+(-Q(i)/(A0*a2axpointwise(i)^2))^2*sigmaa2pt(i)^2));
end
kptwavgax=sum(kpointwiseax(:))/5;
kptsigavg=sum(kaxsigma(:))/5;
figure;
for i=1:5
    Tavg(i)=sum(Tax(i,:))/5; 
end
plot(Tavg(:),kaxcalc(:));
ylabel('Thermal conductivity');
xlabel('Average Temperature');
title('k vs T avg for axial experiment');
