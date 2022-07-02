function dY = modelA2(t, Y)

% function for model 3.5 with A2 (unphosphorylated dimer) promoting
% expression of AB (autoreg) & sigma E

global kpdeg kmdeg ktlnA ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1 f2 f3 K1 K2 K3 kb kd krdeg
global ktlnB k1 km1 k2 km2 k3 km3 k4 k5 km5 k6
global f5b K5b

A=Y(1);
Ap=Y(2);
E=Y(3);
R=Y(4);
ER=Y(5);
mAB=Y(6);
mE=Y(7);
Bs = Y(8);
B = Y(9);
Bp = Y(10);
kpa = Y(11);
pap = Y(12);
sigB = Y(13);
% K3 = K1;
% A Ap E R ER mAB mE Bs B Bp kpa pap
%---------------

dY(1) = ktlnA*mAB-k3*Bp*A + km3*kpa +k6*pap - kpdeg*A; % A
dY(2) = k4*kpa - k5*B*Ap + km5*pap-kpdeg*Ap; % Ap
dY(3)=  ktlnE*mE-kb*E*R+kd*ER+krdeg*ER-kpdeg*E;%E
dY(4)=  ktlnR-kb*E*R+kd*ER-(kpdeg+krdeg)*R;%R
dY(5)=  kb*E*R-kd*ER-(kpdeg+krdeg)*ER;%ER
dY(8) = k1*B - km1*Bs - k2*Bs + km2*Bp + k4*kpa - kpdeg*Bs; % Bs
dY(9) = ktlnB*mAB + km1*Bs - k1*B - k5*B*Ap + km5*pap + k6*pap - kpdeg*B; % B
dY(10) = k2*Bs - km2*Bp + km3*kpa - k3*Bp*A-kpdeg*Bp; % Bp
dY(11) = k3*Bp*A - km3*kpa - k4*kpa-kpdeg*kpa; % KpA
dY(12) = k5*B*Ap - km5*pap - k6*pap-kpdeg*pap; % PAp

dY(6)= kbtpn1*(1+ f1*Ap^2/K1)/(1+Ap^2/K1)+kbtpn2*(f2*E/K2)/(1+E/K2)-1*kmdeg*mAB;%mAB
dY(7)= kbtpn3*(1+f3*Ap^2/K1+ f5b*A^2/K5b)/(1+Ap^2/K1 + A^2/K5b)-kmdeg*mE;%mE
% dY(13) = kb4*(1+f6*E*Ap^2/(KB) + 5*E/(K2))/(1+E*Ap^2/(KB)+E/(K2))-kmdeg*sigB;
dY(13) = kbtpn1*(E/(0.00001*K2))/(1+E/(0.00001*K2)) + kbtpn1*25*(E/(K2))/(1+E/(K2))-4*kmdeg*sigB;
dY=dY';
end