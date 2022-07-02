function dY = modelA2_DnaK(t,Y)

% MprAB-SigE model with DnaK. Different activation mechanism compared 
% to modelA2. MprB:DnaK complex is phosphatase active, free MprB is
% kinase active. S represents unfolded protein load. ST is a parameter
% modeled explicitly to bind DnaK.

global kpdeg kmdeg ktlnA ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1 f2 f3 K1 K2 K3 kb kd krdeg
global ktlnB k2 km2 k3 km3 k4 k5 km5 k6
global f5b K5b f1b f3b
global ktlnD kbdf kbdb kdsf kdsb dnakoverex

A=Y(1);
Ap=Y(2);
E=Y(3);
R=Y(4);
ER=Y(5);
mAB=Y(6);
mE=Y(7);
B = Y(8);
D = Y(9);
Bp = Y(10);
kpa = Y(11);
pap = Y(12);
sigB = Y(13);
BD = Y(14);
DS = Y(15);
S = Y(16);

    % MprAB-DnaK
dY(1) = ktlnA*mAB-k3*Bp*A + km3*kpa +k6*pap - kpdeg*A; % A
dY(2) = k4*kpa - k5*BD*Ap + km5*pap-kpdeg*Ap; % Ap
dY(8) = ktlnB*mAB - k2*B + km2*Bp + k4*kpa -kbdf*B*D + kbdb*BD - kpdeg*B; % B
dY(9) = + ktlnD + dnakoverex*0.03*mAB -kbdf*B*D + kbdb*BD - kdsf*D*S + kdsb*DS - kpdeg*D; %DnaK +0.03*mAB for DnaK overexpression model
dY(14) = kbdf*B*D - kbdb*BD - k5*BD*Ap + km5*pap + k6*pap - kpdeg*BD; %MprB.DnaK
dY(15) = kdsf*D*S - kdsb*DS - kpdeg*DS; % DnaK.SDS %
dY(16) = +kpdeg*DS -kdsf*D*S + kdsb*DS; % SDS %
dY(10) = k2*B - km2*Bp + km3*kpa - k3*Bp*A-kpdeg*Bp; % Bp
dY(11) = k3*Bp*A - km3*kpa - k4*kpa-kpdeg*kpa; % KpA
dY(12) = k5*BD*Ap - km5*pap - k6*pap-kpdeg*pap; % MprB.DnaK.MprAp
    % RseA-SigE
dY(3)=  ktlnE*mE-kb*E*R+kd*ER+krdeg*ER-kpdeg*E;%E
dY(4)=  ktlnR-kb*E*R+kd*ER-(kpdeg+krdeg)*R;%R
dY(5)=  kb*E*R-kd*ER-(kpdeg+krdeg)*ER;%ER
    % transcription
dY(6)= kbtpn1*(1 + f1*Ap^2/K1 + f1b*A^2/K5b)/(1 + Ap^2/K1 + A^2/K5b)+kbtpn2*(f2*E/K2)/(1+E/K2)-kmdeg*mAB;%mAB
dY(7)= kbtpn3*(1 + f3*Ap^2/K3 + f3b*A^2/K5b)/(1 + Ap^2/K3 + A^2/K5b)-kmdeg*mE;%mE
dY(13) = 0.05*kbtpn1 + kbtpn1*(E/(K2*0.2))/(1+E/(K2*0.2))- kmdeg*sigB;
dY=dY';
end