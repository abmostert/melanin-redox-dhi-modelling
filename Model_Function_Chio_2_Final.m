% Copyright (c) 2025 A. Bernardus Mostert
% Licensed under the MIT License. See the LICENSE file in the repository root.


function F = Model_Function_Chio_2_Final(x)

%Variables are: 
% V - volume of solution in ml
% Vmin - the initial volume added in ml
% Vstep - the size of a the volume step in ml
% Vmax - the final volume added to the solution in ml
% V0 - the initial volume in ml
% C - the concentration of the acid/base being added
% CH2Q0 - The idealised initial concentration of the hydroquinone species before dissociation in Molar
% CQ0 - The idealised initial concentration of the quinone species before dissociation in Molar
% H - concentration of hydronium in Molar
% OH - concentration of hydroxyl ion in Molar
% Q - concentration of quinone molecule in Molar
% Q2 - concentration of doubly deprotonated hydroquinone ion in Molar
% HQ - concentration of deprotonated hydroquinone ion in Molar
% H2Q - concentration of hydroquinone molecule in Molar
% SQ - concentration of semiquinone ion in Molar
% HQI - concentration of protonated quinone methide
% QI - concentration of deprotonated quinone methide
% HSQ - concentration of the protonated semiquinone molecule in Molar
% EQ2SQ - reduction potential for semiquinone to hydroquinone
% ESQQ - reduction potential for quinone to semiquinone
% pQ - -log10Q
% logOH - log10OH
% logQ - log10Q
% logQ2 - log10Q2
% logHQ - log10HQ
% logH2Q - log10H2Q
% logSQ- log10SQ
% logHSQ - log10HSQ
% logHQI - log10HQI
% logQI - log10QI
% pH - pH of the system
% E - potential of the system in V
% Kw - water dissociation constant
% pKw - -log10 Kw
% Ka1 - 1st dissociation constant for hydroquinone
% Ka2 - 2nd dissociation constant for hydroquinone
% Kar - dissociation constant for semiquinone
% Kt - equilibrium constant for tautomer formation
% Kaqi - dissociation contant for quinone methide


%creating a global set of parameters that will be defined in this m file as
%well as the function file.
global V Vmin Vstep Vmax V0 C CH2Q0 CQ0 H OH pH E Kw pKw pQ
global Q Q2 HQ H2Q SQ HSQ EQ2SQ ESQQ Ka1 Ka2 Kar HQI QI Kt Kaqi
global logOH logQ logQ2 logHQ logH2Q logSQ logHSQ logHQI logQI

%defining my 3 unknown parameters
E = x(1) ;
pH = x(2) ;
pQ = x(3) ;

%Determining the concentration of H+ and OH-
H = 10.^-pH ;
pKw = 14 ;
Kw = 10.^-pKw ;
OH = Kw./H ;
Ka1 = 10^(-9.56) ;
Ka2 = 10^(-13.13) ;
Kar = 10^(-7) ;
Kt = 3.94 ;
Kaqi = 10.^(-6) ;
EQ2SQ = -0.0981 ; %-0.2169 ;
ESQQ = 0.165; %0.147 ;

%setting the constraints on parameters
Q = 10.^-pQ ;
HQI = Q./Kt ;
QI = (HQI.*Kaqi)./H ;
SQ = Q.*10.^((ESQQ-E)/0.059) ;
HSQ = (SQ.*H)./Kar ;
Q2 = SQ.*10.^((EQ2SQ-E)/0.059) ;
HQ = (H.*Q2)./Ka2 ;
H2Q = HQ.*H./Ka1 ;

%Writing out the 3 balance equations to solve for.
%charge balance, note, to do an acidic calcaultion, make the last term -ve.
%For base, make it +ve10
F = [(H - OH - SQ - HQ - 2.*Q2 - QI - (C*V)./(V0+V));
    
    %Concentration balance for Q like moeities
    (H2Q + HQ + Q2 + HSQ + SQ + Q + HQI + QI - (V0.*(CH2Q0+CQ0))./(V0+V)  );
    
    %electron balance
    (H2Q + HQ + Q2 - Q - HQI - QI - (V0.*(CH2Q0-CQ0))./(V0+V) )  ] ;

%Converting results to log
logOH = -log10(OH) ;
logQ = -log10(Q) ;
logQ2 = -log10(Q2) ;
logHQ = -log10(HQ) ;
logSQ = -log10(SQ) ;
logH2Q = -log10(H2Q) ;
logHSQ = -log10(HSQ) ;
logHQI = -log10(HQI) ;
logQI = -log10(QI) ;