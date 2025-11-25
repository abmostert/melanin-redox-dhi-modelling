% Copyright (c) 2025 A. Bernardus Mostert
% Licensed under the MIT License. See the LICENSE file in the repository root.


clear all
close all

%estimating the redox standard potentials in DHI melanin.

%First we determing Kt, this is done based upon the numbers given by the
%Farmer 2002 paper at pH = 6.3. Kt depends on the amount of H2Q monomer,
%so, I gave a range of values commensurate with the litrature (Gillette & Serpintini)
H = 10^(-6.3) ;
Ka2 = 10^(-9.54) ;
Ka1 = 10.^(-13.09) ;
Kr1 = 10.^(-6.8) ;


%here I set the fraction of H2Q moeities in the system
H2Q_C = 0.95;%linspace(0.01,0.99,1001) ;




%fully optimised data set requires pH = 6.6, pKa2 = 9.61 (weighted average
%of DHI and DHICA pKa2), pKa1 = 13.12 (weighted average), pKr1 = 7.7, pKaqi
%= 6.6, H2Q_C = 0.01, and 5,6, Indolequinone reduction potential at 0.7 V (assuming full protonation).

%We now estimate Kcompprime, the euqilibiurm constant for the
%comproportionation reaction. This is done at pH=6.3 using the numbers
%given by Farmer paper.

Kcompprime =  (1./H2Q_C.^2) ./ (    (((1./H2Q_C).*(1-1/1500))-(1+Ka2/H)).*(1+Ka2/H).*(1500.^2) )      ;

figure(1)
plot(H2Q_C,-log10(Kcompprime),'b')

%Now I need to determine the true Kcomp, which assumes dissociation of H2Q,
%HSQ HQM and tautormerisation. Doing this assuming the pH = 6.3



Kcomp = Kcompprime  .*  ((Kr1^2)/(Ka2*Ka1)).*( (Ka2*Ka1+Ka2*H+H^2)  ./ ((Kr1+H)^2) ) ;


figure(2)
plot(H2Q_C,-log10(Kcomp),'b')


%now, the standard redox potential of 5,6, Indolequinone, 2
%electron rduction is about E3H = 0.7V (Han2018), accounting for the protonation. We also know that
%E3H = (EHSQQ + EQH2SQ)/2. I can determine at least
%ESQQ + EQ2SQ = ESQQ_EQ2SQ, the non pH dependent deprotonated standard
%redox potentials

E3H = 0.7 ;

ESQQ_EQ2SQ = 2.*E3H + 0.059.*log10(Kr1) + 0.059.*log10(Ka2.*Ka1./Kr1);


figure(3)
plot(H2Q_C,ESQQ_EQ2SQ,'b')

%We also know that ESQQi - EQ2SQi = 0.059log10(Kcomprime) = ESQQi_EQ2SQi_minus at a paritcular
%pH (in this case pH = 6.3), I can then estimate the ESQQ - EQ2SQ =
%ESQQ_EQ2SQ_minus, the non pH dependent deprotonated standard
%redox potentials
H = 10.^(-6.3) ;

ESQQi_EQ2SQi_minus = 0.059.*log10(Kcompprime) ;

ESQQ_EQ2SQ_minus = ESQQi_EQ2SQi_minus + 0.059.*log10(Kr1) - 0.059.*log10((Kr1 + H)) - 0.059.*log10(Ka2.*Ka1./Kr1) + 0.059.*log10( (Ka2.*Ka1 + Ka2.*H + H.^2)./(Kr1 + H) ) ;


figure(4)
plot(H2Q_C,ESQQ_EQ2SQ_minus,'b' )

%I can now determine the pH independent ESQQ and EQ2Q standard redution
%potentials.

ESQQ = (ESQQ_EQ2SQ_minus + ESQQ_EQ2SQ)./2 ;
EQ2SQ = (ESQQ_EQ2SQ - ESQQ_EQ2SQ_minus)./2 ;


figure(5)
hold on
plot(H2Q_C,ESQQ,'b')
plot(H2Q_C,EQ2SQ,'r')
plot(H2Q_C,(ESQQ+EQ2SQ)./2,'k')
hold off



%I can now also estimate the full protonated standard redox potentials in
%order to determine whether it aligns with calculations in the litrature
%Han et al.

ESQQ_protonated = ESQQ - 0.059.*log10(Kr1) ;

EQ2SQ_protonated = EQ2SQ - 0.059.*log10(Ka2.*Ka1./Kr1) ;

E3_protonated = (ESQQ_protonated + EQ2SQ_protonated)./2 ;

figure(6)
hold on
plot(H2Q_C,ESQQ_protonated,'b')
plot(H2Q_C,EQ2SQ_protonated,'r')
plot(H2Q_C,E3_protonated,'k')
hold off



%testing to see if oxidation peak for supposed semiquinone from
%Serpintinity aligns with the reducti0.on potential of SQ calculated here.

H = 10^(-5.6) ;

ESQQi = ESQQ - 0.059.*log10(Kr1) + 0.059.*log10( (Kr1 + H) ) ;

figure(7)
plot(H2Q_C,ESQQi)