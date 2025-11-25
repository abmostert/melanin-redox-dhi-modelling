close all
clear all

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
% HSQ - concentration of the protonated semiquinone molecule in Molar
% HQI - concentration of protonated quinone methide
% QI - concentration of deprotonated quinone methide
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

%The initial concentration of the melanin solution in Molar assumming a
%70/30 mixture of DHI/DHICA
initialC = 0.077 ;
%Setting the initial fraction of H2Q
initialH2Qfrac = 0.01 ;

%Determining the initial parameters for the modelling
CH2Q0 = initialC.*initialH2Qfrac  ; 
CQ0 = initialC.*(1-initialH2Qfrac)  ;
V0 = 10 ;
C = 1 ;

%Setting up the variable V on which the modelling is going to depend
Vmin = 0.00001 ;    
Vstep = 0.001 ;
Vmax = 2.00001 ;
Vlength = round((Vmax-Vmin)/Vstep) ;

%Since there are 3 variables to calculate, I need to supply initial conditions for them. 
%The initial conditions are set up as a vector
x0 = [0 0 0];
%x0(1) is the guess for the potential
x0(1) = 0.5 ;
%x0(1) is the guess for the pH
x0(2) = 6 ;
%x0(1) is the guess for the quinone concentration
x0(3) = -log10(CQ0);


%Creating a series of empty vectors in which the calculated results will be
%dropped
logOHsave =  zeros(1,(Vlength));   
logH2Qsave =  zeros(1,(Vlength));
logHQsave =  zeros(1,(Vlength));
logQ2save =  zeros(1,(Vlength));
logQsave =  zeros(1,(Vlength));
logSQsave =  zeros(1,(Vlength));
logHSQsave =  zeros(1,(Vlength));
logHQIsave =  zeros(1,(Vlength));
logQIsave =  zeros(1,(Vlength));
Esave = zeros(1,(Vlength)) ;
pHsave = zeros(1,(Vlength)) ;
pQsave = zeros(1,(Vlength)) ;


%The actual modelling of the Chio data is accomplished with a for loop
%running through the total voume being added.
for V=Vmin:Vstep:Vmax
    
    %Setting up the volume being added to the system for the given step.
    Vlength = (Vmax-Vmin)/Vstep ;
    Vcurrentstepraw = (V-Vmin)/Vstep +1 ;
    Vcurrentstep = round(Vcurrentstepraw) ;
    
    %setting up the least squares algorith,levenberg-marquardt for determining the unknowns  
    options=optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',5000,'MaxFunEvals',5000,'TolFun',1.0e-10,'TolX',1.0e-10);
    %calling on the fsolve function and function file to model.
    [x,fval,exitflag,output] = fsolve(@Model_Function_Chio_2_Final,x0,options) ;
    
    %since the calculation is complete, I'm writing the new set of initial
    %guesses for the next round of calculations.
    x0=x;               
    
    %writing the results from the calculations to their respective vectors.
    logOHsave(Vcurrentstep) = logOH ;   
    logH2Qsave(Vcurrentstep) = logH2Q;
    logHQsave(Vcurrentstep) = logHQ;
    logQ2save(Vcurrentstep) = logQ2; 
    logQsave(Vcurrentstep) = logQ; 
    logSQsave(Vcurrentstep) = logSQ; 
    logHSQsave(Vcurrentstep) = logHSQ; 
    logSQsave(Vcurrentstep) = logSQ; 
    logHQIsave(Vcurrentstep) = logHQI; 
    logQIsave(Vcurrentstep) = logQI; 
    Esave(Vcurrentstep) =  x(1);
    pHsave(Vcurrentstep) = x(2) ;
    pQsave(Vcurrentstep) = x(3) ;
   
end

%Setting up the volume in order to determine the total SQ as a function of
%pH
Vvector = linspace(1e-5,2.00001,2001) ;
%determining the total number of SQ in number of particles per gram
SQtotal = ((((10.^(-logSQsave)).*((V0+Vvector)/1000))).*6.02e23)/0.125 + ((((10.^(-logHSQsave)).*((V0+Vvector)/1000))).*6.02e23)/0.125  ;

%plotting the equivalent curve as in Chio et al
figure(1)
plot(pHsave,SQtotal) ;

%plotting the concentration vs pH behvaiour of all chemical specicies
figure(2)
hold on
plot(pHsave,logSQsave,'b')
plot(pHsave,logH2Qsave,'r')
plot(pHsave,logHQsave,'k')
plot(pHsave,logQ2save,'g')
plot(pHsave,logQsave,'y')
plot(pHsave,logHSQsave,'m')
plot(pHsave,logQIsave,'b')
plot(pHsave,logHQIsave,'g')
hold off

%plotting the potential vs pH
figure(3)
plot(pHsave,Esave)

%plotting the fraction of semiquinone vs quinone like moeities.
figure(4)
plot(pHsave,(10.^(-logSQsave)+10.^(-logHSQsave))./(10.^(-logH2Qsave) + 10.^(-logHQsave) + 10.^(-logQ2save) + 10.^(-logQsave) + 10.^(-logQIsave) + 10.^(-logHQIsave) ))

%plotting the fraction of semiquinone vs quinone imine.
figure(5)
plot(pHsave,(10.^(-logSQsave)+10.^(-logHSQsave))./( 10.^(-logQIsave) + 10.^(-logHQIsave) ))

%plotting a g value change. I'm assuming that HSQ has a g value of 2.0034,
%SQ a g value of 2.0043

g_average = 2.0034.*( (10.^(-logHSQsave))./((10.^(-logHSQsave)) + (10.^(-logSQsave)))) + 2.0043.*( (10.^(-logSQsave))./((10.^(-logHSQsave)) + (10.^(-logSQsave)))) ;

figure(6)
plot(pHsave,g_average)


%final data matrix for 

datamatrix = [flip(pHsave') flip(SQtotal') flip(logSQsave') flip(logHSQsave') flip(logH2Qsave') flip(logHQsave') flip(logQ2save') flip(logQsave') flip(logHQIsave') flip(logQIsave')  flip(Esave')] ;

%datamatrix = [(pHsave') (SQtotal') (logSQsave') (logHSQsave') (logH2Qsave') (logHQsave') (logQ2save') (logQsave') (logHQIsave') (logQIsave') (Esave')] ;