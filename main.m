clear all
close all
clc

F_MLI=10e3;                    %fréquence de MLI
T_MLI=1/F_MLI;
V_MLI=540;                      %tension MLI à partir d une source continue de 540 v.
p=2;                            %deux paires de pôles
Pn=4000;                        %Puissance nominale
Vn=380;                         %Tension nominale
In=8.3;                         %Fréquence nominale
Fn=50;                          %Courant nominal
ws=2*pi*Fn;
Rs=1.6;                         %Résistance statorique 
Rr=1.15;                        %Résistance rotorique
Ls=0.182;                       %Inductance cyclique stator 
Lr=0.194;                       %Inductance cyclique rotor
Msr=0.182;                      %Mutuelle inductance cyclique
sigma=1-Msr^2/(Ls*Lr);            %Coefficient de dispersion
J=0.025;                        %Inertie globale
frot=0.01;                      %Frottement visqueux
T_meca=J/frot;
T_charge=0;
Trau=T_meca*1e-3;
Rsr=Rs+Rr*(Msr/Lr)^2;

%%%%%%%%%%% moteur %%%%%%%%%%%%%%
abc_alpha_beta = sqrt(2/3) * [1, -1/2, -1/2; 
                             0, sqrt(3)/2, -sqrt(3)/2; 
                             1/sqrt(2), 1/sqrt(2), 1/sqrt(2)];

Inv_abc_alpha_beta = inv(abc_alpha_beta);


Phi_Courant= [Ls,  0,   Msr,  0; 
                0,   Ls,  0,   Msr; 
                Msr, 0,   Lr,  0; 
                0,   Msr, 0,   Lr];

Inv_Phi_Courant = inv(Phi_Courant);

Flux_1 = 1/(Msr^2 - Ls * Lr )  * [Rs*Lr, 0, -Rs*Msr, 0; 
                                  0, Rs*Lr, 0, -Rs*Msr; 
                                  -Rr*Msr, 0, Rr*Ls, 0; 
                                  0, -Rr*Msr,0, Rr*Ls];

Flux_2 = 1/(Msr^2 - Ls * Lr )  * [0, 0, 0, 0;
                                  0, 0, 0, 0; 
                                  0, 0, 0, -(Msr^2-Ls*Lr); 
                                  0, 0, Msr^2-Ls*Lr,0];


%%%%%%%%%%%%% MLI/onduleur ********

onduleur= V_MLI/3*[2,-1,-1;
            -1,2,-1;
            -1,-1,2];
    

%%%%%
xi_reg=1;
wo_reg_courant=F_MLI/sqrt(100);
ki= 2*wo_reg_courant*sigma*Ls*xi_reg-Rsr;
Ti=ki/(wo_reg_courant^2*Ls*sigma);

%%%%%
xi_reg_flux=1;
wo_reg_flux=wo_reg_courant/10;
K_reg_flux=2*xi_reg_flux*wo_reg_flux*Lr/Rr-1/Msr;
T_reg_flux=Msr*K_reg_flux*Lr/(Rr*wo_reg_flux^2);


%%%%%%%

%%%% pour xi = 1, wo*t5% = 4.5 et on veux t5% = 0.4s
xi_reg_vitesse=1;
wo_reg_vitesse=45/4;
K_reg_vitesse=(2*xi_reg_vitesse*wo_reg_vitesse*T_meca-1)*frot;
T_reg_vitesse=K_reg_vitesse/(T_meca*frot*wo_reg_vitesse^2);

