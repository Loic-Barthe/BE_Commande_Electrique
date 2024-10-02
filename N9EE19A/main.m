clear all
close all
clc

%%%%%%%% Paramètre du système %%%%%%%%
A1=28e-4;
A2=32e-4;
A3=28e-4;
A4=32e-4;
a1=0.071e-4;
a2=0.057e-4;
a3=0.071e-4;
a4=0.057e-4;
g=9.81;
gamma1=0.7;
gamma2=0.6;

k1=3.33e-6;
k2=3.35e-6;

%%%%%%%%% conditions initale %%%%%%%%

V01=3;
V02=3;

Matrix_V0=[V01;V02];

%%%%%%% équilibre %%%%%%


H01= (gamma1*k1*V01+(1-gamma2)*k2*V02)/(a1*sqrt(2*g));
H02= (gamma2*k2*V02+(1-gamma1)*k1*V01)/(a2*sqrt(2*g));
H03= (1-gamma2)*k2*V02/(a3*sqrt(2*g));
H04= (1-gamma1)*k1*V01/(a4*sqrt(2*g));

Matrix_H0=[H01^2;H02^2;H03^2;H04^2];

%%%%%% Coef pour les matrice %%%%%%%%%

T1=A1/a1*sqrt(2/g)*H01;
T2=A2/a2*sqrt(2/g)*H02;
T3=A3/a3*sqrt(2/g)*H03;
T4=A4/a4*sqrt(2/g)*H04;

%%%%%% Matrice %%%%%%%%%

Matrix_dX_X =[-1/T1, 0, -A3/(A1*T3), 0; 
                0, -1/T2, 0, -A4/(A2*T4); 
                0, 0, -1/T3, 0; 
                0, 0, 0, -1/T4 ];

Matrix_dX_U = [gamma1*k1/A1, 0; 
                0, gamma2*k2/A2; 
                0, (1-gamma2)*k2/A3; 
                (1-gamma1)*k1/A4, 0 ];


%%%% PI %%%%%

K11=T1*gamma1*k1/A1;
K22=T2*gamma2*k2/A2;


xi=0.7;
wo1=100;
wo2=200;

Tc1= 2*xi/wo1-1/(T1*wo1^2);
kc1= (2*xi*wo1*T1-1)/K11;
Tc2= 2*xi/wo2-1/(T2*wo2^2);
kc2= (2*xi*wo1*T2-1)/K22;



