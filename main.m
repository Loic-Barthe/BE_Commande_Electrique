clear all
close all
clc

F_MLI=10e3;                    %fréquence de MLI
T_MLI=1/F_MLI;
V_MLI=540;                      %tension MLI à partir d une source continue de 540 v.
p=2;                            %deux paires de pôles
Pn=4000;                        %Puissance nominale
Vn=220;                         %Tension nominale
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

%%%%%%% Wreau ***********
Kpi=ws*20;
Tin=sqrt(10)/Kpi;

%%%%%


%%%%%
xi_reg=1;
wo_reg_courant=F_MLI/10;

ki= 2*wo_reg_courant*sigma*Ls*xi_reg-Rsr;
Ti=ki/(wo_reg_courant^2*Ls*sigma);

kpd= ki;
kid= ki/Ti;
kpq=ki;
kiq=ki/Ti;

%%%%%
xi_reg_flux=1;
wo_reg_flux=wo_reg_courant/10;
K_reg_flux=(2*xi_reg_flux*wo_reg_flux*Lr/Rr-1)/Msr;
T_reg_flux=Msr*K_reg_flux*Lr/(Rr*wo_reg_flux^2);

kp_reg_flux=K_reg_flux;
ki_reg_flux=K_reg_flux/T_reg_flux;

%%%%%%%
xi_reg_vitesse=1;
wo_reg_vitesse= wo_reg_flux/20; 
K_reg_vitesse=(2*xi_reg_vitesse*wo_reg_vitesse*T_meca-1)*frot;
T_reg_vitesse=K_reg_vitesse/(T_meca*frot*wo_reg_vitesse^2);

ki_reg_vitesse=K_reg_vitesse;
kp_reg_vitesse=K_reg_vitesse/T_reg_vitesse;



% %xi_reg=1;
% wo_reg_courant=F_MLI/10;
% 
% kpd= Ls*sigma*wo_reg_courant;
% kid= kpd*wo_reg_courant/sqrt(10);
% 
% kpq=Ls*sigma*wo_reg_courant;
% kiq=kpq*wo_reg_courant/sqrt(10);
% 
% 
% %%%%%
% wo_reg_flux=wo_reg_courant/10;
% 
% % ki_reg_flux=wo_reg_flux/Msr;
% % kp_reg_flux=Lr/(Rr*Msr)*wo_reg_flux;
% 
% kp_reg_flux=Lr*wo_reg_flux/(Msr*Rr);
% ki_reg_flux=kp_reg_flux*wo_reg_flux/sqrt(10);
% 
% %%%%
% xi_reg_vitesse=10;
% wo_reg_vitesse= wo_reg_flux/20;
% % kp_reg_vitesse=J*wo_reg_vitesse;
% % ki_reg_vitesse=kp_reg_vitesse*wo_reg_vitesse/sqrt(10);
% ki_reg_vitesse=J*wo_reg_vitesse^2;
% kp_reg_vitesse=2*xi_reg_vitesse*J*wo_reg_vitesse;


%%%%% observateur

A_observateur=[-Rsr/(sigma*Ls), 0, Msr*Rr/(sigma*Ls*Lr^2), 0; 
                0,-Rsr/(sigma*Ls), 0, Msr*Rr/(sigma*Ls*Lr^2); 
                Msr*Rr/Lr, 0, -Rr/Lr, 0; 
                0, Msr*Rr/Lr,0, -Rr/Lr];

A_ob_w=[0, 0, 0, Msr*Rr/(sigma*Ls*Lr); 
        0, 0, -Msr*Rr/(sigma*Ls*Lr), 0; 
        0, 0, 0, -1; 
        0, 0 1, 0];

B_observateur=[1/(sigma*Ls),0;
            0, 1/(sigma*Ls);
            0, 0;
            0, 0];

C_observateur=[1,0,0,0;
                0,1,0,0];

%%%placement des poles
Pole=10;

K_observateur=[-Rs/Ls-2*Pole, 0;
                0, -Rs/Ls-2*Pole;
                -Ls*Pole^2, 0;
                0, -Ls*Pole^2];



%% appel de la simulation %%
out = sim('MAS.slx');

% Extraction des paramètres depuis la structure de simulation
consigne = out.consigne.Data; % Extraction des données du signal
time = out.consigne.Time;        % Extraction du vecteur temps

consigne_triangle_MLI=out.consigne_triangle_MLI.Data;
Alpha_MLI=out.Alpha_MLI.Data;
commande_Modulateur=out.commande_Modulateur.Data;
consigne_MLI=out.consigne_MLI.Data;

w_elec = out.w_elec.Data;
w_mec = out.w_mec.Data;
Puis = out.Puis.Data;
glis = out.glis.Data;
Tem=out.Tem.Data;

ISabc=out.ISabc.Data;
IRabc=out.IRabc.Data;
ISAlpha_Beta=out.ISAlpha_Beta.Data;
VAlpha_Beta_O=out.VAlpha_Beta_O.Data;
Phi=out.Phi.Data;

wrau_result=out.wrau_result.Data;
Phir_result=out.Phir_result.Data;
rau_result=out.rau_result.Data;

%% Create a new figure
% Create a new figure for the plots
figure(1);

% First subplot: Consignes and consigne_triangle_MLI
subplot(2, 2, 1); % 2x2 grid, first subplot
plot(time, consigne(:, 1), 'r', 'LineWidth', 1.5); % Consigne 1 in red
hold on;
plot(time, consigne(:, 2), 'b', 'LineWidth', 1.5); % Consigne 2 in blue
plot(time, consigne(:, 3), 'g', 'LineWidth', 1.5); % Consigne 3 in green
plot(time, consigne_triangle_MLI, 'm-o', 'LineWidth', 1.5); % Consigne triangle in yellow
hold off;

% Add labels, title, and legend
title('Consignes et Consigne Triangle MLI');
xlabel('Temps (s)');
ylabel('Valeur');
legend('Consigne 1', 'Consigne 2', 'Consigne 3', 'Consigne Triangle MLI');
grid on;

% Second subplot: Alpha_MLI
subplot(2, 2, 2); % 2x2 grid, second subplot
plot(time, Alpha_MLI(:, 1), 'r', 'LineWidth', 1.5); % Alpha_MLI 1 in red
hold on;
plot(time, Alpha_MLI(:, 2), 'b', 'LineWidth', 1.5); % Alpha_MLI 2 in blue
plot(time, Alpha_MLI(:, 3), 'g', 'LineWidth', 1.5); % Alpha_MLI 3 in green
hold off;

% Add labels, title, and legend
title('Alpha MLI');
xlabel('Temps (s)');
ylabel('Valeur');
legend('Alpha MLI 1', 'Alpha MLI 2', 'Alpha MLI 3');
grid on;

% Third subplot: Commande Modulateur
subplot(2, 2, 3); % 2x2 grid, third subplot
plot(time, commande_Modulateur(:, 1), 'r', 'LineWidth', 1.5); % Commande Modulateur 1 in red
hold on;
plot(time, commande_Modulateur(:, 2), 'b', 'LineWidth', 1.5); % Commande Modulateur 2 in blue
plot(time, commande_Modulateur(:, 3), 'g', 'LineWidth', 1.5); % Commande Modulateur 3 in green
hold off;

% Add labels, title, and legend
title('Commande Modulateur');
xlabel('Temps (s)');
ylabel('Valeur');
legend('Commande Modulateur 1', 'Commande Modulateur 2', 'Commande Modulateur 3');
grid on;

% Fourth subplot: Consigne MLI
subplot(2, 2, 4); % 2x2 grid, fourth subplot
plot(time, consigne_MLI(:, 1), 'r', 'LineWidth', 1.5); % Consigne MLI 1 in red
hold on;
plot(time, consigne_MLI(:, 2), 'b', 'LineWidth', 1.5); % Consigne MLI 2 in blue
plot(time, consigne_MLI(:, 3), 'g', 'LineWidth', 1.5); % Consigne MLI 3 in green
hold off;

% Add labels, title, and legend
title('Consigne MLI');
xlabel('Temps (s)');
ylabel('Valeur');
legend('Consigne MLI 1', 'Consigne MLI 2', 'Consigne MLI 3');
grid on;

% Set a global title for the figure
sgtitle('Analyse des Consignes et Modulation MLI');

% Create a new figure for the plot
figure(2);

% Plot the first element of consigne and consigne_MLI
plot(time, consigne_MLI(:, 1), 'b', 'LineWidth', 1.5); % First element of consigne_MLI in blue
hold on;
plot(time, consigne(:, 1), 'r', 'LineWidth', 1.5); % First element of consigne in red
hold off;

% Add labels, title, and legend
title('Premier Élément de la Consigne et Consigne MLI');
xlabel('Temps (s)');
ylabel('Valeur');
legend('Consigne', 'Consigne après MLI');
grid on;

%%%%
figure(3);

% First subplot: w_elec and w_mec
subplot(2, 2, 1); % 2x2 grid, first subplot
plot(time, w_elec, 'r', 'LineWidth', 1.5); % Plot w_elec in red
hold on;
plot(time, w_mec, 'b', 'LineWidth', 1.5); % Plot w_mec in blue
hold off;
title('Vitesse Électrique et Mécanique');
xlabel('Temps (s)');
ylabel('Vitesse (tr/min)');
legend('\Omega_{elec}', '\Omega_{mec}');
grid on;

% Second subplot: Puissance
subplot(2, 2, 2); % 2x2 grid, second subplot
plot(time, Puis, 'g', 'LineWidth', 1.5); % Plot Puissance in green
title('Puissance');
xlabel('Temps (s)');
ylabel('Puissance (W)');
grid on;

% Third subplot: Glissement
subplot(2, 2, 3); % 2x2 grid, third subplot
plot(time, glis, 'm', 'LineWidth', 1.5); % Plot glissement in magenta
title('Glissement');
xlabel('Temps (s)');
ylabel('Glissement');
grid on;

% Fourth subplot: Tem (Couple)
subplot(2, 2, 4); % 2x2 grid, fourth subplot
plot(time, Tem, 'c', 'LineWidth', 1.5); % Plot Tem in cyan
title('Couple Électromagnétique (Tem)');
xlabel('Temps (s)');
ylabel('Couple (N·m)');
grid on;

% Adjust layout
sgtitle('Analyse des Signaux'); % Title for the entire figure


% Create a new figure for the currents and flux
figure(4);

% First subplot: Stator currents (ISabc)
subplot(2, 2, 1); % 2x2 grid, first subplot
plot(time, ISabc(:, 1), 'r', 'LineWidth', 1.5); % ISa in red
hold on;
plot(time, ISabc(:, 2), 'g', 'LineWidth', 1.5); % ISb in green
plot(time, ISabc(:, 3), 'b', 'LineWidth', 1.5); % ISc in blue
hold off;
title('Courants Statoriques');
xlabel('Temps (s)');
ylabel('Courant (A)');
legend('ISa', 'ISb', 'ISc');
grid on;

% Second subplot: Rotor currents (IRabc)
subplot(2, 2, 2); % 2x2 grid, second subplot
plot(time, IRabc(:, 1), 'm', 'LineWidth', 1.5); % IRa in magenta
hold on;
plot(time, IRabc(:, 2), 'c', 'LineWidth', 1.5); % IRb in cyan
plot(time, IRabc(:, 3), 'k', 'LineWidth', 1.5); % IRc in black
hold off;
title('Courants Rotoriques');
xlabel('Temps (s)');
ylabel('Courant (A)');
legend('IRa', 'IRb', 'IRc');
grid on;

% Third subplot: Currents in the alpha-beta frame (ISAlpha_Beta)
subplot(2, 2, 3); % 2x2 grid, third subplot
plot(time, ISAlpha_Beta(:, 1), 'r', 'LineWidth', 1.5); % IS_alpha in red
hold on;
plot(time, ISAlpha_Beta(:, 2), 'b', 'LineWidth', 1.5); % IS_beta in blue
hold off;
title('Courants Alpha-Beta');
xlabel('Temps (s)');
ylabel('Courant (A)');
legend('IS\_alpha', 'IS\_beta');
grid on;

% Fourth subplot: Flux components (Phi)
subplot(2, 2, 4); % 2x2 grid, fourth subplot
plot(time, Phi(:, 1), 'r', 'LineWidth', 1.5); % phiS_alpha in red
hold on;
plot(time, Phi(:, 2), 'b', 'LineWidth', 1.5); % phiS_beta in blue
plot(time, Phi(:, 3), 'g', 'LineWidth', 1.5); % phiR_alpha in green
plot(time, Phi(:, 4), 'm', 'LineWidth', 1.5); % phiR_beta in magenta
hold off;
title('Composantes du Flux');
xlabel('Temps (s)');
ylabel('Flux (Wb)');
legend('phiS\_alpha', 'phiS\_beta', 'phiR\_alpha', 'phiR\_beta');
grid on;

% Adjust layout for the first figure
sgtitle('Analyse des Courants et du Flux'); % Title for the entire figure

% Create a separate figure for the voltages (VAlpha_Beta_O)
figure(5);
plot(time, VAlpha_Beta_O(:, 1), 'r', 'LineWidth', 1.5); % V_alpha in red
hold on;
plot(time, VAlpha_Beta_O(:, 2), 'b', 'LineWidth', 1.5); % V_beta in blue
plot(time, VAlpha_Beta_O(:, 3), 'g', 'LineWidth', 1.5); % V_0 in green
hold off;
title('Tensions Alpha-Beta-0');
xlabel('Temps (s)');
ylabel('Tension (V)');
legend('V\_alpha', 'V\_beta', 'V\_0');
grid on;

figure(6); % Create a new figure
% Plot the third curve
subplot(3, 1, 1); % Create the third subplot in a 3x1 grid
plot(time,rau_result, 'g'); % Plot rau_result in green
title('rau\_result'); % Title for the third subplot
xlabel('Time'); % X-axis label
ylabel('rau\_result'); % Y-axis label
grid on; % Enable grid

% Plot the first curve
subplot(3, 1, 2); % Create the first subplot in a 3x1 grid
plot(time,wrau_result, 'b'); % Plot wrau_result in blue
title('wrau\_result'); % Title for the first subplot
xlabel('Time'); % X-axis label
ylabel('wrau\_result'); % Y-axis label
grid on; % Enable grid

% Plot the second curve
subplot(3, 1, 3); % Create the second subplot in a 3x1 grid
plot(time,Phir_result, 'r'); % Plot Phir_result in red
title('Phir\_result'); % Title for the second subplot
xlabel('Time'); % X-axis label
ylabel('Phir\_result'); % Y-axis label
grid on; % Enable grid

% Adjust figure layout for better visualization
sgtitle('Results for wrau, Phir, and rau'); % Overall title for the figure

