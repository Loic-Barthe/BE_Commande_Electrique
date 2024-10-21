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


xi=1;
wo1=9;
wo2=9;

Tc1= 2*xi/wo1-1/(T1*wo1^2);
kc1= (2*xi*wo1*T1-1)/K11;
Tc2= 2*xi/wo2-1/(T2*wo2^2);
kc2= (2*xi*wo1*T2-1)/K22;



%%%%%% commande %%%%%%

Com_H1=1.3*H01^2;
delay_Com_H1=1.5;

Com_H2=1.3*H02^2;
delay_Com_H2=0.5;

Com_H1_sched=[ 1 2 3 2 4 5 2 6 2 ]*H01^2;
Com_H2_sched=[ 1 2 3 2 4 5 2 6 2 ]*H02^2;
T_plage=2;

%% gain scheduling%%

%[ 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 ]
%[ 1 2 3 4 5 6 ]
H1_sched=[ 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 ]*H01^2;
H2_sched=[ 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 ]*H02^2;

% Préallocation des matrices de sortie
num_H1 = length(H1_sched);
num_H2 = length(H2_sched);
Tc1new = zeros(num_H1, num_H2);
kc1new = zeros(num_H1, num_H2);
Tc2new = zeros(num_H1, num_H2);
kc2new = zeros(num_H1, num_H2);

% Gain scheduling
for i = 1:num_H1
    for j = 1:num_H2
        % Calcul des nouvelles valeurs pour T1 et T2
        T1new = A1 / a1 * sqrt(2 / g) * H1_sched(i);
        T2new = A2 / a2 * sqrt(2 / g) * H2_sched(j);

        % Calcul des nouveaux gains K11 et K22
        K11new = T1new * gamma1 * k1 / A1;
        K22new = T2new * gamma2 * k2 / A2;

        % Calcul des nouveaux paramètres Tc1, kc1, Tc2, kc2
        Tc1new(i, j) = 2 * xi / wo1 - 1 / (T1new * wo1^2);
        kc1new(i, j) = (2 * xi * wo1 * T1new - 1) / K11new;
        Tc2new(i, j) = 2 * xi / wo2 - 1 / (T2new * wo2^2);
        kc2new(i, j) = (2 * xi * wo2 * T2new - 1) / K22new;
    end
end

disp('Correcteur proportionnel kc1new :');
disp(kc1new);
disp('Correcteur proportionnel Tc1new :');
disp(Tc1new);
disp('Correcteur proportionnel kc2new :');
disp(kc2new);
disp('Correcteur proportionnel Tc2new :');
disp(Tc2new);


%% appel de la simulation %%
out = sim('Sim3.slx');

% Extraction des paramètres depuis la structure de simulation
consigne_H1 = out.consigne_H1.Data; % Extraction des données du signal
consigne_H2 = out.consigne_H2.Data;
time = out.consigne_H1.Time;        % Extraction du vecteur temps

h_PI_lineaire = out.h_PI_lineaire.Data;
dh_PI_lineaire = out.dh_PI_lineaire.Data;
h_PI_non_lineaire = out.h_PI_non_lineaire.Data;
dh_PI_non_lineaire = out.dh_PI_non_lineaire.Data;

%% Tracer les résultats %%
colors = [
    1 0 1;        % Magenta
    1 0 0;        % Rouge
    0.6 0.1 0.4   % Rose foncé
    0 0 1;        % Bleu

];

figure(1)

% Tracé des courbes de hauteur en boucle ouverte linéaire
subplot(2,2,1)
hold on;
for k = 1:4  % Débuter la boucle à 1, car les indices commencent à 1 en MATLAB
    % Tracé de chaque courbe avec une couleur différente
    plot(time, h_PI_lineaire(:,k), 'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', ['h' num2str(k)]);
end
% Tracé des consignes H1 et H2
plot(time, consigne_H1, '--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Consigne H1');
plot(time, consigne_H2, '-.', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Consigne H2');
hold off
xlabel('Temps (s)');
ylabel('Hauteur (m)');
title('Hauteur en boucle PI linéaire');
legend show;
grid on;

subplot(2,2,3)
hold on;
for k = 1:4  % Débuter la boucle à 1, car les indices commencent à 1 en MATLAB
    % Tracé de chaque courbe avec une couleur différente
    plot(time, dh_PI_lineaire(:,k), 'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', ['dh/dt' num2str(k)]);
end
hold off
xlabel('Temps (s)');
ylabel('dH/dt (m.s^{-1})');
title('dH/dt linéaire');
legend show;
grid on;

subplot(2,2,2)
hold on;
for k = 1:4  % Débuter la boucle à 1, car les indices commencent à 1 en MATLAB
    % Tracé de chaque courbe avec une couleur différente
    plot(time, h_PI_non_lineaire(:,k), 'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', ['h' num2str(k)]);
end
% Tracé des consignes H1 et H2
plot(time, consigne_H1, '--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Consigne H1');
plot(time, consigne_H2, '-.', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Consigne H2');
hold off
xlabel('Temps (s)');
ylabel('Hauteur (m)');
title('Hauteur en boucle PI non linéaire');
legend show;
grid on;

subplot(2,2,4)
hold on;
for k = 1:4  % Débuter la boucle à 1, car les indices commencent à 1 en MATLAB
    % Tracé de chaque courbe avec une couleur différente
    plot(time, dh_PI_non_lineaire(:,k), 'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', ['dh/dt' num2str(k)]);
end
hold off
xlabel('Temps (s)');
ylabel('dH/dt (m.s^{-1})');
title('dH/dt non linéaire');
legend show;
grid on;





