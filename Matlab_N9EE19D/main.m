clear all
close all
clc

tic;% Début de la mesure du temps

Ka=0.1;
Vm=13;                      
Cn=22e-2;
Kc=5.7e-2;
Jmot=4.7e-5;%Inertie globale
pas=5e-3;
m=1.3;
masse = [ 0 5 10 15 20 25 30 35 40 45 50 55 ];
f=1.52e-4;
Cr=0;%Cn/4;
Jtot=Jmot+m*(pas/(2*pi))^2;

Tauf = Jtot / f;
disp(['Tauf = ', num2str(Tauf)]);

%%%%%%%% proportionelle %%%%%%%
xi=1;
K1=f^2*(2*pi)/(Kc*Ka*pas*Jtot*4*xi^2);
disp(['Correcteur proportionelle K1 = ', num2str(K1)]);


%%%%%%%%% méthode de Kessler l optimum symétrique %%%%%%%

Phimax=49*pi/180;

Kgain=Ka*Kc*pas/(2*pi*f);

a_Kessler = (1+sin(Phimax))/(1-sin(Phimax));
wo_Kessler = 1/(sqrt(a_Kessler)*Tauf);
Taui_Kessler = a_Kessler*Tauf;
Ki_Kessler = wo_Kessler^2 *Taui_Kessler*sqrt(1+(wo_Kessler*Tauf)^2)/(Kgain*sqrt(1+(wo_Kessler*a_Kessler*Tauf)^2));


%%%%%%%% PID ITAE %%%%%%%%%%

Temps5_ITAE=1;
wn_ITAE = 3/Temps5_ITAE;
N_ITAE=2.15*wn_ITAE-1/Tauf;
Kp_ITAE=(3.14*Tauf*wn_ITAE^2-N_ITAE)/Kgain;
KI_ITAE=2.7*wn_ITAE^3*Tauf/Kgain-Kp_ITAE;
Kd_ITAE=wn_ITAE^4*Tauf/Kgain-N_ITAE*KI_ITAE;



%%%%%%%%%%%%%%

% Initialisation des matrices pour stocker les résultats
positions = [];   % Matrice pour stocker les positions
commandes = [];   % Matrice pour stocker les commandes
consignes = [];   % Matrice pour stocker la consigne (même taille que les autres)
time_vector = []; % Vecteur temporel commun à toutes les simulations
Temps_de_reponse=zeros(1,length(masse));

% Choisir un colormap pour les courbes
colors = jet(length(masse)); % Colormap 'jet' qui fournit un dégradé de couleurs en arc-en-ciel

for k = 1 : length(masse)
    % Calcul du nouveau moment d'inertie
    Jnew = Jmot + (m + masse(k)) * (pas / (2 * pi))^2;

    % Exécution de la simulation
    out = sim('Sim_moteur.slx');

    % Extraction des données de la simulation
    position = out.position.Data;    % Extraction des données numériques de position
    time = out.position.Time;        % Extraction des données temporelles
    consigne = out.consigne.Data;    % Extraction des données de consigne
    commande = out.commande.Data;    % Extraction des données de commande
    
    % Stockage des résultats dans les matrices
    positions = [positions, position];  % Ajouter la position en colonne
    commandes = [commandes, commande];  % Ajouter la commande en colonne
    consignes = [consignes, consigne];  % Ajouter la consigne en colonne

    % Stocker le temps une seule fois (supposé identique pour toutes les simulations)
    if isempty(time_vector)
        time_vector = time;
    end
end

% Calcul de la valeur finale de la consigne
valeur_finale = consignes(end, 1);

% Détermination du seuil de réponse à 95% de la valeur finale
seuil = 0.95 * valeur_finale;

% Affichage des résultats en une seule fois après la boucle
figure(1)
% Tracé de la réponse en position pour chaque masse
subplot(2,1,1);
hold on;
for k = 1 : length(masse)
    % Calcul du temps de réponse pour la masse actuelle
    indice_reponse = find(positions(:,k) >= seuil, 1);
    Temps_de_reponse(1,k)=time_vector(indice_reponse);
    % Vérifier si l'indice de réponse n'est pas vide avant de l'utiliser
    if ~isempty(indice_reponse)
        % Ajouter le temps de réponse à la légende
        label = ['Masse = ' num2str(masse(k)) ' kg, Tr_{5%} = ' num2str(time_vector(indice_reponse), '%.2f') ' s'];
    else
        % Si l'indice est vide, on ne peut pas calculer le temps de réponse
        label = ['Masse = ' num2str(masse(k)) ' kg, Tr_{5%} non défini'];
    end
    
    % Tracé de chaque courbe avec une couleur différente
    plot(time_vector, positions(:,k), 'Color', colors(k, :), 'LineWidth', 1.5, 'DisplayName', label);
end

% Tracé de la courbe de consigne (assumée identique pour toutes les masses)
plot(time_vector, consignes(:,1), '-','LineWidth', 1.5, 'DisplayName', 'Consigne'); % Consigne identique donc on prend la première colonne
hold off;
xlabel('Temps (s)');
ylabel('Position (rad)');
title('Réponse en position pour différentes masses');
legend show;
grid on;

% Tracé de la commande pour chaque masse
subplot(2,1,2);
hold on;
for k = 1 : length(masse)
    % Tracé de chaque courbe de commande avec une couleur différente
    plot(time_vector, commandes(:,k), 'Color', colors(k, :), 'LineWidth', 1.5, 'DisplayName', ['Masse = ' num2str(masse(k)) ' kg']);
end
hold off;
xlabel('Temps (s)');
ylabel('Commande (V)');
title('Réponse en commande pour différentes masses');
legend show;
grid on;

figure(2)
plot(masse,Temps_de_reponse,'x-','LineWidth', 1.5);
xlabel('Masse (kg)');
ylabel('Temps de réponse à 5% (s)');
title('Temps de réponse à 5% pour différentes masses');
legend show;
grid on;

