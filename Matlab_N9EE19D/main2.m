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
masse = [ 0 20 40 ];
f=1.52e-4;
Cr=0;%Cn/4;
Jtot=Jmot+m*(pas/(2*pi))^2;
Te=1e-4;

Tauf = Jtot / f;
disp(['Tauf = ', num2str(Tauf)]);

Kgain=Ka*Kc*pas/(2*pi);

%% Discret %%

Fonction_Transfert_moteur = tf(1,[Tauf, 1, 0]);
Discret_Transfert_moteur = c2d(Fonction_Transfert_moteur,Te,'impulse');


%% espace d etat %%
A_etat=[0, 2*pi/pas ; 0, -1/Tauf];
B_etat=[0 ; Kgain*pas/(Jtot*2*pi)];
C_etat=[1, 0];

%%%% espace augmenté %%%

% Dimensions de chaque composant
[nA, mA] = size(A_etat); % Taille de A_etat (2x2)
[nC, mC] = size(C_etat); % Taille de C_etat (1x2)

% Définition des matrices du système augmenté
A_etat_augment = [A_etat, zeros(nA, 1); -C_etat, 0]; % A_etat et -C_etat doivent avoir le même nombre de colonnes
B_etat_augment = [B_etat; 0];                       % B_etat doit avoir le bon nombre de lignes
C_etat_augment = [C_etat, 0];                       % C_etat doit avoir le bon nombre de colonnes

% Définition des pôles désirés
P1 = -3;
P2 = P1*10;
P = [P1, P1, P2];

% Calcul du gain de retour d'état avec acker
K_ac = acker(A_etat_augment, B_etat_augment, P);

% Extraction des gains individuels
K1 = K_ac(1);
K2 = K_ac(2);
Kr = K_ac(3);
Kt=[K1, K2];
g_comp=Kr/P2;

g_annule=1/(-C_etat*inv(A_etat-B_etat*Kt)*B_etat);

% Affichage des résultats
disp(['K1 = ', num2str(K1)]);
disp(['K2 = ', num2str(K2)]);
disp(['Kr = ', num2str(Kr)]);
disp([' compensation du pole 3eme pôle, g_comp = ', num2str(g_comp)]);
disp([' Anulation de régulation, g_annule = ', num2str(g_annule)]);

%%% espace augmenté corrigé %%%
A_etat_augment_corrige = [A_etat-B_etat*Kt , -B_etat*Kr; -C_etat, 0];
B_etat_augment_corrige = [B_etat*g_comp; 1];
C_etat_augment_corrige = [1, 0, 0];



sys_augment=ss(A_etat_augment_corrige,B_etat_augment_corrige,C_etat_augment_corrige,0);




%%

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
    out = sim('Sim_moteur2.slx');

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


% Affichage des résultats en une seule fois après la boucle
figure(1)
% Tracé de la réponse en position pour chaque masse
subplot(2,1,1);
hold on;
for k = 1 : length(masse)
    % Calcul de la valeur finale de la consigne
    valeur_finale = positions(end, k);
    % Détermination du seuil de réponse à 95% de la valeur finale
    seuil = 0.95 * valeur_finale;
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

figure(3)

impulse(Fonction_Transfert_moteur,Discret_Transfert_moteur);


figure(4)
iopzmap(sys_augment);
hold on 
cercle(0,0,1);
hold off
title('Système 0');
grid on;
