%% Partie SIS/SIR
%Affichage en sortie :

%Figure 1 : état sur la trajectoire
%etat : xrobot,yrobot,xamer1,yamer1,xamer2,yamer2 (6 premieres courbes)
%bleu : état caché sur la trajectoire
%rouge : intervalle +/- 3 sigma autour de l'estimée de l'état
% 2 dernieres courbes : Neff (repartition des poids des particules) en echelle lineaire et log

%Figure 2: points pleins : état caché (position robot, position amer 1, position amer 2)
%points creux : estimées de l'etat, avec ellipses de confiance à plus ou moins 3*ecart type

clear all;
close all;
clc;

%% donnees probleme (partie modifiable)
plot_p = 0; %0 pour ne rien afficher, 1 pour voir la simulation/generation de donnees
Npart = 10; %nombre de particules
Ns = 0.3*Npart; %seuil de reechantillonnage pour SIR (si Neff < Ns, reechantillonnage)
b = 0; %mettre a 1 pour voir les particules ! (taille marqueur prop. au poids)
%(Ne pas en mettre trop, sinon affichage long...)
algo = 0; %choix de l'algo :
%0 : SIS avec importance = dyn a priori
%1 : SIR avec importance = dyn a priori
%2 : SIS avec importance prenant en compte la mesure
%3 : SIR avec importance prenant en compte la mesure

%% Simulation donnees et lancement du filtrage

[N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,Xreel] = simulationDonnees(plot_p);
X = nan(6,Npart,N);
W = nan(1,Npart,N);

if algo == 0 %algo SIS
    %k=0 :
    [x,w] = SIS(nan(6,Npart),nan(1,Npart),[],0); %[] : mesure vide, 0 : k 
    X(:,:,1) = x;
    W(:,:,1) = w;

    %k>0 :
    for k = 1:N-1
        [x,w] = SIS(X(:,:,k),W(:,:,k),Z(:,k),k); %x : (6,Part) w : (1,Part)
        X(:,:,k+1)=x;
        W(:,:,k+1)=w;
    end
    
elseif algo == 1 %algo SIR 
    %k=0 :
    [x,w] = SIR(nan(6,Npart),nan(1,Npart),[],0,1); %[] : mesure vide, 0 : k, 1 : pas resamp. 
    X(:,:,1) = x;
    W(:,:,1) = w;
    %k>0 :
    for k = 1:N-1
        [x,w] = SIR(X(:,:,k),W(:,:,k),Z(:,k),k,Ns); %x : (6,Part) w : (1,Part)
        X(:,:,k+1)=x;
        W(:,:,k+1)=w;
    end
    
elseif algo == 2 %algo SIS importance mesure
    %k=0 :
    [x,w] = SIS2(nan(6,Npart),nan(1,Npart),[],0); %[] : mesure vide, 0 : k 
    X(:,:,1) = x;
    W(:,:,1) = w;

    %k>0 :
    for k = 1:N-1
        [x,w] = SIS2(X(:,:,k),W(:,:,k),Z(:,k),k); %x : (6,Part) w : (1,Part)
        X(:,:,k+1)=x;
        W(:,:,k+1)=w;
    end
    
elseif algo == 3 %algo SIR importance mesure
    %k=0 :
    [x,w] = SIR2(nan(6,Npart),nan(1,Npart),[],0,1); %[] : mesure vide, 0 : k, 1 : pas resamp. 
    X(:,:,1) = x;
    W(:,:,1) = w;
    %k>0 :
    for k = 1:N-1
        [x,w] = SIR2(X(:,:,k),W(:,:,k),Z(:,k),k,Ns); %x : (6,Part) w : (1,Part)
        X(:,:,k+1)=x;
        W(:,:,k+1)=w;
    end
end

%% estimation des moments a posteriori :

Ex = zeros(6,N); %esperance etat
for k=1:N
    for i = 1:Npart
        Ex(:,k) = Ex(:,k)+X(:,i,k)*W(1,i,k);
    end
end

Px = zeros(6,6,N); %covariance etat
for k = 1:N
    for i = 1:Npart
        Px(:,:,k) = Px(:,:,k)+((X(:,i,k)-Ex(:,k))*(X(:,i,k)-Ex(:,k))')*W(1,i,k);
    end
end

%% Affichage des moments (et particules)

figure(1)
courbes(Xreel,Ex,Px,W);

figure(2)
plotEllipse(Xreel,Ex,Px,X,W,b);