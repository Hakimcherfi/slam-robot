clear all;
close all;
clc;

%Affichage :
%Figure 1 : etat sur la trajectoire
%etat : xrobot,yrobot,xamer1,yamer1,xamer2,yamer2 (6 premieres courbes)
%bleu : vrai etat sur la trajectoire
%rouge : intervalle +/- 3 sigma autour de l'estimée de l'état
% 2 dernieres courbes : Neff (repartition poids) echelle lineaire et log
%Figure 2: points pleins : vrai etat (position robot, position amer 1,
%position amer 2
%points creux : estimées etat, avec ellipses de confiance à 3*ecart type

%donnees
plot_p = 0; %0 pour ne rien afficher, 1 pour voir la simulation/generation de donnees
Part = 100; %nombre de particules
Ns = 0.3*Part; %pour SIR, seuil du resampling
b = 0; %mettre a 1 pour voir les particules ! (Ne pas en mettre trop, sinon long...)
algo = 1; %choix de l'algo. 0 : SIR, autre : SIS

%debut algo
[N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,Xreel] = simulationDonnees(plot_p);
X = nan(6,Part,N);
W = nan(1,Part,N);
if algo == 0 %algo SIR
    %k=0 :
    [x,w] = SIR(nan(6,Part),nan(1,Part),[],0,1); %[] : mesure vide, 0 : k, 1 : pas resamp. 
    X(:,:,1) = x;
    W(:,:,1) = w;
    %k>0 :
    for k = 1:N-1
        [x,w] = SIR(X(:,:,k),W(:,:,k),Z(:,k),k,Ns); %x : (6,Part) w : (1,Part)
        X(:,:,k+1)=x;
        W(:,:,k+1)=w;
    end
    
else %algo SIS
    %k=0 :
    [x,w] = SIS(nan(6,Part),nan(1,Part),[],0); %[] : mesure vide, 0 : k 
    X(:,:,1) = x;
    W(:,:,1) = w;

    %k>0 :
    for k = 1:N-1
        [x,w] = SIS(X(:,:,k),W(:,:,k),Z(:,k),k); %x : (6,Part) w : (1,Part)
        X(:,:,k+1)=x;
        W(:,:,k+1)=w;
    end
    
end

%estimation des moments :
Ex = zeros(6,N); %esperance etat
for k=1:N
    for i = 1:Part
        Ex(:,k) = Ex(:,k)+X(:,i,k)*W(1,i,k);
    end
end

Px = zeros(6,6,N); %covariance etat
for k = 1:N
    for i = 1:Part
        Px(:,:,k) = Px(:,:,k)+((X(:,i,k)-Ex(:,k))*(X(:,i,k)-Ex(:,k))')*W(1,i,k);
    end
end

figure(1)
courbes(Xreel,Ex,Px,W);

figure(2)
plotEllipse(Xreel,Ex,Px,X,W,b);