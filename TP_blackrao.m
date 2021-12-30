clear all;
close all;
clc;

%Affichage :
%Figure 1 : état sur la trajectoire
%etat : xrobot,yrobot,xamer1,yamer1,xamer2,yamer2 (6 premieres courbes)
%bleu : vrai etat sur la trajectoire
%rouge : intervalle +/- 3 sigma autour de l'estimée de l'état
% 2 dernieres courbes : Neff (repartition poids) echelle lineaire et log
%Figure 2: points pleins : vrai etat (position robot, position amer 1,
%position amer 2
%points creux : estimées etat, avec ellipses de confiance à 3*ecart type

%% donnees
plot_p = 0; %0 pour ne rien afficher, 1 pour voir la simulation/generation de donnees
Part = 100; %nombre de particules
%Ns = 0.3*Part; %pour SIR, seuil du resampling
%b = 0; %mettre a 1 pour voir les particules ! (Ne pas en mettre trop, sinon long...)

[N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,Xreel] = simulationDonnees(plot_p);

%% definition structures de donnees
w = nan(1,Part,N);
r = nan(2,Part,N);
m1 = nan(2,Part,N);
m2 = nan(2,Part,N);
P1 = nan(2,2,Part,N);
P2 = nan(2,2,Part,N);

%% initialisation
[wp,rp,m1p,m2p,P1p,P2p] = blackrao(nan(2,Part),nan(2,Part),nan(2,Part),nan(2,2,Part),nan(2,2,Part),[],[],0);
w(:,:,1) = wp;
r(:,:,1) = rp;
m1(:,:,1) = m1p;
m2(:,:,1) = m2p;
P1(:,:,:,1) = P1p;
P2(:,:,:,1) = P2p;

%% filtre
for k = 1:N-1
    [wp,rp,m1p,m2p,P1p,P2p] = blackrao(r(:,:,k),m1(:,:,k),m2(:,:,k),P1(:,:,:,k),P2(:,:,:,k),w(:,:,k),Z(:,k),k);
    w(:,:,k+1) = wp;
    r(:,:,k+1) = rp;
    m1(:,:,k+1) = m1p;
    m2(:,:,k+1) = m2p;
    P1(:,:,:,k+1) = P1p;
    P2(:,:,:,k+1) = P2p;
end

%% calculs moments a posteriori :
rpost = zeros(2,N);
Prpost = zeros(2,2,N);
m1post = zeros(2,N);
m2post = zeros(2,N);
P1post = zeros(2,2,N);
P2post = zeros(2,2,N);
for k=1:N
    for i = 1:Part
        m1post(:,k) = m1post(:,k)+ w(:,i,k)*m1(:,i,k);
        m2post(:,k) = m2post(:,k)+ w(:,i,k)*m2(:,i,k);
        rpost(:,k)  = rpost(:,k)+ w(:,i,k)*r(:,i,k);
    end
    for i = 1:Part
        P1post(:,:,k) = P1post(:,:,k) + ((m1(:,i,k)-m1post(:,k))*((m1(:,i,k)-m1post(:,k))')+P1(:,:,i,k))*w(1,i,k);
        P2post(:,:,k) = P2post(:,:,k) + ((m2(:,i,k)-m2post(:,k))*((m2(:,i,k)-m2post(:,k))')+P2(:,:,i,k))*w(1,i,k);
        Prpost(:,:,k) = Prpost(:,:,k) + (r(:,i,k)-rpost(:,k))*((r(:,i,k)-rpost(:,k))')*w(1,i,k);
    end
end

%% Affichage

%reformatage donnees pour s'interfacer avec les fonctions d'affichage
Eest = zeros(6,N);
Pest = zeros(6,6,N);
for t =1:N
    Eest(1:2,t) = rpost(:,t);
    Eest(3:4,t) = m1post(:,t);
    Eest(5:6,t) = m2post(:,t);
    Pest(1:2,1:2,t) = Prpost(:,:,t);
    Pest(3:4,3:4,t) = P1post(:,:,t);
    Pest(5:6,5:6,t) = P2post(:,:,t);
end

figure(1)
courbes(Xreel,Eest,Pest,w);

figure(2)
plotEllipse(Xreel,[rpost;m1post;m2post],Pest,nan(6,Part,N),[],0);