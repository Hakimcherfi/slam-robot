clear all;
close all;
clc;
plot_p = 0; %0 pour ne rien afficher, 1 pour voir la simu
[N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,X] = simulationDonnees(plot_p);
Part = 1000;
xresamp = 0.5; %entre 0 et 1 (0 : jamais, 1 : toujours) 
X = zeros(6,Part,N);
W = zeros(1,Part,N);
%k=0 :
[x,w] = SIR(zeros(6,Part),zeros(1,Part),[],0,1); %[] : mesure vide, 0 : k, 1 : pas resamp. 
X(:,:,1) = x;
W(:,:,1) = w;

%k>0 :
for k = 1:N-1
    [x,w] = SIR(X(:,:,k),W(:,:,k),Z(:,k),k,1+(Part-1)*xresamp); %x : (6,Part) w : (1,Part)
    X(:,:,k+1)=x;
    W(:,:,k+1)=w;
end
