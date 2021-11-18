clear all;
close all;
clc;
plot_p = 0; %0 pour ne rien afficher, 1 pour voir la simu
[N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,X] = simulationDonnees(plot_p);
Part = 2;
X = zeros(6,Part,N);
W = zeros(1,Part,N);
%k=0 :
[x,w] = SIR(zeros(6,Part),zeros(1,Part),[],0,Part/3);
X(:,:,1) = x;
W(:,:,1) = w;

%k>0 :
for k = 1:N-1
    [x,w] = SIR(X(:,:,k),W(:,:,k),Z(:,k),k,Part/3); %x : (6,Part) w : (1,Part)
    X(:,:,k+1)=x;
    W(:,:,k+1)=w;
end
