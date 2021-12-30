function [x,w] = SIS2(x1,w1,z,k)
% Fonction SIS
%     Tire de nouveaux échantillons et évalue leurs poids associés, à partir
%     des échantillons à l'instant précédent, de leurs poids à l'instant précédent
%     et de la mesure à assimiler. La fonction d'importance tient compte de la mesure
%     pour la partie amer 1 et amer 2, et uniquement la dynamique a priori pour la partie
%     robot
%     Si l'instant k vaut 0 alors les données en entrée sont ignorées et des
%     échantillons sont tirés selon la distribution initiale.
%     
% Entrée :
%     x1 : tableau de (6 lignes x "Part" colonnes), correspondant aux échantillons
%     à l'instant précédent
%     w1 : tableau de (1 ligne x "Part" colonnes), correspondant aux poids des
%     échantillons à l'instant précédent
%     z : tableau de 4 lignes, correspondant à la mesure à assimiler,
%     pouvant contenir des NaN
%     k : instant courant
% Sorties :
%     x : tableau de (6 lignes x "Part" colonnes), correspondant aux nouveaux
%     échantillons, tirés selon la loi d'importance
%     w : tableau de (1 ligne x "Part" colonnes), correspondant aux poids des
%     échantillons générés à l'instant courant par assimilation de la mesure
    Part = size(x1);
    Part = Part(2);
    deltaT = 1;
    w=pi/4;
    mX0 = [2;0;-1;4;4;-1];
    PX0 = diag([.1 .2 .2 .2 .2 .2]);
    Qw = diag([.05^2 .05^2 (1e-10)^2 (1e-10)^2 (1e-10)^2 (1e-10)^2]);
    Rv = diag([.1^2 .1^2 .1^2 .1^2]);
    
    F = blkdiag([cos(w*deltaT) -sin(w*deltaT) ; sin(w*deltaT) cos(w*deltaT)], eye(2), eye(2));
    H1=[-1 0 1 0 0 0 ; 0 -1 0 1 0 0];
    H2=[-1 0 0 0 1 0 ; 0 -1 0 0 0 1];
    H=[H1;H2];
    Hfull=H;
    x = nan(6,Part);
    w = nan(1,Part);
    if (k==0)
        for i=1:Part
            x(:,i) = mX0+chol(PX0)'*randn(6,1);
            w(:,i)=1/Part;
        end
    else
        for i = 1:Part
            if     ~isnan(z(1,:)) && ~isnan(z(3,:)) %2 amers visibles
                %echantillonnage
                x(1:2,i) = F(1:2,1:2)*x1(1:2,i)+chol(Qw(1:2,1:2))'*randn(2,1);%partie robot
                x(3:4,i) = F(1:2,1:2)*x1(1:2,i)+z(1:2,:)+chol(Rv(1:2,1:2))'*randn(2,1);%partie amer 1
                x(5:6,i) = F(1:2,1:2)*x1(1:2,i)+z(3:4,:)+chol(Rv(3:4,3:4))'*randn(2,1);%partie amer 2
                
                %evaluation
                pzkxk = 1/sqrt(det(2*pi*Rv))*exp(-0.5*(z-Hfull*x(:,i))'/Rv*(z-Hfull*x(:,i)));
                pm1mk1 = 1;%/sqrt(det(2*pi*Qw(3:4,3:4)))*exp(-0.5*(x(3:4,i)-x1(3:4,i))'/(Qw(3:4,3:4))*(x(3:4,i)-x1(3:4,i)));
                pm2mk1 = 1;%/sqrt(det(2*pi*Qw(5:6,5:6)))*exp(-0.5*(x(5:6,i)-x1(5:6,i))'/(Qw(5:6,5:6))*(x(5:6,i)-x1(5:6,i)));
                q1 = 1/sqrt(det(2*pi*Rv(1:2,1:2)))*exp(-0.5*(z(1:2,:)-(x(3:4,i)-x(1:2,i)))'/(Rv(1:2,1:2))*(z(1:2,:)-(x(3:4,i)-x(1:2,i))));
                q2 = 1/sqrt(det(2*pi*Rv(3:4,3:4)))*exp(-0.5*(z(3:4,:)-(x(5:6,i)-x(1:2,i)))'/(Rv(3:4,3:4))*(z(3:4,:)-(x(5:6,i)-x(1:2,i))));
                w(:,i) = w1(:,i)*pzkxk*pm1mk1*pm2mk1/(q1*q2);
                assert(~isnan(w(:,i)));
            
            elseif ~isnan(z(1,:)) &&  isnan(z(3,:)) %premier amer visible
                %echantillonnage
                x(1:2,i) = F(1:2,1:2)*x1(1:2,i)+chol(Qw(1:2,1:2))'*randn(2,1);%partie robot
                x(3:4,i) = F(1:2,1:2)*x1(1:2,i)+z(1:2,:)+chol(Rv(1:2,1:2))'*randn(2,1);%partie amer 1
                x(5:6,i) = x1(5:6,i)+chol(Qw(5:6,5:6))'*randn(2,1);%partie amer 2
                
                %evaluation
                pzkxk = 1/sqrt(det(2*pi*Rv(1:2,1:2)))*exp(-0.5*(z(1:2,:)-H1*x(:,i))'/(Rv(1:2,1:2))*(z(1:2,:)-H1*x(:,i)));
                pm1mk1 = 1;%/sqrt(det(2*pi*Qw(3:4,3:4)))*exp(-0.5*(x(3:4,i)-x1(3:4,i))'/(Qw(3:4,3:4))*(x(3:4,i)-x1(3:4,i)));
                q1 = 1/sqrt(det(2*pi*Rv(1:2,1:2)))*exp(-0.5*(z(1:2,:)-(x(3:4,i)-x(1:2,i)))'/(Rv(1:2,1:2))*(z(1:2,:)-(x(3:4,i)-x(1:2,i))));
                w(:,i) = w1(:,i)*pzkxk*pm1mk1/q1;
                assert(~isnan(w(:,i)));
            
            elseif  isnan(z(1,:)) && ~isnan(z(3,:))%deuxieme amer visible
                %echantillonnage
                x(1:2,i) = F(1:2,1:2)*x1(1:2,i)+chol(Qw(1:2,1:2))'*randn(2,1);%partie robot
                x(3:4,i) = x1(3:4,i)+chol(Qw(3:4,3:4))'*randn(2,1);%partie amer 1
                x(5:6,i) = F(1:2,1:2)*x1(1:2,i)+z(3:4,:)+chol(Rv(3:4,3:4))'*randn(2,1);%partie amer 2
                
                %evaluation
                pzkxk = 1/sqrt(det(2*pi*Rv(3:4,3:4)))*exp(-0.5*(z(3:4,:)-H2*x(:,i))'/(Rv(3:4,3:4))*(z(3:4,:)-H2*x(:,i)));
                pm2mk1 = 1;%/sqrt(det(2*pi*Qw(5:6,5:6)))*exp(-0.5*(x(5:6,i)-x1(5:6,i))'/(Qw(5:6,5:6))*(x(5:6,i)-x1(5:6,i)));
                q2 = 1/sqrt(det(2*pi*Rv(3:4,3:4)))*exp(-0.5*(z(3:4,:)-(x(5:6,i)-x(1:2,i)))'/(Rv(3:4,3:4))*(z(3:4,:)-(x(5:6,i)-x(1:2,i))));
                w(:,i) = w1(:,i)*pzkxk*pm2mk1/q2;
                assert(~isnan(w(:,i)));
            
            else %aucun amer visible
                x(:,i)=F*x1(:,i)+chol(Qw)'*randn(6,1);
                w(:,i)=w1(:,i);
                assert(~isnan(w(:,i)));
            end
        end
        assert(sum(w)~=0);
        w = w/(sum(w));
    end
end