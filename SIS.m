function [x,w] = SIS(x1,w1,z,k)
% Fonction SIS
%     Tire de nouveaux échantillons et évalue leurs poids associés, à partir
%     des échantillons à l'instant précédent, de leurs poids à l'instant précédent
%     et de la mesure à assimiler. La fonction d'importance est la dynamique
%     a priori.
%     Si l'instant k vaut 0 alors les données en entrée sont ignorées et des
%     échantillons sont tirés selon la distribution initiale.
% Entrée :
%     x1 : tableau de (6 lignes x "Part" colonnes), correspondant aux échantillons
%     à l'instant précédent
%     w1 : tableau de (1 ligne x "Part" colonnes), correspondant aux poids des
%     échantillons à l'instant précédent
%     z : tableau de 4 lignes, correspondant à la mesure à assimiler, pouvant
%     des NaN
%     k : instant courant
% Sorties :
%     x : tableau de (6 lignes x "Part" colonnes), correspondant aux nouveaux
%     échantillons, échantillonnés selon la dynamique a priori
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
    %Qw = diag([.1^2 .1^2 (1e-2)^2 (1e-2)^2 (1e-2)^2 (1e-2)^2]);
    %Rv = diag([.1^2 .1^2 .1^2 .1^2])*10;
    
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
            x(:,i)=F*x1(:,i)+chol(Qw)'*randn(6,1);
            if     ~isnan(z(1,:)) && ~isnan(z(3,:)) %2 amers visibles
                w(:,i)=w1(:,i)*exp(-0.5*(z-Hfull*x(:,i))'/Rv*(z-Hfull*x(:,i))); % cas 2 amers visibles et importance = dynamique a priori
                assert(~isnan(w(:,i)));
            elseif ~isnan(z(1,:)) &&  isnan(z(3,:)) %premier amer visible
                w(:,i)=w1(:,i)*exp(-0.5*(z(1:2,:)-H1*x(:,i))'/(Rv(1:2,1:2))*(z(1:2,:)-H1*x(:,i))); % cas premier amer visible et importance = dynamique a priori
                assert(~isnan(w(:,i)));
            elseif  isnan(z(1,:)) && ~isnan(z(3,:))%deuxieme amer visible
                w(:,i)=w1(:,i)*exp(-0.5*(z(3:4,:)-H2*x(:,i))'/(Rv(3:4,3:4))*(z(3:4,:)-H2*x(:,i))); % cas deuxieme amer visible et importance = dynamique a priori
                assert(~isnan(w(:,i)));
            else %aucun amer visible
                w(:,i)=w1(:,i);
                assert(~isnan(w1(:,i)));
            end
        end
        assert(sum(w)~=0);
        w = w/(sum(w));
    end
end