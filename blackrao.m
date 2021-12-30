function [w,r,m1,m2,P1,P2] = blackrao(r1,m11,m21,P11,P21,w1,z,k)
    %entrees :
    %r1 : particules partie robot
    %m11 : particules moment 1 amer 1 instant precedent
    %P11 : particules moment 2 amer 1 instant precedent
    %w1 : poids particules instant precedent
    %z : mesure
    %k : instant

    %sorties :
    %w : poids
    %r : particules etat robot
    %m1 : particules moment 1 amer 1
    %P1 : particules moment 2 amer 1
    %m2
    %P2
    
    Part = size(r1);
    Part = Part(2);
    deltaT = 1;
    w=pi/4;
    mX0 = [2;0;-1;4;4;-1];
    PX0 = diag([.1 .2 .2 .2 .2 .2]);
    Qw = diag([.05^2 .05^2 (1e-10)^2 (1e-10)^2 (1e-10)^2 (1e-10)^2]);
    Rv = diag([.1^2 .1^2 .1^2 .1^2]);    
    F = blkdiag([cos(w*deltaT) -sin(w*deltaT) ; sin(w*deltaT) cos(w*deltaT)], eye(2), eye(2));
    %H1=[-1 0 1 0 0 0 ; 0 -1 0 1 0 0];
    %H2=[-1 0 0 0 1 0 ; 0 -1 0 0 0 1];
    %H=[H1;H2];
    %Hfull=H;
    w = nan(1,Part);
    r = nan(2,Part);
    m1 = nan(2,Part);
    m2 = nan(2,Part);
    P1 = nan(2,2,Part);
    P2 = nan(2,2,Part);
    if (k==0)
        for i=1:Part
            r(:,i) = mX0(1:2)+chol((PX0(1:2,1:2))')*randn(2,1);
            m1(:,i) = mX0(3:4);
            m2(:,i) = mX0(5:6);
            P1(:,:,i) = PX0(3:4,3:4);
            P2(:,:,i) = PX0(5:6,5:6);
            w(:,i)=1/Part;
        end
   else
         for i = 1:Part
             r(:,i)=F(1:2,1:2)*r1(:,i)+chol(Qw(1:2,1:2))'*randn(2,1);
             %Gain Kalman: K
             if     ~isnan(z(1,:)) && ~isnan(z(3,:)) %2 amers visibles
                  d1 = z(1:2)-(m11(:,i)-r(:,i));
                  d2 = z(3:4)-(m21(:,i)-r(:,i));
                  w(:,i)=w1(:,i)*exp(-0.5*d1'/(P11(:,:,i)+Rv(1:2))*d1)*exp(-0.5*d2'/(P21(:,:,i)+Rv(3:4))*d2); % cas 2 amers visibles et importance = dynamique a priori
                  K1 = P11(:,:,i)/(Rv(1:2,1:2)+P11(:,:,i));
                  m1(:,i) = m11(:,i)+K1*(z(1:2)+r(:,i)-m11(:,i));
                  P1(:,:,i) = (eye(2)-K1)*P11(:,:,i);
                  K2 = P21(:,:,i)/(Rv(3:4,3:4)+P21(:,:,i));
                  m2(:,i) = m21(:,i)+K2*(z(3:4)+r(:,i)-m21(:,i));
                  P2(:,:,i) = (eye(2)-K2)*P21(:,:,i);
                  assert(~isnan(w(:,i)));
                  
             elseif ~isnan(z(1,:)) &&  isnan(z(3,:)) %premier amer visible
                  d1 = z(1:2)-(m11(:,i)-r(:,i));
                  w(:,i)=w1(:,i)*exp(-0.5*d1'/(P11(:,:,i)+Rv(1:2,1:2))*d1);
                  K1 = P11(:,:,i)/(Rv(1:2,1:2)+P11(:,:,i));
                  m1(:,i) = m11(:,i)+K1*d1;
                  P1(:,:,i) = (eye(2)-K1)*P11(:,:,i);
                  m2(:,i) = m21(:,i);
                  P2(:,:,i) = P21(:,:,i);
                  assert(~isnan(w(:,i)));
                  
             elseif  isnan(z(1,:)) && ~isnan(z(3,:))%deuxieme amer visible
                  d2 = z(3:4)-(m21(:,i)-r(:,i));
                  w(:,i)=w1(:,i)*exp(-0.5*d2'/(P21(:,:,i)+Rv(3:4))*d2);
                  K2 = P21(:,:,i)/(Rv(3:4,3:4)+P21(:,:,i));
                  m2(:,i) = m21(:,i)+K2*d2;
                  P2(:,:,i) = (eye(2)-K2)*P21(:,:,i);
                  m1(:,i) = m11(:,i);
                  P1(:,:,i) = P11(:,:,i);
                  assert(~isnan(w(:,i)));
                  
             else %aucun amer visible
                  w(:,i)=w1(:,i);
                  m1(:,i) = m11(:,i);
                  P1(:,:,i) = P11(:,:,i);
                  m2(:,i) = m21(:,i);
                  P2(:,:,i) = P21(:,:,i);
                  assert(~isnan(w1(:,i)));
             end
         end
         assert(sum(w)~=0);
         w = w/(sum(w));
    end
end