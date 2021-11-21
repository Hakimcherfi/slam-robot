function plotEllipse(x,xest,pest,xpart,wpart,b)
    T = size(x);
    T = T(2);
    Npart = size(xpart);
    Npart = Npart(2);
    for k = 1:T
        axis([-3 6 -3 6])
        for l = 1:3
           MatricePestk = pest(:,:,k);
           figure(2);
           ellipse(xest(2*l-1:2*l,k),MatricePestk(2*l-1:2*l,2*l-1:2*l),'r');
           hold on;
           scatter(x(2*l-1,k),x(2*l,k),25,'filled');
           hold on;
           if b == 1
               wmax = max(wpart(1,:,k));
               for i = 1:Npart
                    scatter(xpart(2*l-1,i,k),xpart(2*l,i,k),wpart(1,i,k)/wmax*50,'d')
               end
           end
        end
         
        input(['Actuellement, k = ' num2str(k-1) '. Taper entree pour avancer']);
        clf;
    end
end