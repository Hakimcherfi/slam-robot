function courbes(x,xest,pest)
    T = size(x);
    T = T(2);
    for a = 1:6
        subplot(6,1,a);
        plot(1:T,x(a,:),'b',1:T,xest(a,:)+3*(reshape(pest(a,a,:),[1,T]).^(0.5)),'r',1:T,xest(a,:)-3*(reshape(pest(a,a,:),[1,T]).^(0.5)),'r');
        hold on
    end
end