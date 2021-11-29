function courbes(x,xest,pest,w)
    T = size(x);
    T = T(2);
    for a = 1:6
        subplot(8,1,a);
        plot(1:T,x(a,:),'b',1:T,xest(a,:)+3*(reshape(pest(a,a,:),[1,T]).^(0.5)),'r',1:T,xest(a,:)-3*(reshape(pest(a,a,:),[1,T]).^(0.5)),'r');
        hold on
    end
    weff = nan(1,T);
    for t = 1:T
        weff(1,t) = 1/sum(w(:,:,t).^2);
    end
    subplot(8,1,7);
    plot(1:T,weff);
    hold on
    subplot(8,1,8);
    plot(1:T,weff);
    set(gca, 'YScale', 'log')
end