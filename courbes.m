function courbes(x,xest,pest,w)
    T = size(x);
    T = T(2);
    titles = ["x robot", "y robot", "x amer 1","y amer 1", "x amer 2","y amer 2","Neff lin.","Neff log"];
    for a = 1:6
        subplot(8,1,a);
        plot(1:T,x(a,:),'b',1:T,xest(a,:)+3*(reshape(pest(a,a,:),[1,T]).^(0.5)),'r',1:T,xest(a,:)-3*(reshape(pest(a,a,:),[1,T]).^(0.5)),'r');
        title(titles(a));
        hold on
    end
    weff = nan(1,T);
    for t = 1:T
        weff(1,t) = 1/sum(w(:,:,t).^2);
    end
    subplot(8,1,7);
    plot(1:T,weff);
    title(titles(7));
    hold on
    subplot(8,1,8);
    plot(1:T,weff);
    title(titles(8));
    set(gca, 'YScale', 'log')
end