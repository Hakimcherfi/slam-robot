function [x,w] = resampling(xi,wi)
    Part = size(xi);
    Part = Part(2);
    x = nan(size(xi));
    ci = cumsum(wi);
    u = rand()/Part;
    i = 1;
    for j = 1:Part
        while u>ci(:,i)
            i = i+1;
        end
        x(:,j)=xi(:,i);
        u = u+(1/Part);
    end
    w = ones(size(wi))/Part;
end