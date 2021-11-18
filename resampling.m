function [x,w] = resampling(xi,wi)
    Part = size(xi);
    Part = Part(2);
    x = zeros(size(xi));
    ci = cumsum(wi);
    u = rand()/Part;
    i = 1;
    for j = 1:Part
        u = u+(j-1)/Part;
        while u>ci(i)
            i = i+1;
        end
        x(:,j)=xi(:,i);
    end
    w = ones(size(wi))/Part;
end