function [x,w] = SIR(x1,w1,z,k,Ns)
    [x,w] = SIS(x1,w1,z,k);
    Neff = 1/(sum(w.^2));
    if Neff < Ns
        [x,w] = resampling(x,w);
    end
end