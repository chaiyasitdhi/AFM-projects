function F = hertz_evaluate(depth,p)
    F = (1/sqrt(2))*p(1)/(1-(p(2).^2))*(tan(p(3))).*(depth.^2);
end