function F = new_hertz_evaluate(contact_point,depth,p)
    F = (1/sqrt(2))*p(1)/(1-(p(2).^2))*(tan(p(3))).*(depth-contact_point).^2;
end