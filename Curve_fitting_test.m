temps = logspace(0.001, 4, 1000);
rates = original(temps, 8.514e-10, 0, 0, 9.5666e-04, -4.4040e-05, 2.3496e-06);
function original_fit = original(T, a_0, a_12, a_1, b_12, b_1, b_32)
    upper  = a_0 + (a_12)*(T.^(1/2)) + a_1*T;
    lower = (T.^(1/6)) + (b_12)*(T.^(1/2)) + (b_1)*(T) + (b_32)*(T.^(3/2));
    original_fit = upper./lower;
end
function new_fit = new(T, a, b, g)
    new_fit = a * ((T./300)^b) * exp(-g/T);
end