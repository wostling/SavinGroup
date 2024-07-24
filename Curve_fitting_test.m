%https://www.mathworks.com/matlabcentral/fileexchange/25872-free-knot-spline-approximation?s_tid=srchtitle
temps = logspace(0.001, 4, 1000);
rates = original(temps, 8.514e-10, 0, 0, 9.5666e-04, -4.4040e-05, 2.3496e-06);
newrates = new(temps, 1, 1, 1);


%x(1) = a(1);
%x(2) = b(1);
%x(3) = g(1); 
F = @(x,xdata)x(1) * ((xdata/300).^x(2)) .* exp(-x(3)./xdata)*100000000000;
x0 = [3.18673748e-10 -1.78641920e-01 4.97568696e-02];

opt=optimoptions('lsqcurvefit');
opt.StepTolerance= 1E-14
opt.FunctionTolerance= 1e-14

[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,temps,rates)

function original_fit = original(T, a_0, ...
    a_12, a_1, b_12, b_1, b_32)
    upper  = a_0 + (a_12)*(T.^(1/2)) + a_1*T;
    lower = (T.^(1/6)) + (b_12)*(T.^(1/2)) + (b_1)*(T) + (b_32)*(T.^(3/2));
    original_fit = 100000000000*upper./lower;
end
function new_fit = new(T, a, b, g)
    new_fit = a * ((T/300).^b) .* exp(-g./T);
end