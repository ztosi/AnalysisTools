function [transform, lambda] = autoBoxCox(dat)

global data;

data = dat;

if any(dat <=0 )
    lambdas = fminsearch(@bcReal, [0, 1]);
    if lambdas(1) == 0
        dat = log(dat + lambdas(2));
    else
        dat = ((dat+lambdas(2)).^lambdas(1) - 1)/lambdas(1);
    end
    
else
    lambdas = fminsearch(@bcPos,0);
    if lambdas(1) == 0
        dat = log(dat);
    else
        dat = ((dat).^lambdas(1) - 1)/lambdas(1);
    end
end
%transform = (dat-mean(dat))./std(dat);
lambda = lambdas;


end