function [nlogL]  = bcReal(lamb1, lamb2)

global data;

    dat = data;
    if lamb1 == 0
        dat = log(dat+lamb2);
    else
        dat = ((dat+lamb2).^lamb1 - 1) ./ lamb1;      
    end
    mn = mean(dat);
    st = std(dat);
    nlogL = normlike([mn st], dat);

end
