function [nlogL]  = bcPos(lamb)

global data;

    dat = data;
    if lamb == 0
        dat = log(dat);
    else
        dat = (dat.^lamb - 1) ./ lamb;      
    end
%     dat = (dat-mean(dat)) ./ std(dat);
    mn = mean(dat);
    st = std(dat);
    nlogL = normlike([mn st], dat);


end

