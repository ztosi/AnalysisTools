function [dly_dist] = delay_distribution(src, tar)
ii=1;
jj=1;
kk=1;
dly_dist=zeros(length(tar),1);
while ii < length(src)-1
    while jj < length(tar)-1 && tar(jj)<=src(ii)
        jj = jj+1;
    end
    while ii < length(src)-1 && src(ii+1)<tar(jj)
        ii=ii+1;
    end
    if tar(jj) < src(ii)
        break;
    end
    dly_dist(kk) = tar(jj)-src(ii);
    ii=ii+1;
    kk=kk+1;
end
dly_dist = nonzeros(dly_dist);

end