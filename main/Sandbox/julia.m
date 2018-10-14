function julia (N, c, iter, cmap)

pxw = 4/N;
[X, Y] = meshgrid(-2:pxw:2, -1:pxw:1);
mxL = log(iter);
Z = (X+(Y.*1i)).^2 + c;
clear X Y;
[m,n]=size(Z);
[I, J]= find((real(Z).^2 + imag(Z).^2) >= 4);
linind = sub2ind([m,n], I, J);
Zp = zeros(size(Z));
Zp(linind) = 1 - log2(log(abs(Z(linind)))/mxL);
[I, J]= find((real(Z).^2 + imag(Z).^2) < 4);
linind = sub2ind([m,n], I, J);
Z = sparse(I,J, Z(linind));
figure; imagesc(Zp);
for ii=1:iter
    ii
    if mod(ii, 10) == 0
        figure; imagesc(Zp);
    end
Z = Z.^2 + c;
linind = (real(nonzeros(Z)).^2 + imag(nonzeros(Z)).^2) < 4;
[ip, jp, vp] = find(Z);
lp = sub2ind([m, n], ip, jp);
lpp = lp(linind);
Zp(lpp) = iter - log2(log(abs(Z(lpp)))/mxL);
Z = sparse(ip(linind),jp(linind), vp(linind));
end

colormap(cmap);
figure; imagesc(Zp);


end