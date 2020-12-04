function prefFRs = metaHP(adjMat, estFRs, eta)
    global N dt mu sig rt2;
    % Adj is supposed to be transposed

   estFRs(estFRs<1e-5) = 1e-5;
   lgE = log(estFRs);

   mat = -adjMat .* (lgE-lgE');
    
   p = adjMat .* 0.5 .*(1+erf(mat./(sig.*rt2)));
   p = adjMat.*rand(size(mat)).*2.*((rand(size(mat))<p)-.5);
   
   p = p .* 250 .*exp(-4*mat.^2)./N;
   
   a = 0.5 .*(1+erf((mu-lgE)./(sig.*rt2)));
   a = rand(size(a)).*2.*((rand(size(a))<a)-.5);
   
   dpfr = sum(p)'+a;
   
   prefFRs = exp(log(estFRs) + dt*eta*dpfr);

end