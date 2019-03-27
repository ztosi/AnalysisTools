function prefFRs = metaHP(adjMat, inD, prefFRs, estFRs, eta, b)
    global N alpha lowFR beta dt;
    % Adj is supposed to be transposed

%     max = 100;
%     min = 0.01;
%     
%     mat = adjMat .* (log(estFRs+0.00001)-log(estFRs'+0.0001));
%     
%     mat = 2/b .* mat .* exp(-((mat.^2)./(b))) .* -tanh(adjMat.*(abs(mat)-1));
%     
% %     pfg1=prefFRs.*((prefFRs>1).* -1/(max-1) + (1+1/max)); 
% %     pfl1 = prefFRs.*((prefFRs<=1) .* (1+min) - min);
%     
%     dpfr = (sum(mat)'./1000); %.* sqrt(1-(log(prefFRs)/4.5).^2);
%     
%     prefFRs = exp(log(prefFRs) + dt * eta * dpfr);
    
    
  
    dfs = estFRs-estFRs';
    
    dfs = dfs .* adjMat;
    
    ltVals = (exp(-(((dfs>0).*dfs)./prefFRs').^2));
    ltVals = 1.5* sum(ltVals .*adjMat)';
    gtVals = (exp(-( ((dfs<0).*dfs)./prefFRs').^2));
    gtVals = sum(gtVals .* adjMat)';

    f_plus = exp(-prefFRs./(beta*lowFR));
    lf = prefFRs <= lowFR;
    f_minus_L = prefFRs(lf)./lowFR;
    f_minus = (1 + log(1+alpha*(prefFRs(~lf)./lowFR -1))./alpha);
    gtVals(lf) = gtVals(lf) .* f_minus_L;
    gtVals(~lf) = gtVals(~lf) .* f_minus;
    ltVals = ltVals .* f_plus;
    
    prefFRs = prefFRs + eta * dt .* (ltVals - gtVals) ./ inD .* ((0.1.*randn(N,1))+1);
    

end