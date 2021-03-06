function th = homPlast(th, frAvg, prefFR, lambda)

global dt

   dth = dt*(log((frAvg)./(prefFR))) * lambda;
   dth(isnan(dth)) = -dt*lambda;
   th = th + dth;
   
end