figure;
subplot(231); 
hold on;
scatter(log(PrefFRs(~ei)), pc(~ei)); 
scatter(log(PrefFRs(ei)), pc(ei)); 
[rho, pval] = corr([log(PrefFRs)', pc])
title('Coltrane 2 Ln(TFR) vs. Pop. Coupling');
xlabel(['\rho = ' num2str(rho(1, 2)) '; p < ' num2str(pval(1,2))]);
hold off;

subplot(232);
hold on;
scatter(kIn(~ei), pc(~ei));
scatter(kIn(ei), pc(ei));
[rho, pval] = corr([kIn', pc])
title('Coltrane 2 In-Deg vs. Pop. Coupling'); 
xlabel(['\rho = ' num2str(rho(1, 2)) '; p < ' num2str(pval(1,2))]);
hold off;


subplot(233);
hold on;
[rho, pval] = corr([kOut', pc])
scatter(kOut(~ei), pc(~ei)); 
scatter(kOut(ei), pc(ei)); 
title('Coltrane 2 Out-Deg vs. Pop. Coupling');
xlabel(['\rho = ' num2str(rho(1, 2)) '; p < ' num2str(pval(1,2))]);
hold off;

subplot(234);
hold on;
[rho, pval] = corr([log(PrefFRs)', kIn'])
scatter(log(PrefFRs(~ei)), kIn(~ei)); 
scatter(log(PrefFRs(ei)), kIn(ei)); 
title('Coltrane 2 Ln(TFR) vs. In-Deg'); 
xlabel(['\rho = ' num2str(rho(1, 2)) '; p < ' num2str(pval(1,2))]);
hold off;

subplot(235); 
hold on;
[rho, pval] = corr([kIn', kOut'])
scatter(kIn(~ei), kOut(~ei)); 
scatter(kIn(ei), kOut(ei)); 
title('Coltrane 2 In-Deg vs. Out-Deg');
xlabel(['\rho = ' num2str(rho(1, 2)) '; p < ' num2str(pval(1,2))]);
hold off;

subplot(236);
hold on;
[rho, pval] = corr([sum(wtMat)', pc])
scatter(sum(wtMat(:, ~ei)), pc(~ei));
scatter(sum(wtMat(:, ei)), pc(ei));
title('Coltrane 2 tot In (inc. inhib) vs. Pop. Coupling');
xlabel(['\rho = ' num2str(rho(1, 2)) '; p < ' num2str(pval(1,2))]);
hold off;






