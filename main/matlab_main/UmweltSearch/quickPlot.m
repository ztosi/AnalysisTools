function quickPlot(neuNo, map, mapHid, outputs, randS, ce, acc) 

figure('Position', [0 100 1800 450]);
subplot(131);
hold on;

scatter3(map(outputs(randS)==0,1), map(outputs(randS)==0,2),...
    map(outputs(randS)==0,3), 20, repmat([0 0.05 .3],sum(outputs(randS)==0),1), 'filled');
scatter3(map(outputs(randS)==1,1), map(outputs(randS)==1,2),...
    map(outputs(randS)==1,3), 20, repmat([0 .4 0.02],sum(outputs(randS)==1),1), 'filled');
legend('Did not spike', 'Spiked');
title({['Neu: ' num2str(neuNo), ' tSNE Projection Evolved Input States'], ...
    '(color: if the state came before a spike)'});
hold off;
subplot(132);
hold on;
scatter3(mapHid(outputs(randS)==0,1), mapHid(outputs(randS)==0,2),...
    mapHid(outputs(randS)==0,3), 20, repmat([0 0.05 .3], sum(outputs(randS)==0), 1), 'filled');
scatter3(mapHid(outputs(randS)==1,1), mapHid(outputs(randS)==1,2),...
    mapHid(outputs(randS)==1,3), 20, repmat([0 .4 0.02], sum(outputs(randS)==1), 1), 'filled');
legend('Did not spike', 'Spiked');
title({['Neu: ' num2str(neuNo), ' tSNE Projection Evolved Hidden States'], ...
    '(color: if the state came before a spike)'});
hold off;
subplot(133);
hold on;
histogram(ce, 'Normalization', 'probability', 'FaceColor', [0.01 .35 0.02], 'FaceAlpha', 1);
annotation('textbox', [.75 .1, .3 .3], 'EdgeColor', [1 1 1], ...
    'String', ['Acc: ' num2str(acc*100) '%'], 'FitBoxToText', 'on', ...
    'FontWeight', 'bold', 'FontSize', 12);
title(['Neu: ', num2str(neuNo), ' Performance']);
xlabel('Cross-Entropy');
hold off;