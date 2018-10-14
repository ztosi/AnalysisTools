figure;
order = [1 4 7 19 22 25 37 40 43];
for ii=1:8
    [~, ed] = histcounts(log10(PrefFRs));
    ipt = order(ii);
    s=subplot(6,9,[ipt, ipt+1, ipt+9, ipt+10]);
    h=histogram(log10(PrefFRs(group{ii})), ed, 'Normalization', 'pdf');
    h.FaceColor = [0 0 0];
    h.FaceAlpha = 1;
    annotation('textbox', [.2 .5 .3 .3], 'Position', [.1 .1 .3 .1], ...
        'FitBoxToText', 'on', 'String', ...
        {['\mu = ' num2str(mean(log10(PrefFRs(group{ii}))))], ...
        [ '\sigma = ' num2str(std(log10(PrefFRs(group{ii}))))]});
    xlabel('Log Firing Rate');
    title(['Module ' num2str(ii)]);
    s.Box = 'off';
    subplot(6,9, order(ii)+2);
    
    labels = {'', ''};
    p1=pie([sum(sum(inMat(:,group{ii}))), ...
        sum(sum(wtMat(ei, group{ii})))], labels);
    %if ii==1
        title('Inp/Rec Drive');
    %end
    subplot(6,9, order(ii)+11);
    p2=pie([length(intersect(find(ei), group{ii})), ...
        length(intersect(find(~ei), group{ii}))], labels);
    %if ii==1
        title('Exc./Inh. Pop.');
    %end
    p2(1).FaceColor = [0.05 0.7 .2];
    if length(p2)>3
        p2(3).FaceColor = [0.7 0.05 0.1];
    end
end
