figure; 

[bhsz, ed] = histcounts(szes, 'BinMethod', 'integer', 'Normalization', 'pdf');
binsz = (0.5*diff(ed))+ed(1:end-1);

%unqAv = cellfun(@unique, cellfun(@find, str_avas, 'UniformOutput', false), 'UniformOutput', false);
%unqSzs = cellfun(@length, unqAv);
[bhUsz, ed] = histcounts(szes, logspace(floor(log10(min(szes))), ceil(log10(max(szes))), 100), 'Normalization', 'pdf');
binUsz = (0.5*diff(ed))+ed(1:end-1);

[bhLen, ed] = histcounts(lens, 'BinMethod', 'integer', 'Normalization', 'pdf');
binLen = (0.5*diff(ed))+ed(1:end-1);

subplot(311);
scatter(binsz, bhsz);
hold on; 
plot(10:100, (10:100).^(-(2)), 'k--');
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
subplot(312);
scatter(binUsz, bhUsz);
hold on; 
plot(10:100, (10:100).^(-(2)), 'k--');
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
subplot(313);
scatter(binLen, bhLen);
hold on; 
plot(10:100, (10:100).^(-(2)), 'k--');
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');