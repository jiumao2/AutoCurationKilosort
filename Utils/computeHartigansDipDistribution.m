% HartigansDipSignifTestLookUpTable
N = sort(unique([2:10, round(1.2.^(1:64)+1), 1e5, 5e5, 1e6]), 'descend');
nboot = 10000;

percentile_90 = zeros(1,length(N));
percentile_95 = zeros(1,length(N));
percentile_99 = zeros(1,length(N));
percentile_999 = zeros(1,length(N));

for k = 1:length(N)
    dips = zeros(1,nboot);
    for j = 1:nboot
        dips(j) = HartigansDipTest(unifrnd(0,1,1,N(k)));
    end
    
    percentile_90(k) = prctile(dips, 90);
    percentile_95(k) = prctile(dips, 95);
    percentile_99(k) = prctile(dips, 99);
    percentile_999(k) = prctile(dips, 99.9);

    fprintf('%d / %d done!\n', k, length(N));
end

%%
figure;
ax = axes;
idx = 1:length(N)-2;
hold on;
plot(ax, N(idx), percentile_90(idx), '.');
plot(ax, N(idx), percentile_95(idx), '.');
plot(ax, N(idx), percentile_99(idx), '.');
plot(ax, N(idx), percentile_999(idx), '.');

set(gca, 'xscale', 'log', 'yscale', 'log');

%% log-log is a straight line
mdl_90 = polyfit(log(N(idx)), log(percentile_90(idx)), 1);
mdl_95 = polyfit(log(N(idx)), log(percentile_95(idx)), 1);
mdl_99 = polyfit(log(N(idx)), log(percentile_99(idx)), 1);
mdl_999 = polyfit(log(N(idx)), log(percentile_999(idx)), 1);

dip_crit_p90 = @(x) exp(interp1(log(N(idx)), log(percentile_90(idx)), log(x), 'linear', 'extrap'));
dip_crit_p95 = @(x) exp(interp1(log(N(idx)), log(percentile_95(idx)), log(x), 'linear', 'extrap'));
dip_crit_p99 = @(x) exp(interp1(log(N(idx)), log(percentile_99(idx)), log(x), 'linear', 'extrap'));
dip_crit_p999 = @(x) exp(interp1(log(N(idx)), log(percentile_999(idx)), log(x), 'linear', 'extrap'));

save dip_crit.mat dip_crit_p90 dip_crit_p95 dip_crit_p99 dip_crit_p999;

format long
disp(mdl_90);
disp(mdl_95);
disp(mdl_99);
disp(mdl_999);

plot(ax, N(idx), exp(polyval(mdl_90, log(N(idx)))), '-');
plot(ax, N(idx), exp(polyval(mdl_95, log(N(idx)))), '-');
plot(ax, N(idx), exp(polyval(mdl_99, log(N(idx)))), '-');
plot(ax, N(idx), exp(polyval(mdl_999, log(N(idx)))), '-');

