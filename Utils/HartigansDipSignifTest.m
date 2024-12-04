function [h, dip] = HartigansDipSignifTest(xpdf, alpha, dip_crit_all, nboot)
% HartigansDipSignifTest
%
% calculates Hartigan's DIP statistic and its significance for the empirical
% p.d.f  XPDF (vector of sample values).
%
% This routine calls the matlab routine 'HartigansDipTest' that actually
% calculates the DIP NBOOT is the user-supplied sample size of boot-strap
% Code by F. Mechler (27 August 2002)
% Modified by Yue Huang (4 December 2024)
%

% calculate the DIP statistic from the empirical pdf
if nargin < 3
    dip_crit_all = load('dip_crit.mat');
end

if nargin < 4
    nboot = 1000;
end

% sort and normalize to be in 0..1
dip = HartigansDipTest(xpdf);
N = length(xpdf);

if alpha == 0.1
    dip_crit = dip_crit_all.dip_crit_p90(N);
elseif alpha == 0.05
    dip_crit = dip_crit_all.dip_crit_p95(N);
elseif alpha == 0.01
    dip_crit = dip_crit_all.dip_crit_p99(N);
elseif alpha == 0.001
    dip_crit = dip_crit_all.dip_crit_p999(N);
else
    % calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
    boot_dip = zeros(nboot);
    for i = 1:nboot
       unifpdfboot = sort(unifrnd(0,1,1,N));
       boot_dip(i) = HartigansDipTest(unifpdfboot);
    end
    
    dip_crit = prctile(boot_dip, 100*(1-alpha));
%     p_value = sum(boot_dip > dip) / nboot;
end

if dip > dip_crit
    h = 0;
else
    h = 1;
end

