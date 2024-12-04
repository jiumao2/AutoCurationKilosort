function overlap_ratio = overlapRatio(mdl)
normpdf1 = @(x)normpdf(x, mdl.mu(1), mdl.Sigma(1))./normpdf(0,0,mdl.Sigma(1));
normpdf2 = @(x)normpdf(x, mdl.mu(2), mdl.Sigma(2))./normpdf(0,0,mdl.Sigma(2));

x_test = -10000:0.1:10000;

pdf1 = normpdf1(x_test);
pdf2 = normpdf2(x_test);
min_pdf = min(pdf1, pdf2);

overlap_ratio = max(sum(min_pdf)./sum(pdf1), sum(min_pdf)./sum(pdf2));
end

% figure;
% plot(pdf1);
% hold on;
% plot(pdf2);
