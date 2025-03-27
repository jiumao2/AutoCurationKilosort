function xtickangle(axes_all, angle)

if ~iscell(axes_all)
    xtickangle(axes_all, angle);
    return
end

for k = 1:size(axes_all,1)
    for j = 1:size(axes_all,2)
        xtickangle(axes_all{k,j}, angle);
    end
end

end