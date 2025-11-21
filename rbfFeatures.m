function phi = rbfFeatures(z, rbf)
    % z = [VT; T]
    c  = rbf.centers;
    l  = rbf.widths;
    if isscalar(l), l = l*ones(size(c,2),1); end

    diffs = c - z(:);
    d2 = sum(diffs.^2, 1);
    phi = exp(-l(:)'.*d2);
    phi = phi(:);
end