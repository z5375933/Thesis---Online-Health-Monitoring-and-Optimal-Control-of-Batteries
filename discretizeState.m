function s_idx = discretizeState(s, SOCBins, TBins)

    z = s(1);
    T = s(4);

    zBin = find(z >= SOCBins(1:end-1) & z < SOCBins(2:end), 1, 'first');
    if isempty(zBin)
        zBin = numel(SOCBins) - 1;
    end

    TBin = find(T >= TBins(1:end-1) & T < TBins(2:end), 1, 'first');
    if isempty(TBin)
        TBin = numel(TBins) - 1;
    end

    nZ = numel(SOCBins) - 1;

    s_idx = (TBin - 1)*nZ + zBin;
end
