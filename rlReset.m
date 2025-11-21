function x0 = rlReset(x0List)
    idx = randi(numel(x0List));
    x0  = x0List{idx};
end
