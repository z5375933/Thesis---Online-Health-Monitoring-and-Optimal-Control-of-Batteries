function rbf = buildCenters()
    
    VTRange = linspace(2.5, 4.2, 4);  
    TRange  = linspace(280, 318, 4);  
    
    centers = [];
    for vt = VTRange
        for t = TRange
            centers = [centers [vt; t]];
        end
    end
    

    widths = ones(16,1) * 0.05;
    
    rbf.centers = centers;
    rbf.widths  = widths;
end