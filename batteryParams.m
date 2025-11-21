function params = batteryParams()

    SOC = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

    R0_data = [0.0629 0.0614 0.0641 0.0630 0.0613 0.0615 0.0610 0.0612 0.0602];
    R1_data = [0.0392 0.0377 0.0330 0.0312 0.0300 0.0596 0.0436 0.0354 0.0198];
    C1_data = [1075.17 1287.61 1290.36 1268.09 1249.34 905.44 880.62 921.69 1112.36];

    params.R0  = @(z) interp1(SOC, R0_data, z, 'spline', 'extrap');
    params.R1  = @(z) interp1(SOC, R1_data, z, 'spline', 'extrap');
    params.C1  = @(z) interp1(SOC, C1_data, z, 'spline', 'extrap');

    params.Q = 7200;          
    params.eta = 1.0;      
    params.dt = 1;           
    params.Cth = 1.0;        
    params.Rth = 1.3;     
    params.Tamb = 298;   

    params.OCV = @(z) 3.0 + 1.2 * z;
end