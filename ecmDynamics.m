function [x_next, y] = ecmDynamics(x, u, params)

    z = x(1);
    U1 = x(2);
    T = x(3);
    I = u;

    dt = params.dt;
    
    R0 = params.R0(z); 
    R1 = params.R1(z);
    C1 = params.C1(z);
    OCV = params.OCV(z);

    z_next = z + params.eta * dt * I / params.Q;
    U1_next = (U1 - R1*I) * exp(-dt / (R1*C1)) + R1*I;
    T_next = T + ((dt / params.Cth) * (I^2 * (R0 + R1) - (T - params.Tamb)/params.Rth));

    VT = OCV - U1 - R0 * I;

    x_next = [z_next; U1_next; T_next];
    y = VT;
end