x_ic = zeros(12,1); u = 0; tspan = [0 20];
[t_out,x_out] = ode45(@monospinnerDynamics(t,x),tspan,x_ic);