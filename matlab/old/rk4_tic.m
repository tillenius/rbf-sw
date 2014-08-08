function [t,H] = rk4_tic(rhs_tic,tspan,H,tasklib)
% rhs - rhs(t,y)
% tspan - time vector [t0, t0+dt, ...,t1]
% H - solution

for i = 2:length(tspan)
    dt = tspan(i)-tspan(i-1);
    t = tspan(i);
    K = H;

    d1 = dt*rhs_tic(t,K,tasklib);

tasklib_tic(tasklib);
    K = H + 0.5*d1;
tasklib_toc(tasklib, 'matmul1');

    d2 = dt*rhs_tic(t+0.5*dt,K,tasklib);

tasklib_tic(tasklib);
    K = H + 0.5*d2;
tasklib_toc(tasklib, 'matmul2');

    d3 = dt*rhs_tic(t+0.5*dt,K,tasklib);

tasklib_tic(tasklib);
    K = H + d3;
tasklib_toc(tasklib, 'matmul3');

    d4 = dt*rhs_tic(t+dt,K,tasklib);

tasklib_tic(tasklib);
    H = H + 1/6*(d1 + 2*d2 + 2*d3 + d4);
tasklib_toc(tasklib, 'step');
end
t = tspan(end);
H = H(:).'; % Must return row vector

