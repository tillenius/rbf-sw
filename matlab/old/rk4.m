function [t,H] = rk4(rhs,tspan,H)
% rhs - rhs(t,y)
% tspan - time vector [t0, t0+dt, ...,t1]
% H - solution

for i = 2:length(tspan)
    dt = tspan(i)-tspan(i-1);
    t = tspan(i);
    K = H;
    d1 = dt*rhs(t,K);
    K = H + 0.5*d1;
    d2 = dt*rhs(t+0.5*dt,K);
    K = H + 0.5*d2;
    d3 = dt*rhs(t+0.5*dt,K);
    K = H + d3;
    d4 = dt*rhs(t+dt,K);

    H = H + 1/6*(d1 + 2*d2 + 2*d3 + d4);
end
t = tspan(end);
H = H(:).'; % Must return row vector