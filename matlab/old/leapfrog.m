function [t,H] = leapfrog(rhs,tspan,H,gamma)
% rhs - rhs(t,y)
% tspan - time vector [t0, t0+dt, ...,t1]
% H - solution

% Bootstrap with RK4
[~,H1] = rk4(rhs,[tspan(1) tspan(2)],H);
H1 = H1(end,:).';
H0bar = H;

dt = tspan(2)-tspan(1); % Assume equal time steps (except last)

for i = 3:length(tspan)-1
    t = tspan(i-1);
    H = H0bar + 2*dt*rhs(t,H1);
    H0bar = H1 + gamma*(H - 2*H1 + H0bar);
    % Set the value for the Robert filter for the next time step.
    H1 = H;
end

if tspan(end)-tspan(end-1)<dt
    % End with RK4 if last time step is shorter
    [t,H] = rk4(rhs,[tspan(end-1) tspan(end)],H);
else
    H = H0bar + 2*dt*rhs(t,H1);
    H = H(:).';
    t = tspan(end);
end