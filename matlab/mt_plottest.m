clear;
par = struct('test', 'galew', 'n', 6400, 'fd', 31, 'ep', 2.7, 'order', 4, 'gamma_c', -0.05);
mt_plot(par, getsol(par, 600));
figure
par = struct('test', 'tc5', 'n', 6400, 'fd', 31, 'ep', 2.7, 'order', 4, 'gamma_c', -0.05);
mt_plot(par, getsol(par, 900));
