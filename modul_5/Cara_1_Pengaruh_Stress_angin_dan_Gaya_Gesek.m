clear; clc; close all;

identitas = 'Feby Syarief Al A`raaf - 26050124130087 - Oseanografi C';

folder_output = 'Cara_1_Grafik pengaruh stress angin dan gaya gesek dasar';
if ~exist(folder_output, 'dir'); mkdir(folder_output); end

g = 9.81;

% NIM = 26050124130087
% 3 digit terakhir NIM = 087
% SKENARIO (NIM Ganjil: 087 -> A = 1.780)

scenario(1).nama  = 'Skenario 1';
scenario(1).L     = 6000;
scenario(1).dx    = 30;
scenario(1).LS    = 2250;
scenario(1).dt    = 5;
scenario(1).A     = 1.780;
scenario(1).T     = 450;
scenario(1).H     = 3;
scenario(1).lamCD = 3.2e-6;
scenario(1).W     = 5;
scenario(1).r     = 0.02;

scenario(2).nama  = 'Skenario 2';
scenario(2).L     = 6000;
scenario(2).dx    = 60;
scenario(2).LS    = 2250;
scenario(2).dt    = 5;
scenario(2).A     = 1.780;
scenario(2).T     = 450;
scenario(2).H     = 3;
scenario(2).lamCD = 3.2e-6;
scenario(2).W     = -5;
scenario(2).r     = 0.02;

scenario(3).nama  = 'Skenario 3';
scenario(3).L     = 6000;
scenario(3).dx    = 15;
scenario(3).LS    = 2250;
scenario(3).dt    = 5;
scenario(3).A     = 1.780;
scenario(3).T     = 450;
scenario(3).H     = 3;
scenario(3).lamCD = 3.2e-6;
scenario(3).W     = 10;
scenario(3).r     = 0.02;

scenario(4).nama  = 'Skenario 4';
scenario(4).L     = 6000;
scenario(4).dx    = 30;
scenario(4).LS    = 2250;
scenario(4).dt    = 2;
scenario(4).A     = 1.780;
scenario(4).T     = 450;
scenario(4).H     = 3;
scenario(4).lamCD = 3.2e-6;
scenario(4).W     = 5;
scenario(4).r     = 0.005;

scenario(5).nama  = 'Skenario 5';
scenario(5).L     = 6000;
scenario(5).dx    = 30;
scenario(5).LS    = 2250;
scenario(5).dt    = 2;
scenario(5).A     = 1.780;
scenario(5).T     = 450;
scenario(5).H     = 3;
scenario(5).lamCD = 3.2e-6;
scenario(5).W     = 5;
scenario(5).r     = 0.05;

for s = 1:length(scenario)

    L     = scenario(s).L;
    dx    = scenario(s).dx;
    LS    = scenario(s).LS;
    dt    = scenario(s).dt;
    A     = scenario(s).A;
    T     = scenario(s).T;
    H     = scenario(s).H;
    lamCD = scenario(s).lamCD;
    W     = scenario(s).W;
    r     = scenario(s).r;

    x  = 0:dx:L;
    t  = 0:dt:LS;
    Nx = length(x);
    Nt = length(t);

    Z  = zeros(Nx, Nt);
    u  = zeros(Nx, Nt);

    Cs    = sqrt(g * H);
    lam   = Cs * T;
    k     = 2*pi / lam;
    sigma = 2*pi / T;
    CFL   = Cs * dt / dx;

    tau_s = lamCD * W * abs(W);

    fprintf('%s | Cs=%.4f | lambda=%.2f | CFL=%.4f | W=%.1f | r=%.4f | tau_s=%.6f\n', ...
        scenario(s).nama, Cs, lam, CFL, W, r, tau_s);

    for j = 1:Nx
        Z(j,1) = A * cos(k * x(j));
        u(j,1) = (A/H) * Cs * cos(k * (x(j) + 0.5*dx));
    end

    for n = 1:Nt-1

        tn1 = t(n+1);

        for j = 1:Nx-1
            u(j,n+1) = u(j,n) * (1 - r*dt/H * abs(u(j,n))) ...
                       - g * dt/dx * (Z(j+1,n) - Z(j,n)) ...
                       + dt/H * tau_s;
        end
        u(Nx,n+1) = (A/H) * Cs * cos(k * L - sigma * tn1);

        for j = 2:Nx
            Z(j,n+1) = Z(j,n) - H * dt/dx * (u(j,n+1) - u(j-1,n+1));
        end

        Z(1,n+1) = A * cos(k * x(1) - sigma * tn1);

    end

    tau_b = r * u .* abs(u);

    j_plot = round(Nx/2);

    fig = figure('Visible','off', 'Units','centimeters', 'Position',[0 0 24 14]);
    ax = axes('Parent', fig);
    ax.Position = [0.10 0.13 0.82 0.68];
    plot(ax, t, tau_b(j_plot,:), 'b-', 'LineWidth', 1.5);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    yline(tau_s, 'r--', 'LineWidth', 1.2);
    xlabel(ax, 't (s)', 'FontSize', 11);
    ylabel(ax, '\tau_b', 'FontSize', 11);
    title(ax, {
        sprintf('\\tau_b vs Waktu - %s (j=%d, x=%.0f m)', ...
            scenario(s).nama, j_plot, x(j_plot))
        identitas
    }, 'FontSize', 10, 'FontWeight', 'normal');
    legend(ax, {'\tau_b (gesekan dasar)', '\tau_s (stress angin)'}, ...
        'Location','northeast', 'FontSize', 10);

    saveas(fig, fullfile(folder_output, ...
        sprintf('taub_vs_waktu_skenario_%d.png', s)));
    close(fig);

end

disp('Seluruh grafik tau_b selesai dibuat dan disimpan.');