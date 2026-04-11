clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi difusi 1 dimensi metode Crank-Nicolson
% Sumber polutan kontinu
% Grafik: konsentrasi polutan terhadap ruang
% NIM = 087, sehingga:
% x = 0, y = 8, z = 7

x_digit = 0;
y_digit = 8;
z_digit = 7;

C_source = 10 * z_digit;

L = 3000;
T = 7200;

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 15;
scenario(1).dt   = 6;
scenario(1).G    = 0.5;
scenario(1).m    = 15 + z_digit;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 30;
scenario(2).dt   = 6;
scenario(2).G    = 0.05;
scenario(2).m    = 30 + z_digit;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 15;
scenario(3).dt   = 6;
scenario(3).G    = 1;
scenario(3).m    = 15 + z_digit;

folder_output = 'output_Kontinu_CN_Polutan_Terhadap_Ruang';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    G        = scenario(s).G;
    m_source = scenario(s).m;

    M_grid = round(L / dx) + 1;
    Mmax   = max(M_grid, m_source);
    Nmax   = round(T / dt);

    alpha = (G * dt) / (dx^2);

    fprintf('dx = %.3f m, dt = %.3f s, G = %.3f m^2/s, m = %d\n', dx, dt, G, m_source);
    fprintf('Mmax = %d, Nmax = %d, alpha = %.6f\n', Mmax, Nmax, alpha);

    F = zeros(Nmax + 1, Mmax);

    t_step = 0:Nmax;
    t_time = t_step * dt;
    x_grid = 0:(Mmax - 1);
    x_pos  = x_grid * dx;

    tv = [60 180 360 500 800 1160];
    tv_valid = tv(tv >= 0 & tv <= T);
    n_plot = round(tv_valid / dt) + 1;

    if length(tv_valid) < length(tv)
        fprintf('%s: beberapa waktu berada di luar domain sehingga tidak digrafikkan.\n', scenario(s).nama);
    end

    a = -alpha;
    b = 2 + 2 * alpha;
    c = -alpha;

    H = zeros(1, Mmax);
    Q = zeros(1, Mmax);
    D = zeros(1, Mmax);

    for n = 1:Nmax

        F(n, m_source) = C_source;

        for i = 2:(Mmax-1)
            D(i) = alpha * F(n, i+1) + (2 - 2 * alpha) * F(n, i) + alpha * F(n, i-1);
        end

        D(2)      = D(2) + alpha * F(n+1, 1);
        D(Mmax-1) = D(Mmax-1) + alpha * F(n+1, Mmax);

        H(2) = c / b;
        Q(2) = D(2) / b;

        for i = 3:(Mmax-1)
            denom = b - a * H(i-1);
            H(i) = c / denom;
            Q(i) = (D(i) - a * Q(i-1)) / denom;
        end

        F(n+1, Mmax-1) = Q(Mmax-1);

        for j = (Mmax-2):-1:2
            F(n+1, j) = H(j) * F(n+1, j+1) + Q(j);
        end

        F(n+1, 1)    = F(n+1, 2);
        F(n+1, Mmax) = F(n+1, Mmax-1);
        F(n+1, m_source) = C_source;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(n_plot)
        plot(x_pos, F(n_plot(k), :), 'LineWidth', 1.2);
        legend_text{end+1} = sprintf('t = %d s', tv_valid(k));
    end

    grid on;
    xlabel('Ruang x (meter)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap ruang', scenario(s).nama), ...
        'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([x_pos(1), x_pos(end)]);

    y_max = max(F(n_plot, :), [], 'all');
    if isempty(y_max) || isnan(y_max) || y_max <= 0
        ylim([0 1]);
    else
        ylim([0, max(C_source * 1.2, y_max * 1.1)]);
    end

    set(gca, 'FontSize', 11);

    legend(legend_text, ...
        'Location', 'northwest', ...
        'FontSize', 8, ...
        'NumColumns', 2, ...
        'Box', 'on');

    hold off;

    nama_file = sprintf('cn_kontinu_ruang_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');