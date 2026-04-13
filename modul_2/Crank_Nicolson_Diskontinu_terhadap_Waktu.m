clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi difusi 1 dimensi metode Crank-Nicolson
% Sumber polutan diskontinu
% Grafik: konsentrasi polutan terhadap waktu
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

folder_output = 'output_Diskontinu_CN_Polutan_Terhadap_Waktu';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    G        = scenario(s).G;
    m_source = scenario(s).m;
    n_source = 20 + z_digit;

    M_grid = round(L / dx) + 1;
    Mmax   = max(M_grid, m_source);
    Nmax   = round(T / dt);

    alpha = (G * dt) / (dx^2);

    fprintf('dx = %.3f m, dt = %.3f s, G = %.3f m^2/s, m = %d, n = %d\n', dx, dt, G, m_source, n_source);
    fprintf('Mmax = %d, Nmax = %d, alpha = %.6f\n', Mmax, Nmax, alpha);

    F = zeros(Nmax + 1, Mmax);

    t_step = 0:Nmax;
    t_time = t_step * dt;

    rv = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41];
    rv_valid = rv(rv >= 1 & rv <= Mmax);

    if length(rv_valid) < length(rv)
        fprintf('%s: beberapa grid berada di luar domain sehingga tidak digrafikkan.\n', scenario(s).nama);
    end

    a = -alpha;
    b = 2 + 2 * alpha;
    c = -alpha;

    H = zeros(1, Mmax);
    Q = zeros(1, Mmax);
    D = zeros(1, Mmax);

    for n = 1:Nmax

        if n == n_source
            F(n, m_source) = C_source;
        end

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
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(rv_valid)
        plot(t_time, F(:, rv_valid(k)), 'LineWidth', 1.2);
        legend_text{end+1} = sprintf('m = %d', rv_valid(k));
    end

    grid on;
    xlabel('Waktu t (detik)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap waktu', scenario(s).nama), ...
        'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([t_time(1), t_time(end)]);

    y_max = max(F(:, rv_valid), [], 'all');
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

    nama_file = sprintf('cn_diskontinu_waktu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');