clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi difusi/adveksi 1 dimensi metode FTCS
% Sumber polutan diskontinu
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

folder_output = 'output_Diskontinu_FTCS_Polutan_Terhadap_Ruang';
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

    if alpha > 0.5
        fprintf('Peringatan: alpha > 0.5, solusi FTCS dapat tidak stabil.\n');
    end

    F0 = zeros(1, Mmax);
    F  = zeros(1, Mmax);

    t_step = 0:Nmax;
    t_time = t_step * dt;
    x_grid = 0:(Mmax-1);
    x_pos  = x_grid * dx;

    t_plot = [0, 60, 180, 360, 500, 800, 1160];
    t_plot = unique(t_plot, 'stable');
    t_plot_valid = t_plot(t_plot >= 0 & t_plot <= T);

    if length(t_plot_valid) < length(t_plot)
        fprintf('%s: beberapa waktu berada di luar domain sehingga tidak digrafikkan.\n', scenario(s).nama);
    end

    n_plot = round(t_plot_valid / dt) + 1;

    profile = NaN(length(t_plot_valid), Mmax);

    F0(:) = 0;

    for q = 1:length(n_plot)
        if n_plot(q) == 1
            profile(q, :) = F0;
        end
    end

    for n = 1:Nmax

        F0(1)    = F0(2);
        F0(Mmax) = F0(Mmax-1);

        if n == n_source
            F0(m_source) = C_source;
        end

        for i = 2:(Mmax-1)
            F(i) = (1 - 2 * alpha) * F0(i) + alpha * (F0(i+1) + F0(i-1));
        end

        F(1)    = F(2);
        F(Mmax) = F(Mmax-1);

        idx_now = n + 1;
        for q = 1:length(n_plot)
            if idx_now == n_plot(q)
                profile(q, :) = F;
            end
        end

        F0 = F;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(t_plot_valid)
        plot(x_pos, profile(k, :), 'LineWidth', 1.2);
        legend_text{end+1} = sprintf('t = %d s', t_plot_valid(k));
    end

    grid on;
    xlabel('Ruang x (meter)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap ruang', scenario(s).nama), ...
          'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([x_pos(1), x_pos(end)]);

    y_max = max(profile(:), [], 'omitnan');
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

    nama_file = sprintf('ftcs_diskontinu_ruang_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');