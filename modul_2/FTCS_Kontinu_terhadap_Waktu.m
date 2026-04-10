clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi difusi/adveksi 1 dimensi metode FTCS
% Sumber polutan kontinu
% NIM = 087, sehingga:
% x = 0, y = 8, z = 7

x_digit = 0;
y_digit = 8;
z_digit = 7;

C_source = 10 * z_digit;   % konsentrasi sumber polutan = 70 mg/L

L = 3000;   % panjang kanal (m)
T = 7200;   % lama simulasi (detik)

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

folder_output = 'output_Kontinu_FTCS_Polutan_Terhadap_Waktu';
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

    if alpha > 0.5
        fprintf('Peringatan: alpha > 0.5, solusi FTCS dapat tidak stabil.\n');
    end

    F0 = zeros(1, Mmax);
    F  = zeros(1, Mmax);

    t_step = 0:Nmax;
    t_time = t_step * dt;

    m_plot = [1, ...
              m_source - 20, m_source - 10, m_source - 5, m_source - 3, m_source - 1, ...
              m_source, ...
              m_source + 1, m_source + 3, m_source + 5, m_source + 10, m_source + 20];

    m_plot = unique(m_plot, 'stable');
    m_plot_valid = m_plot(m_plot >= 1 & m_plot <= Mmax);

    if length(m_plot_valid) < length(m_plot)
        fprintf('%s: beberapa grid berada di luar domain sehingga tidak digrafikkan.\n', scenario(s).nama);
    end

    history = NaN(length(t_step), length(m_plot_valid));

    F0(:) = 0;

    for p = 1:length(m_plot_valid)
        history(1, p) = F0(m_plot_valid(p));
    end

    for n = 1:Nmax

        % syarat batas pada kondisi lama
        F0(1)    = F0(2);
        F0(Mmax) = F0(Mmax-1);

        % sumber polutan kontinu
        F0(m_source) = C_source;

        % skema FTCS sesuai flowchart
        for i = 2:(Mmax-1)
            F(i) = (1 - 2 * alpha) * F0(i) + alpha * (F0(i+1) + F0(i-1));
        end

        % syarat batas pada kondisi baru
        F(1)    = F(2);
        F(Mmax) = F(Mmax-1);

        % sumber polutan tetap kontinu
        F(m_source) = C_source;

        for p = 1:length(m_plot_valid)
            history(n + 1, p) = F(m_plot_valid(p));
        end

        F0 = F;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(m_plot_valid)
        plot(t_time, history(:, k), 'LineWidth', 1.2);
        legend_text{end+1} = sprintf('m = %d', m_plot_valid(k));
    end

    grid on;
    xlabel('Waktu t (detik)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap waktu', scenario(s).nama), ...
          'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([t_time(1), t_time(end)]);

    y_max = max(history(:), [], 'omitnan');
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

    nama_file = sprintf('ftcs_kontinu_waktu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');