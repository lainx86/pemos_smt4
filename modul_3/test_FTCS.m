clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi adveksi-difusi 1 dimensi dengan metode FTCS
% Sumber polutan kontinu
% Grafik: konsentrasi polutan terhadap waktu
% NIM = 087, sehingga:
% x = 0, y = 8, z = 7
% u_xyz = 0.780 m/s
% Konsentrasi sumber = 10*z = 70 mg/L

x_digit = 0;
y_digit = 8;
z_digit = 7;

u_xyz    = 0.780;
C_source = 10 * z_digit;

L = 3000;   % panjang domain (m)
T = 10800;  % waktu simulasi (s)

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 30;
scenario(1).dt   = 6;
scenario(1).u    = u_xyz;
scenario(1).D    = 0.05;
scenario(1).m    = 30 + z_digit;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 60;
scenario(2).dt   = 6;
scenario(2).u    = u_xyz;
scenario(2).D    = 0.05;
scenario(2).m    = 60 + z_digit;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 15;
scenario(3).dt   = 6;
scenario(3).u    = u_xyz;
scenario(3).D    = 0.05;
scenario(3).m    = 15 + z_digit;

scenario(4).nama = 'Skenario 4';
scenario(4).dx   = 30;
scenario(4).dt   = 6;
scenario(4).u    = 0.5 * u_xyz;
scenario(4).D    = 0.05;
scenario(4).m    = 30 + z_digit;

scenario(5).nama = 'Skenario 5';
scenario(5).dx   = 30;
scenario(5).dt   = 6;
scenario(5).u    = 2 * u_xyz;
scenario(5).D    = 0.05;
scenario(5).m    = 30 + z_digit;

scenario(6).nama = 'Skenario 6';
scenario(6).dx   = 30;
scenario(6).dt   = 6;
scenario(6).u    = -u_xyz;
scenario(6).D    = 0.05;
scenario(6).m    = 30 + z_digit;

scenario(7).nama = 'Skenario 7';
scenario(7).dx   = 30;
scenario(7).dt   = 3;
scenario(7).u    = u_xyz;
scenario(7).D    = 0.05;
scenario(7).m    = 30 + z_digit;

scenario(8).nama = 'Skenario 8';
scenario(8).dx   = 30;
scenario(8).dt   = 12;
scenario(8).u    = u_xyz;
scenario(8).D    = 0.05;
scenario(8).m    = 30 + z_digit;

scenario(9).nama = 'Skenario 9';
scenario(9).dx   = 30;
scenario(9).dt   = 6;
scenario(9).u    = u_xyz;
scenario(9).D    = 0.50;
scenario(9).m    = 30 + z_digit;

scenario(10).nama = 'Skenario 10';
scenario(10).dx   = 30;
scenario(10).dt   = 6;
scenario(10).u    = u_xyz;
scenario(10).D    = 1.00;
scenario(10).m    = 30 + z_digit;

folder_output = 'output_Kontinu_FTCS_Adveksi_Difusi_Terhadap_Waktu';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    u        = scenario(s).u;
    D        = scenario(s).D;
    m_source = scenario(s).m;

    x_grid = 0:dx:L;
    t_grid = 0:dt:T;

    Mmax = length(x_grid);
    Nmax = length(t_grid);

    courant = u * dt / dx;
    alpha   = D * dt / dx^2;

    fprintf('dx = %.3f m, dt = %.3f s, u = %.3f m/s, D = %.3f m^2/s, m = %d\n', ...
            dx, dt, u, D, m_source);
    fprintf('Mmax = %d, Nmax = %d, Courant = %.4f, alpha = %.6f\n', ...
            Mmax, Nmax, courant, alpha);

    if m_source > Mmax
        fprintf('%s tidak dapat diproses karena grid sumber di luar domain.\n', scenario(s).nama);
        continue;
    end

    if alpha > 0.5
        fprintf('%s: peringatan, alpha > 0.5 sehingga skema FTCS berpotensi tidak stabil.\n', scenario(s).nama);
    end

    F0 = zeros(1, Mmax);
    F  = zeros(1, Mmax);

    m_plot = [1, ...
              m_source - 20, m_source - 10, m_source - 5, m_source - 3, m_source - 1, ...
              m_source, ...
              m_source + 1, m_source + 3, m_source + 5, m_source + 10, m_source + 20];

    m_plot = unique(m_plot, 'stable');
    m_plot_valid = m_plot(m_plot >= 1 & m_plot <= Mmax);

    if length(m_plot_valid) < length(m_plot)
        fprintf('%s: beberapa grid berada di luar domain sehingga tidak digrafikkan.\n', scenario(s).nama);
    end

    history = NaN(Nmax, length(m_plot_valid));

    F0(:) = 0;
    F0(m_source) = C_source;

    for p = 1:length(m_plot_valid)
        history(1, p) = F0(m_plot_valid(p));
    end

    for n = 1:(Nmax - 1)

        % kondisi batas Neumann sederhana
        F0(1)    = F0(2);
        F0(Mmax) = F0(Mmax - 1);

        % sumber kontinu
        F0(m_source) = C_source;

        for i = 2:(Mmax - 1)
            if i == m_source
                F(i) = C_source;
            else
                F(i) = F0(i) ...
                     - (u * dt / (2 * dx)) * (F0(i + 1) - F0(i - 1)) ...
                     + (D * dt / dx^2) * (F0(i + 1) - 2 * F0(i) + F0(i - 1));
            end
        end

        % kondisi batas sesudah update
        F(1)    = F(2);
        F(Mmax) = F(Mmax - 1);

        % sumber tetap kontinu setiap waktu
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
        plot(t_grid, history(:, k), 'LineWidth', 1.2);
        legend_text{end + 1} = sprintf('m = %d', m_plot_valid(k));
    end

    grid on;
    xlabel('Waktu t (detik)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap waktu', scenario(s).nama), ...
          'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([t_grid(1), t_grid(end)]);

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

    nama_file = sprintf('ftcs_adveksi_difusi_kontinu_waktu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');
