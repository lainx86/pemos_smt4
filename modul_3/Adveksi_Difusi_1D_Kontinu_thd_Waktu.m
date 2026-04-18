clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi persamaan adveksi-difusi 1 dimensi dengan metode FTCS
% Sumber grid polutan kontinu
% Grafik terhadap ruang

x_digit = 0;
y_digit = 8;
z_digit = 7;

C_source = 10 * z_digit;

L = 3000;
T = 7200;

m_kontinu = 15 + z_digit;

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 30;
scenario(1).dt   = 6;
scenario(1).u    = 0.5;
scenario(1).D    = 4.0;
scenario(1).m    = m_kontinu;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 30;
scenario(2).dt   = 6;
scenario(2).u    = 0.1;
scenario(2).D    = 2.0;
scenario(2).m    = m_kontinu;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 30;
scenario(3).dt   = 6;
scenario(3).u    = 0.5;
scenario(3).D    = 2.0;
scenario(3).m    = m_kontinu;

scenario(4).nama = 'Skenario 4';
scenario(4).dx   = 30;
scenario(4).dt   = 6;
scenario(4).u    = 0.1;
scenario(4).D    = 4.0;
scenario(4).m    = m_kontinu;

folder_output = 'output_Kontinu_FTCS_Adveksi_Difusi_Terhadap_Waktu';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\n==================================================\n');
    fprintf('Memproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    u        = scenario(s).u;
    D        = scenario(s).D;
    m_source = scenario(s).m;

    x_space = 0:dx:L;
    t_grid  = 0:dt:T;

    Mmax = length(x_space);
    Nmax = length(t_grid) - 1;

    courant = u * dt / dx;
    alpha   = D * dt / dx^2;

    fprintf('dx = %.3f m, dt = %.3f s, u = %.3f m/s, D = %.3f m^2/s, m = %d\n', ...
        dx, dt, u, D, m_source);
    fprintf('Mmax = %d, Nmax = %d\n', Mmax, Nmax);
    fprintf('Courant = %.6f, alpha = %.6f\n', courant, alpha);

    if m_source < 2 || m_source > Mmax - 1
        fprintf('%s tidak dapat diproses karena titik sumber terlalu dekat dengan batas/domain tidak valid.\n', ...
            scenario(s).nama);
        continue;
    end

    stabil_diffusi = (alpha <= 0.5);
    stabil_gabungan = (courant^2 <= 2 * alpha);

    if ~stabil_diffusi
        fprintf(['PERINGATAN: alpha = %.6f > 0.5, sehingga skema eksplisit ', ...
                 'untuk komponen difusi berpotensi tidak stabil.\n'], alpha);
    end

    if ~stabil_gabungan
        fprintf(['PERINGATAN: Courant^2 = %.6f > 2*alpha = %.6f, ', ...
                 'sehingga kestabilan gabungan adveksi-difusi berpotensi tidak terpenuhi.\n'], ...
                 courant^2, 2*alpha);
    end

    if stabil_diffusi && stabil_gabungan
        fprintf('Syarat kestabilan numerik utama terpenuhi.\n');
    end

    C_old = zeros(1, Mmax);
    C_new = zeros(1, Mmax);

    idx_grid_plot = round(linspace(2, Mmax - 1, 5));
    x_plot = x_space(idx_grid_plot);

    snapshot = NaN(length(idx_grid_plot), Nmax + 1);

    S = zeros(1, Mmax);
    S(m_source) = C_source / dt;

    snapshot(:, 1) = C_old(idx_grid_plot);

    for n = 1:Nmax

        C_old(1)    = C_old(2);
        C_old(Mmax) = C_old(Mmax - 1);

        C_new = zeros(1, Mmax);

        for i = 2:(Mmax - 1)
            C_new(i) = C_old(i) ...
                     - (courant / 2) * (C_old(i + 1) - C_old(i - 1)) ...
                     + alpha * (C_old(i + 1) - 2 * C_old(i) + C_old(i - 1)) ...
                     + dt * S(i);
        end

        C_new(1)    = C_new(2);
        C_new(Mmax) = C_new(Mmax - 1);

        snapshot(:, n + 1) = C_new(idx_grid_plot);

        C_old = C_new;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(idx_grid_plot)
        plot(t_grid, snapshot(k, :), 'LineWidth', 1.2);
        legend_text{end + 1} = sprintf('x = %.0f m', x_plot(k));
    end

    grid on;
    xlabel('Waktu t (s)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap waktu', scenario(s).nama), ...
        'FontSize', 14, 'FontWeight', 'bold');

    xlim([t_grid(1), t_grid(end)]);

    y_max = max(snapshot(:), [], 'omitnan');
    if isempty(y_max) || isnan(y_max) || y_max <= 0
        ylim([0 1]);
    else
        ylim([0, y_max * 1.1]);
    end

    set(gca, 'FontSize', 11);

    if ~isempty(legend_text)
        legend(legend_text, 'Location', 'northwest', ...
            'FontSize', 8, 'NumColumns', 2, 'Box', 'on');
    end

    hold off;

    nama_file = sprintf('ftcs_adveksi_difusi_waktu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');
