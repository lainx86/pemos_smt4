clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C

% Simulasi adveksi 1 dimensi dengan metode FTCS
% Sumber polutan diskontinu
% NIM ganjil, xyz = 087, sehingga z = 7 dan 0.zyx = 0.780

x_digit = 0;
y_digit = 8;
z_digit = 7;

u_xyz    = 0.780;
C_source = 10 * z_digit;

L = 3000;
T = 10800;

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 30;
scenario(1).dt   = 8;
scenario(1).u    = u_xyz;
scenario(1).n    = 20 + z_digit;
scenario(1).m    = 30 + z_digit;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 60;
scenario(2).dt   = 8;
scenario(2).u    = u_xyz;
scenario(2).n    = 20 + z_digit;
scenario(2).m    = 60 + z_digit;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 15;
scenario(3).dt   = 8;
scenario(3).u    = u_xyz;
scenario(3).n    = 20 + z_digit;
scenario(3).m    = 15 + z_digit;

scenario(4).nama = 'Skenario 4';
scenario(4).dx   = 30;
scenario(4).dt   = 8;
scenario(4).u    = 0.5 * u_xyz;
scenario(4).n    = 20 + z_digit;
scenario(4).m    = 30 + z_digit;

scenario(5).nama = 'Skenario 5';
scenario(5).dx   = 30;
scenario(5).dt   = 8;
scenario(5).u    = 2 * u_xyz;
scenario(5).n    = 20 + z_digit;
scenario(5).m    = 30 + z_digit;

scenario(6).nama = 'Skenario 6';
scenario(6).dx   = 30;
scenario(6).dt   = 8;
scenario(6).u    = -u_xyz;
scenario(6).n    = 20 + z_digit;
scenario(6).m    = 30 + z_digit;

scenario(7).nama = 'Skenario 7';
scenario(7).dx   = 30;
scenario(7).dt   = 2;
scenario(7).u    = u_xyz;
scenario(7).n    = 100 + z_digit;
scenario(7).m    = 30 + z_digit;

scenario(8).nama = 'Skenario 8';
scenario(8).dx   = 30;
scenario(8).dt   = 16;
scenario(8).u    = u_xyz;
scenario(8).n    = 10 + z_digit;
scenario(8).m    = 30 + z_digit;

folder_output = 'output_Diskontinu_FTCS_Polutan_Terhadap_Waktu';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    u        = scenario(s).u;
    n_source = scenario(s).n;
    m_source = scenario(s).m;

    M_grid = round(L / dx);
    Mmax   = max(M_grid, m_source);
    Nmax   = round(T / dt);

    courant = u * dt / dx;

    fprintf('dx = %.3f m, dt = %.3f s, u = %.3f m/s, n = %d, m = %d\n', ...
            dx, dt, u, n_source, m_source);
    fprintf('Mmax = %d, Nmax = %d, C = %.4f\n', Mmax, Nmax, courant);

    % F0 = kondisi pada waktu sekarang (n)
    % F  = kondisi pada waktu berikutnya (n+1)
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

    % Kondisi awal n = 0
    F0(:) = 0;

    % Sumber diskontinu pada n = 0
    if n_source == 0
        F0(m_source) = C_source;
    end

    for p = 1:length(m_plot_valid)
        history(1, p) = F0(m_plot_valid(p));
    end

    % Loop waktu FTCS
    for J = 1:Nmax

        % Boundary condition
        F0(1)    = F0(2);
        F0(Mmax) = F0(Mmax-1);

        % Skema FTCS
        for i = 2:(Mmax-1)
            F(i) = F0(i) - (u * dt / (2 * dx)) * (F0(i+1) - F0(i-1));
        end

        % Boundary condition untuk hasil baru
        F(1)    = F(2);
        F(Mmax) = F(Mmax-1);

        % Sumber diskontinu hanya diberikan sekali pada (n_source, m_source)
        if J == n_source
            F(m_source) = C_source;
        end

        % Simpan riwayat waktu
        for p = 1:length(m_plot_valid)
            history(J + 1, p) = F(m_plot_valid(p));
        end

        % Update ke waktu berikutnya
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
    title(sprintf('%s - Konsentrasi polutan terhadap waktu (FTCS)', scenario(s).nama), ...
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

    nama_file = sprintf('ftcs_diskontinu_waktu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');