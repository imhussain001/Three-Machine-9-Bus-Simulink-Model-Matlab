%WSCC9_LFC_PLOT  Simulate wscc9_lfc_governor.slx and plot the load-frequency response.
%   Figure 0  Overview: frequency of all SGs in one plot + RoCoF, dPm, dPe, angles, key numbers
%   Figure 1  Frequency: f of each SG (full run and zoom on the event), RoCoF
%   Figure 2  Power: Pm, Pe, governor response dPm, accelerating power Pm - Pe
%   Figure 3  Rotor angles relative to SG1
%   Command window: per-SG table of nadir, RoCoF, settling, steady state, dPm.
%   Parameters (load step size/bus/time, R, Tg, Tt, D) are in wscc9_lfc_init.m.

out = sim('wscc9_lfc_governor');            % InitFcn runs wscc9_lfc_init in the base workspace

to_cols = @(ts) reshape(squeeze(ts.Data), numel(ts.Time), []);   % N x 3
t     = out.f_SG_Hz.Time;
f_SG  = to_cols(out.f_SG_Hz);               % [Hz]
Pm_MW = to_cols(out.Pm_SG_MW);              % [MW]
Pe_MW = to_cols(out.Pe_SG_MW);              % [MW]
dl    = to_cols(out.delta_SG_rad);          % [rad]

% RoCoF from the swing equation (exact, no numerical differentiation):
% d(dw)/dt = (pi f / H) (Pm - Pe - D kS dw/ws)  ->  RoCoF [Hz/s] = d(dw)/dt / (2 pi)
dw_pu = (f_SG - f) / f;
rocof = (f ./ (2*H')) .* ((Pm_MW - Pe_MW)/S_base - D * kS' .* dw_pu);   % [Hz/s]

% Reference values
after   = t >= t_step;
k0      = find(t < t_step, 1, 'last');      % last sample before the step
dP_load = sum(Pe_MW(end,:) - Pe_MW(k0,:));  % actual extra electrical load incl. losses [MW]
df_ss   = -(dP_load/S_base) / beta * f;     % predicted steady-state deviation [Hz]
f_ss    = f + df_ss;
band    = 0.05 * abs(df_ss);                % settling band: 5 % of the steady-state deviation [Hz]

names  = {'SG1 (247.5 MVA)', 'SG2 (192 MVA)', 'SG3 (128 MVA)'};
colors = lines(3);
t_zoom = [t_step - 0.5, t_step + 10];

%% Figure 1: frequency
fig1 = figure('Name', 'Frequency of each SG', 'Color', 'w');
tl1  = tiledlayout(fig1, 3, 1, 'TileSpacing', 'compact');
title(tl1, sprintf('Load step %+g MW at bus %d, t = %g s  (R = %g, T_g = %g s, T_t = %g s, D = %g)', ...
    dP_L_MW, step_bus, t_step, R, Tg, Tt, D));

ax = nexttile(tl1); hold(ax, 'on');
for k = 1:3, plot(ax, t, f_SG(:,k), 'Color', colors(k,:), 'LineWidth', 1.1); end
yline(ax, f_ss, '--k', sprintf('predicted %.4f Hz', f_ss), 'LabelHorizontalAlignment', 'left');
ylabel(ax, 'Frequency [Hz]'); title(ax, 'Full run'); grid(ax, 'on');
legend(ax, names, 'Location', 'southeast');

ax = nexttile(tl1); hold(ax, 'on');
for k = 1:3
    plot(ax, t, f_SG(:,k), 'Color', colors(k,:), 'LineWidth', 1.1);
    [fn, in] = min(f_SG(after,k)); ta = t(after);
    plot(ax, ta(in), fn, 'v', 'Color', colors(k,:), 'MarkerFaceColor', colors(k,:));
end
yline(ax, f_ss, '--k');
xlim(ax, t_zoom); ylabel(ax, 'Frequency [Hz]'); title(ax, 'Zoom on the event (markers: nadir)'); grid(ax, 'on');

ax = nexttile(tl1); hold(ax, 'on');
for k = 1:3, plot(ax, t, rocof(:,k), 'Color', colors(k,:), 'LineWidth', 1.1); end
xlim(ax, t_zoom); ylabel(ax, 'RoCoF [Hz/s]'); xlabel(ax, 'Time [s]'); title(ax, 'Rate of change of frequency'); grid(ax, 'on');

%% Figure 2: power
fig2 = figure('Name', 'Power of each SG', 'Color', 'w');
tl2  = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact');
title(tl2, 'Mechanical and electrical power');
dPm_MW = Pm_MW - Pm_MW(k0,:);
dPm_expected = -(kS/R) * (df_ss/f) * S_base;      % droop share at the steady state [MW], 3 x 1
panels = {Pm_MW, 'P_m [MW]', 'Mechanical power (governor + turbine)'
          Pe_MW, 'P_e [MW]', 'Electrical power'
          dPm_MW, '\DeltaP_m [MW]', 'Governor response (dotted: droop share)'
          Pm_MW - Pe_MW - D*kS'.*dw_pu*S_base, 'P_a [MW]', 'Accelerating power P_m - P_e - D k_S \Delta\omega'};
for p = 1:4
    ax = nexttile(tl2); hold(ax, 'on');
    for k = 1:3, plot(ax, t, panels{p,1}(:,k), 'Color', colors(k,:), 'LineWidth', 1.1); end
    if p == 3
        for k = 1:3, yline(ax, dPm_expected(k), ':', 'Color', colors(k,:), 'LineWidth', 1.2); end
    end
    ylabel(ax, panels{p,2}); title(ax, panels{p,3}); grid(ax, 'on'); xlabel(ax, 'Time [s]');
    if p == 4, xlim(ax, t_zoom); end
    if p == 1, legend(ax, names, 'Location', 'east'); end
end

%% Figure 3: rotor angles
fig3 = figure('Name', 'Rotor angles', 'Color', 'w');
ax = axes(fig3); hold(ax, 'on');
plot(ax, t, rad2deg(dl(:,2) - dl(:,1)), 'Color', colors(2,:), 'LineWidth', 1.1);
plot(ax, t, rad2deg(dl(:,3) - dl(:,1)), 'Color', colors(3,:), 'LineWidth', 1.1);
ylabel(ax, 'Angle [deg]'); xlabel(ax, 'Time [s]'); grid(ax, 'on');
title(ax, 'Rotor angles relative to SG1'); legend(ax, {'\delta_{21}', '\delta_{31}'}, 'Location', 'east');

%% Summary table
nadir = zeros(3,1); t_nadir = nadir; rocof_max = nadir; t_settle = nadir;
for k = 1:3
    fa = f_SG(after,k); ta = t(after);
    [nadir(k), in] = min(fa); t_nadir(k) = ta(in);
    [~, ir] = max(abs(rocof(after,k))); ra = rocof(after,k); rocof_max(k) = ra(ir);
    outside = find(abs(fa - f_SG(end,k)) > band, 1, 'last');
    if isempty(outside), t_settle(k) = 0; else, t_settle(k) = ta(min(outside + 1, numel(ta))) - t_step; end
end
summary = table(nadir, t_nadir, rocof_max, t_settle, f_SG(end,:)', dPm_MW(end,:)', dPm_expected, ...
    (Pe_MW(end,:) - Pe_MW(k0,:))', 'RowNames', {'SG1', 'SG2', 'SG3'}, 'VariableNames', ...
    {'Nadir_Hz', 't_nadir_s', 'RoCoF_max_Hz_per_s', 'Settling_s', 'f_final_Hz', 'dPm_MW', 'dPm_droop_MW', 'dPe_MW'});
fprintf('\nLoad step %+g MW at bus %d (t = %g s); actual extra load incl. losses %.2f MW\n', ...
    dP_L_MW, step_bus, t_step, dP_load);
fprintf('Predicted steady state %.4f Hz (df = %.1f mHz); settling band +/- %.2f mHz\n\n', ...
    f_ss, 1e3*df_ss, 1e3*band);
disp(summary);

%% Figure 0: one-page overview (frequency of every SG + the other key quantities)
fig0 = figure('Name', 'Load-frequency overview', 'Color', 'w');
tl0  = tiledlayout(fig0, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl0, sprintf('WSCC 9-bus: frequency response to a %+g MW load step at bus %d, t = %g s', ...
    dP_L_MW, step_bus, t_step), 'FontWeight', 'bold');
t_view = [max(0, t_step - 1), min(t(end), t_step + 30)];

ax = nexttile(tl0, 1, [2 2]); hold(ax, 'on');
hf = gobjects(3, 1);
for k = 1:3
    hf(k) = plot(ax, t, f_SG(:,k), 'Color', colors(k,:), 'LineWidth', 1.4, 'DisplayName', ...
        sprintf('%s, nadir %.4f Hz at %.2f s', names{k}, nadir(k), t_nadir(k)));
    plot(ax, t_nadir(k), nadir(k), 'v', 'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), 'MarkerSize', 7);
end
yline(ax, f_ss, '--k', sprintf('predicted steady state %.4f Hz', f_ss), ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
xlim(ax, t_view); grid(ax, 'on'); xlabel(ax, 'Time [s]'); ylabel(ax, 'Frequency [Hz]');
title(ax, 'Frequency of each SG (markers: nadir)');
legend(ax, hf, 'Location', 'southeast');

ax = nexttile(tl0, 3); hold(ax, 'on');
for k = 1:3, plot(ax, t, rocof(:,k), 'Color', colors(k,:), 'LineWidth', 1.1); end
xlim(ax, [t_view(1), t_step + 8]); grid(ax, 'on'); ylabel(ax, 'RoCoF [Hz/s]'); title(ax, 'Rate of change of frequency');

ax = nexttile(tl0, 6); hold(ax, 'on');
for k = 1:3
    plot(ax, t, dPm_MW(:,k), 'Color', colors(k,:), 'LineWidth', 1.1);
    yline(ax, dPm_expected(k), ':', 'Color', colors(k,:), 'LineWidth', 1.2);
end
xlim(ax, t_view); grid(ax, 'on'); ylabel(ax, '\DeltaP_m [MW]'); title(ax, 'Governor response (dotted: droop share)');

ax = nexttile(tl0, 7); hold(ax, 'on');
for k = 1:3, plot(ax, t, Pe_MW(:,k) - Pe_MW(k0,k), 'Color', colors(k,:), 'LineWidth', 1.1); end
xlim(ax, t_view); grid(ax, 'on'); xlabel(ax, 'Time [s]'); ylabel(ax, '\DeltaP_e [MW]');
title(ax, 'Change of electrical power');

ax = nexttile(tl0, 8); hold(ax, 'on');
plot(ax, t, rad2deg(dl(:,2) - dl(:,1)), 'Color', colors(2,:), 'LineWidth', 1.1);
plot(ax, t, rad2deg(dl(:,3) - dl(:,1)), 'Color', colors(3,:), 'LineWidth', 1.1);
xlim(ax, t_view); grid(ax, 'on'); xlabel(ax, 'Time [s]'); ylabel(ax, 'Angle [deg]');
title(ax, 'Rotor angles relative to SG1'); legend(ax, {'\delta_{21}', '\delta_{31}'}, 'Location', 'best');

ax = nexttile(tl0, 9); axis(ax, 'off');
lines_txt = {
    sprintf('Load step %+g MW at bus %d, t = %g s', dP_L_MW, step_bus, t_step)
    sprintf('Actual extra load (incl. losses) %.2f MW', dP_load)
    sprintf('R = %g, T_g = %g s, T_t = %g s, D = %g', R, Tg, Tt, D)
    ''
    '        nadir [Hz]  t [s]  RoCoF [Hz/s]  t_{set} [s]  \DeltaP_m [MW]'};
for k = 1:3
    lines_txt{end+1} = sprintf('SG%d    %9.4f  %5.2f  %10.3f  %9.2f  %10.2f', ...
        k, nadir(k), t_nadir(k), rocof_max(k), t_settle(k), dPm_MW(end,k)); %#ok<SAGROW>
end
lines_txt{end+1} = '';
lines_txt{end+1} = sprintf('Final f: %.4f Hz (predicted %.4f Hz)', mean(f_SG(end,:)), f_ss);
text(ax, 0, 1, lines_txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Consolas', 'FontSize', 9);
title(ax, 'Key numbers');

%% Save the figures as PNG in .\figures, named after the case
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
case_tag = sprintf('dP%+gMW_bus%d_t%gs', dP_L_MW, step_bus, t_step);
figs = {fig0, 'overview', [30 30 1500 900]; fig1, 'frequency', [50 50 1000 850]
        fig2, 'power', [50 50 1100 700]; fig3, 'rotor_angles', [50 50 800 350]};
for k = 1:size(figs, 1)
    set(figs{k,1}, 'Position', figs{k,3});
    png = fullfile(fig_dir, sprintf('wscc9_lfc_%s_%s.png', figs{k,2}, case_tag));
    try
        exportgraphics(figs{k,1}, png, 'Resolution', 200);
        fprintf('Saved %s\n', png);
    catch err
        warning('wscc9_lfc_plot:save', 'Could not write %s (close it if it is open in another program): %s', ...
            png, err.message);
    end
end
