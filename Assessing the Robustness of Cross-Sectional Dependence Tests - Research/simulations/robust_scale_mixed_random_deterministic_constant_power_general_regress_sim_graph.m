
clear all
clc

% Set seed.

rng(36006963);

%% MONTE CARLO SETTINGS

% Time dimensions.

T_vals = [20 50 100];

% Cross-sectional dimensions.

N_vals = [20 50 100];

% Pervasiveness parameter:
% 0 <= a <= 1
%
% Number of deterministic-scale individuals:
% N_deterministic = floor(N^a)

a_vals = 0:0.05:1;

% Number of Monte Carlo replications.

reps = 5000;

% Significance level.

alpha = 0.05;

%% OUTPUT FOLDER

% Store exported figures in their own output subfolder.

output_folder = fullfile('output', 'deterministic');

% Create the folder if it does not already exist.

if ~exist(output_folder, 'dir')

    mkdir(output_folder);

end

%% PREALLOCATE EMPIRICAL SIZE STORAGE

% Dimensions:
% T x N x a

CD_size = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

CD_robust_size = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

CD_general_size = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

%% PREALLOCATE NUMBER OF DETERMINISTIC UNITS

% Dimensions:
% N x a
%
% N_deterministic depends only on N and a, not T.

N_deterministic_storage = ...
    zeros(length(N_vals), length(a_vals));

%% MONTE CARLO SIZE EXPERIMENT

for k = 1:length(a_vals)

    a = a_vals(k);

    for i = 1:length(T_vals)

        for j = 1:length(N_vals)

            T = T_vals(i);
            N = N_vals(j);

            % Size DGP (gamma_i = 0):
            % y_it = beta_0 + beta_1 x_it + u_it
            % u_it = gamma_i f_t + h_it epsilon_it
            %
            % Deterministic-scale units:
            % h_it = (3/T) floor(t/10)
            %
            % Remaining units:
            % h_it = 1
            %
            % Number of deterministic-scale units:
            % N_deterministic = floor(N^a)
            %
            % Identities of these units are randomly selected
            % inside every Monte Carlo replication.

            result = ...
                robust_scale_mixed_random_deterministic_constant_power_general_regress_mc( ...
                    N, T, reps, a, 'size', [0 0], alpha);

            % Store number of deterministic-scale units.

            N_deterministic_storage(j,k) = ...
                result.N_deterministic;

            % Store empirical sizes.

            CD_size(i,j,k) = ...
                result.cd_rejection_rate;

            CD_robust_size(i,j,k) = ...
                result.cd_robust_rejection_rate;

            CD_general_size(i,j,k) = ...
                result.cd_general_rejection_rate;

        end

    end

end

%% DISPLAY NUMBER OF DETERMINISTIC UNITS

N_deterministic_table = array2table( ...
    N_deterministic_storage, ...
    'VariableNames', strcat("a", string(a_vals)), ...
    'RowNames', strcat("N", string(N_vals)));

disp(' ')
disp('NUMBER OF DETERMINISTIC-SCALE UNITS')
disp(N_deterministic_table)

%% COMMON Y-AXIS RANGE

% Combine the empirical sizes from all three tests.

all_sizes = [ ...
    CD_size(:); ...
    CD_robust_size(:); ...
    CD_general_size(:); ...
    alpha];

% Find overall minimum and maximum.

y_min = min(all_sizes);
y_max = max(all_sizes);

% Add a small amount of space above and below the plotted values.

y_range = y_max - y_min;

y_lower = y_min - 0.05*y_range;
y_upper = y_max + 0.05*y_range;

%% PLOT EMPIRICAL SIZE AGAINST a

% One figure is produced for each (N,T) combination.

for i = 1:length(T_vals)

    for j = 1:length(N_vals)

        T = T_vals(i);
        N = N_vals(j);

        % Extract empirical sizes across all a values
        % for this particular (N,T) combination.

        cd_plot = squeeze(CD_size(i,j,:));
        cd_robust_plot = squeeze(CD_robust_size(i,j,:));
        cd_general_plot = squeeze(CD_general_size(i,j,:));

        % Create new figure.

        fig = figure;

        % Plot standard CD empirical size.

        plot(a_vals, cd_plot, '-o', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'CD');

        hold on

        % Plot robust CD empirical size.

        plot(a_vals, cd_robust_plot, '-s', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Robust CD');

        % Plot general CD empirical size.

        plot(a_vals, cd_general_plot, '-^', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'General CD');

        % Plot nominal significance level.

        yline(alpha, '--', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Nominal size');

        % Axis labels.

        xlabel('a')
        ylabel('Empirical size')

        % Figure title.

        title(['N = ', num2str(N), ...
               ', T = ', num2str(T)])

        % Display legend.

        legend('Location', 'best')

        % Add grid.

        grid on

        % Keep same horizontal and vertical range for every figure.

        xlim([0 1])
        ylim([y_lower y_upper])

        hold off

        % Export figure.

        filename = ['deterministic_N', num2str(N), ...
                    '_T', num2str(T), '.fig'];

        savefig(fig, ...
            fullfile(output_folder, filename));

        % exportgraphics(fig, ...
        %     fullfile(output_folder, filename), ...
        %     'Resolution', 300);

    end

end

%% PLOT ALL N,T PANELS IN ONE FIGURE

% Create one 9 x 9 inch figure containing all N,T combinations.

figure( ...
    'Units', 'inches', ...
    'Position', [1 1 9 9]);

tiledlayout( ...
    length(T_vals), ...
    length(N_vals), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:length(T_vals)

    for j = 1:length(N_vals)

        T = T_vals(i);
        N = N_vals(j);

        % Extract empirical sizes across all a values
        % for this particular (N,T) combination.

        cd_plot = squeeze(CD_size(i,j,:));
        cd_robust_plot = squeeze(CD_robust_size(i,j,:));
        cd_general_plot = squeeze(CD_general_size(i,j,:));

        nexttile

        % Plot standard CD empirical size.

        plot(a_vals, cd_plot, '-', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'CD');

        hold on

        % Plot robust CD empirical size.

        plot(a_vals, cd_robust_plot, '-', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Robust CD');

        % Plot general CD empirical size.

        plot(a_vals, cd_general_plot, '-', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'General CD');

        % Plot nominal significance level.

        yline(alpha, '--', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Nominal size');

        % Create one common legend using the first panel.

        if i == 1 && j == 1

            lgd = legend( ...
                'CD', 'Robust CD', 'General CD', 'Nominal size', ...
                'Orientation', 'horizontal');

            lgd.Layout.Tile = 'south';

        end

        % Axis labels.
        
        if i == length(T_vals), xlabel('a'), end
        if j == 1, ylabel('Empirical size'), end

        % Panel title.

        title(['N = ', num2str(N), ...
               ', T = ', num2str(T)])

        % Add grid.

        grid on

        % Keep same horizontal and vertical range for every panel.

        xlim([0 1])
        ylim([y_lower y_upper])

        hold off

    end

end

% Export figure.

savefig(gcf, ...
    fullfile(output_folder, 'deterministic_all_panels.fig'));

% exportgraphics(gcf, ...
%     fullfile(output_folder, 'deterministic_all_panels.png'), ...
%     'Resolution', 300);
