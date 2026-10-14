
%% ROLLING WINDOW CD ANALYSIS

clear
clc
close all

%% USER SETTINGS

% ABP data file.

abp_file = 'abp.csv';

% FRED U.S. real GDP data file.

gdpc1_file = 'gdpc1.csv';

% Significance level.

alpha = 0.05;

% Rolling window size.

window_length = 20;

% Overall analysis range.

analysis_start = "1975Q2";
analysis_end = "2014Q4";

%% CONVERT START AND END QUARTERS TO NUMERIC INDICES

% Parse starting quarter.

start_year = str2double(extractBefore(analysis_start, "Q"));
start_quarter = str2double(extractAfter(analysis_start, "Q"));

% Parse ending quarter.

end_year = str2double(extractBefore(analysis_end, "Q"));
end_quarter = str2double(extractAfter(analysis_end, "Q"));

% Convert YYYYQ# into a numeric index.

start_index = 4*start_year + start_quarter;
end_index = 4*end_year + end_quarter;

%% CONSTRUCT ROLLING WINDOWS

% Number of quarterly changes in the full analysis range.

number_changes = end_index - start_index + 1;

% Number of rolling windows.

number_windows = number_changes - window_length + 1;

% Store start and end quarter of every window.

rolling_start = strings(number_windows, 1);
rolling_end = strings(number_windows, 1);

for w = 1:number_windows

    % Numeric start and end index for this window.

    window_start_index = start_index + w - 1;
    window_end_index = window_start_index + window_length - 1;

    % Convert numeric start and end indices back to YYYYQ# labels.

    window_start_year = floor((window_start_index - 1)/4);
    window_start_quarter = mod(window_start_index - 1, 4) + 1;

    window_end_year = floor((window_end_index - 1)/4);
    window_end_quarter = mod(window_end_index - 1, 4) + 1;

    % Store readable labels.

    rolling_start(w) = ...
        string(window_start_year) + "Q" + string(window_start_quarter);

    rolling_end(w) = ...
        string(window_end_year) + "Q" + string(window_end_quarter);

end

%% PREALLOCATE STORAGE

% column 1 = Standard CD
% column 2 = Robust CD
% column 3 = General CD

deseasonalised_stats = NaN(number_windows, 3);
local_controls_stats = NaN(number_windows, 3);
national_local_controls_stats = NaN(number_windows, 3);
national_factors_stats = NaN(number_windows, 3);
national_regional_factors_stats = NaN(number_windows, 3);
national_regional_factors_local_controls_stats = NaN(number_windows, 3);
national_regional_factors_national_local_controls_stats = NaN(number_windows, 3);

% Store factor loadings from national + regional factors + local controls model.

rolling_factor_loading_tables = cell(number_windows, 1);

%% RUN ROLLING WINDOWS

for w = 1:number_windows

    % Current rolling window.

    start_quarter = rolling_start(w);
    end_quarter = rolling_end(w);

    fprintf( ...
        'Rolling window %d of %d: %s to %s\n', ...
        w, number_windows, start_quarter, end_quarter ...
        )

    % De-seasonalised real house-price changes only.

    evalc( ...
        'result_deseasonalised = main_deseasonalised(abp_file, start_quarter, end_quarter, alpha);' ...
        );

    deseasonalised_stats(w,:) = [ ...
        result_deseasonalised.cd_stat, ...
        result_deseasonalised.cd_robust_stat, ...
        result_deseasonalised.cd_general_stat ...
        ];

    % Local controls.

    evalc( ...
        'result_local_controls = main_local_controls(abp_file, start_quarter, end_quarter, alpha);' ...
        );

    local_controls_stats(w,:) = [ ...
        result_local_controls.cd_stat, ...
        result_local_controls.cd_robust_stat, ...
        result_local_controls.cd_general_stat ...
        ];

    % National + local controls.

    evalc( ...
        'result_controls_gdp = main_national_local_controls(abp_file, gdpc1_file, start_quarter, end_quarter, alpha);' ...
        );

    national_local_controls_stats(w,:) = [ ...
        result_controls_gdp.cd_stat, ...
        result_controls_gdp.cd_robust_stat, ...
        result_controls_gdp.cd_general_stat ...
        ];

    % National factor.

    evalc( ...
        'result_national = main_national_factors(abp_file, start_quarter, end_quarter, alpha);' ...
        );

    national_factors_stats(w,:) = [ ...
        result_national.cd_stat, ...
        result_national.cd_robust_stat, ...
        result_national.cd_general_stat ...
        ];

    % National + regional factors.

    evalc( ...
        'result_regional = main_national_regional_factors(abp_file, start_quarter, end_quarter, alpha);' ...
        );

    national_regional_factors_stats(w,:) = [ ...
        result_regional.cd_stat, ...
        result_regional.cd_robust_stat, ...
        result_regional.cd_general_stat ...
        ];

    % National + regional factors + local controls.

    evalc( ...
        'result_controls = main_national_regional_factors_local_controls(abp_file, start_quarter, end_quarter, alpha);' ...
        );

    national_regional_factors_local_controls_stats(w,:) = [ ...
        result_controls.cd_stat, ...
        result_controls.cd_robust_stat, ...
        result_controls.cd_general_stat ...
        ];

    N = result_controls.N;

    rolling_factor_loading_tables{w} = table( ...
        repmat(start_quarter, N, 1), ...
        repmat(end_quarter, N, 1), ...
        result_controls.msa_codes(:), ...
        result_controls.msa_region_codes(:), ...
        result_controls.pi_national_factor_loadings(:), ...
        result_controls.pi_regional_factor_loadings(:), ...
        result_controls.gpop_national_factor_loadings(:), ...
        result_controls.gpop_regional_factor_loadings(:), ...
        result_controls.ginc_national_factor_loadings(:), ...
        result_controls.ginc_regional_factor_loadings(:), ...
        'VariableNames', { ...
            'WindowStart', ...
            'WindowEnd', ...
            'MSACode', ...
            'RegionCode', ...
            'Pi_NationalLoading', ...
            'Pi_RegionalLoading', ...
            'GPOP_NationalLoading', ...
            'GPOP_RegionalLoading', ...
            'GINC_NationalLoading', ...
            'GINC_RegionalLoading' ...
            } ...
        );

    % National + regional factors + national + local controls.

    evalc( ...
        'result_full_controls = main_national_regional_factors_national_local_controls(abp_file, gdpc1_file, start_quarter, end_quarter, alpha);' ...
        );

    national_regional_factors_national_local_controls_stats(w,:) = [ ...
        result_full_controls.cd_stat, ...
        result_full_controls.cd_robust_stat, ...
        result_full_controls.cd_general_stat ...
        ];

end

%% STORE ROLLING FACTOR LOADINGS

% Stack the factor loadings from every window into one table.

rolling_factor_loadings = vertcat(rolling_factor_loading_tables{:});

%% STORE ROLLING RESULTS

rolling_results = table( ...
    rolling_start, ...
    rolling_end, ...
    deseasonalised_stats(:,1), ...
    deseasonalised_stats(:,2), ...
    deseasonalised_stats(:,3), ...
    local_controls_stats(:,1), ...
    local_controls_stats(:,2), ...
    local_controls_stats(:,3), ...
    national_local_controls_stats(:,1), ...
    national_local_controls_stats(:,2), ...
    national_local_controls_stats(:,3), ...
    national_factors_stats(:,1), ...
    national_factors_stats(:,2), ...
    national_factors_stats(:,3), ...
    national_regional_factors_stats(:,1), ...
    national_regional_factors_stats(:,2), ...
    national_regional_factors_stats(:,3), ...
    national_regional_factors_local_controls_stats(:,1), ...
    national_regional_factors_local_controls_stats(:,2), ...
    national_regional_factors_local_controls_stats(:,3), ...
    national_regional_factors_national_local_controls_stats(:,1), ...
    national_regional_factors_national_local_controls_stats(:,2), ...
    national_regional_factors_national_local_controls_stats(:,3), ...
    'VariableNames', { ...
        'WindowStart', ...
        'WindowEnd', ...
        'Deseasonalised_CD', ...
        'Deseasonalised_RobustCD', ...
        'Deseasonalised_GeneralCD', ...
        'LocalControls_CD', ...
        'LocalControls_RobustCD', ...
        'LocalControls_GeneralCD', ...
        'NationalLocalControls_CD', ...
        'NationalLocalControls_RobustCD', ...
        'NationalLocalControls_GeneralCD', ...
        'NationalFactors_CD', ...
        'NationalFactors_RobustCD', ...
        'NationalFactors_GeneralCD', ...
        'NationalRegionalFactors_CD', ...
        'NationalRegionalFactors_RobustCD', ...
        'NationalRegionalFactors_GeneralCD', ...
        'NationalRegionalFactorsLocalControls_CD', ...
        'NationalRegionalFactorsLocalControls_RobustCD', ...
        'NationalRegionalFactorsLocalControls_GeneralCD', ...
        'NationalRegionalFactorsNationalLocalControls_CD', ...
        'NationalRegionalFactorsNationalLocalControls_RobustCD', ...
        'NationalRegionalFactorsNationalLocalControls_GeneralCD' ...
        } ...
    );

%% STANDARD-NORMAL CRITICAL VALUE

critical_value = norminv(1 - alpha/2);

%% CONSTRUCT DATES FOR PLOTS

% Use the ending quarter of each rolling window as the date attached
% to the corresponding statistics.

plot_dates = NaT(number_windows, 1);

for w = 1:number_windows

    end_string = rolling_end(w);

    end_year = str2double(extractBefore(end_string, "Q"));
    end_quarter = str2double(extractAfter(end_string, "Q"));

    % Convert quarter number to ending month:
    %
    % Q1 -> March
    % Q2 -> June
    % Q3 -> September
    % Q4 -> December

    end_month = 3*end_quarter;

    plot_dates(w) = datetime(end_year, end_month, 1);

end

% Display one x-axis label every two years.

tick_positions = plot_dates(1):calyears(2):plot_dates(end);

%% COMMON Y-AXIS RANGE

% Combine all rolling statistics from all seven models.

all_stats = [ ...
    deseasonalised_stats(:); ...
    local_controls_stats(:); ...
    national_local_controls_stats(:); ...
    national_factors_stats(:); ...
    national_regional_factors_stats(:); ...
    national_regional_factors_local_controls_stats(:); ...
    national_regional_factors_national_local_controls_stats(:); ...
    critical_value; ...
    -critical_value ...
    ];

% Find overall minimum and maximum.

y_min = min(all_stats);
y_max = max(all_stats);

% Add a small amount of space above and below the plotted values.

y_range = y_max - y_min;

y_lower = y_min - 0.05*y_range;
y_upper = y_max + 0.05*y_range;

%% PLOT THE THREE CD STATISTICS

figure('Name', 'Rolling Window CD Statistics', 'Color', 'w')

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% FIGURE 1: DE-SEASONALISED REAL HOUSE-PRICE CHANGES ONLY

ax1 = nexttile;

plot(plot_dates, deseasonalised_stats(:,1), 'LineWidth', 1.5)
hold on
plot(plot_dates, deseasonalised_stats(:,2), 'LineWidth', 1.5)
plot(plot_dates, deseasonalised_stats(:,3), 'LineWidth', 1.5)

% Two-sided standard-normal critical values.

yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')

% Zero reference line.

yline(0, ':', 'HandleVisibility', 'off')
hold off

title('De-seasonalised Real House-Price Changes Only')

% xlabel('End quarter')
ylabel('Statistic')

% ylim([y_lower, y_upper])
xlim([plot_dates(1), plot_dates(end)])

xticks(tick_positions)
xticklabels([])

grid on

legend( ...
    'Standard CD', ...
    'Robust CD', ...
    'General CD', ...
    'Location', 'northwest' ...
    )

% %% FIGURE 2: LOCAL CONTROLS ONLY
% 
% ax2 = nexttile;
% 
% plot(plot_dates, local_controls_stats(:,1), 'LineWidth', 1.5)
% hold on
% plot(plot_dates, local_controls_stats(:,2), 'LineWidth', 1.5)
% plot(plot_dates, local_controls_stats(:,3), 'LineWidth', 1.5)
% 
% yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% 
% yline(0, ':', 'HandleVisibility', 'off')
% hold off
% 
% title('Local Controls')
% 
% % xlabel('End quarter')
% ylabel('Statistic')
% 
% % ylim([y_lower, y_upper])
% xlim([plot_dates(1), plot_dates(end)])
% 
% xticks(tick_positions)
% xticklabels([])
% 
% grid on

%% FIGURE 3: NATIONAL + LOCAL CONTROLS

ax3 = nexttile;

plot(plot_dates, national_local_controls_stats(:,1), 'LineWidth', 1.5)
hold on
plot(plot_dates, national_local_controls_stats(:,2), 'LineWidth', 1.5)
plot(plot_dates, national_local_controls_stats(:,3), 'LineWidth', 1.5)

yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')

yline(0, ':', 'HandleVisibility', 'off')
hold off

title('National + Local Controls')

% xlabel('End quarter')
ylabel('Statistic')

% ylim([y_lower, y_upper])
xlim([plot_dates(1), plot_dates(end)])

xticks(tick_positions)
xticklabels([])

grid on

% %% FIGURE 4: NATIONAL FACTOR
% 
% ax4 = nexttile;
% 
% plot(plot_dates, national_factors_stats(:,1), 'LineWidth', 1.5)
% hold on
% plot(plot_dates, national_factors_stats(:,2), 'LineWidth', 1.5)
% plot(plot_dates, national_factors_stats(:,3), 'LineWidth', 1.5)
% 
% % Two-sided standard-normal critical values.
% 
% yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% 
% % Zero reference line.
% 
% yline(0, ':', 'HandleVisibility', 'off')
% hold off
% 
% title('National Factor')
% 
% % xlabel('End quarter')
% ylabel('Statistic')
% 
% % ylim([y_lower, y_upper])
% xlim([plot_dates(1), plot_dates(end)])
% 
% xticks(tick_positions)
% xticklabels([])
% 
% grid on

% %% FIGURE 5: NATIONAL + REGIONAL FACTORS
% 
% ax5 = nexttile;
% 
% plot(plot_dates, national_regional_factors_stats(:,1), 'LineWidth', 1.5)
% hold on
% plot(plot_dates, national_regional_factors_stats(:,2), 'LineWidth', 1.5)
% plot(plot_dates, national_regional_factors_stats(:,3), 'LineWidth', 1.5)
% 
% yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% 
% yline(0, ':', 'HandleVisibility', 'off')
% hold off
% 
% title('National + Regional Factors')
% 
% % xlabel('End quarter')
% ylabel('Statistic')
% 
% % ylim([y_lower, y_upper])
% xlim([plot_dates(1), plot_dates(end)])
% 
% xticks(tick_positions)
% xticklabels([])
% 
% grid on

%% FIGURE 6: NATIONAL + REGIONAL FACTORS + LOCAL CONTROLS

ax6 = nexttile;

plot(plot_dates, national_regional_factors_local_controls_stats(:,1), 'LineWidth', 1.5)
hold on
plot(plot_dates, national_regional_factors_local_controls_stats(:,2), 'LineWidth', 1.5)
plot(plot_dates, national_regional_factors_local_controls_stats(:,3), 'LineWidth', 1.5)

yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')

yline(0, ':', 'HandleVisibility', 'off')
hold off

title('National + Regional Factors + Local Controls')

xlabel('End quarter')
ylabel('Statistic')

% ylim([y_lower, y_upper])
xlim([plot_dates(1), plot_dates(end)])

xticks(tick_positions)
xtickformat('yyyy')
% xticklabels([])

grid on

% %% FIGURE 7: NATIONAL + REGIONAL FACTORS + NATIONAL + LOCAL CONTROLS
% 
% ax7 = nexttile;
% 
% plot(plot_dates, national_regional_factors_national_local_controls_stats(:,1), 'LineWidth', 1.5)
% hold on
% plot(plot_dates, national_regional_factors_national_local_controls_stats(:,2), 'LineWidth', 1.5)
% plot(plot_dates, national_regional_factors_national_local_controls_stats(:,3), 'LineWidth', 1.5)
% 
% yline(critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% yline(-critical_value, '--', 'LineWidth', 1, 'HandleVisibility', 'off')
% 
% yline(0, ':', 'HandleVisibility', 'off')
% hold off
% 
% title('National + Regional Factors + National + Local Controls')
% 
% xlabel('End quarter')
% ylabel('Statistic')
% 
% % ylim([y_lower, y_upper])
% xlim([plot_dates(1), plot_dates(end)])
% 
% xticks(tick_positions)
% xtickformat('yyyy')
% % xtickangle(45)
% 
% grid on

%% LINK AXES

% Ensure all seven plots keep the same horizontal scale.

% linkaxes([ax1, ax2, ax3, ax4, ax5, ax6, ax7], 'xy')
linkaxes([ax1, ax3, ax6], 'x')

%% OVERALL TITLE

sgtitle('Rolling Window Cross-Sectional Dependence Statistics')

%% SAVE RESULTS

% Store exported figures in their own output subfolder.

output_folder = fullfile('output', 'rolling');

% Create the folder if it does not already exist.

if ~exist(output_folder, 'dir')

    mkdir(output_folder);

end

% Export results.

writetable( ...
    rolling_factor_loadings, ...
    fullfile(output_folder, 'rolling_factor_loadings_5y.csv') ...
    );

% Export factor loadings.

writetable( ...
    rolling_results, ...
    fullfile(output_folder, 'rolling_results_three_models_5y.csv') ...
    );

% Export figure.

exportgraphics( ...
    gcf, ...
    fullfile(output_folder, 'rolling_figure_three_models_5y.png'), ...
    'Resolution', 300 ...
    );