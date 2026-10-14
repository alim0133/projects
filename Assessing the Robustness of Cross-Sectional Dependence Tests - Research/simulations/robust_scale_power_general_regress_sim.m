
clear all
clc

% Set seed.

rng(36006963);

%% MONTE CARLO SETTINGS

% Time dimensions.

T_vals = [10 20 50 100];

% Cross-sectional dimensions.

N_vals = [10 20 50 100];

% Number of Monte Carlo replications.

reps = 5000;

% Significance level.

alpha = 0.05;

% Factor loading ranges.
%
% Weak factor:
% gamma_i ~ U(0.1,0.3)
%
% Strong factor:
% gamma_i ~ U(0.3,0.5)

gamma_weak = [0.1 0.3];

gamma_strong = [0.3 0.5];

%% PREALLOCATE SIZE STORAGE

% Standard CD.

CD_size_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_size_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_size_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% Robust CD.

CD_robust_size_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_size_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_size_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% General CD.

CD_general_size_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_general_size_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_general_size_stochastic = ...
    zeros(length(T_vals), length(N_vals));

%% PREALLOCATE CD MEAN STORAGE

% Standard CD.

CD_mean_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_mean_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_mean_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% Robust CD.

CD_robust_mean_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_mean_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_mean_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% General CD.

CD_general_mean_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_general_mean_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_general_mean_stochastic = ...
    zeros(length(T_vals), length(N_vals));

%% PREALLOCATE CD VARIANCE STORAGE

% Standard CD.

CD_var_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_var_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_var_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% Robust CD.

CD_robust_var_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_var_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_var_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% General CD.

CD_general_var_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_general_var_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_general_var_stochastic = ...
    zeros(length(T_vals), length(N_vals));

%% PREALLOCATE PAIRWISE DIAGNOSTIC STORAGE

% Each cell corresponds to one (T,N) combination.
%
% Within each cell:
% rows/entries = individual pairs (i,j)
%
% For example, when N = 10 there are:
% N(N-1)/2 = 45 pairwise statistics.

%% UNIT SCALE

% Pearson pairwise statistic:
% sqrt(T) * rho_hat_ij

pearson_pair_mean_iid = ...
    cell(length(T_vals), length(N_vals));

pearson_pair_var_iid = ...
    cell(length(T_vals), length(N_vals));

% Giraitis pairwise statistic:
% t_tilde_ij

giraitis_pair_mean_iid = ...
    cell(length(T_vals), length(N_vals));

giraitis_pair_var_iid = ...
    cell(length(T_vals), length(N_vals));

%% DETERMINISTIC SCALE

pearson_pair_mean_deterministic = ...
    cell(length(T_vals), length(N_vals));

pearson_pair_var_deterministic = ...
    cell(length(T_vals), length(N_vals));

giraitis_pair_mean_deterministic = ...
    cell(length(T_vals), length(N_vals));

giraitis_pair_var_deterministic = ...
    cell(length(T_vals), length(N_vals));

%% STOCHASTIC SCALE

pearson_pair_mean_stochastic = ...
    cell(length(T_vals), length(N_vals));

pearson_pair_var_stochastic = ...
    cell(length(T_vals), length(N_vals));

giraitis_pair_mean_stochastic = ...
    cell(length(T_vals), length(N_vals));

giraitis_pair_var_stochastic = ...
    cell(length(T_vals), length(N_vals));

%% SIZE EXPERIMENT

for i = 1:length(T_vals)

    for j = 1:length(N_vals)

        T = T_vals(i);
        N = N_vals(j);

        %% UNIT SCALE

        % Panel regression:
        % y_it = 0.2 + 0.5 x_it + u_it
        % u_it = h_t epsilon_it
        %
        % where:
        % h_t = 1

        result = ...
            robust_scale_power_general_regress_mc(N, T, reps, 'iid', 'size', [0 0], alpha);

        % Store empirical sizes.

        CD_size_iid(i,j) = ...
            result.cd_rejection_rate;

        CD_robust_size_iid(i,j) = ...
            result.cd_robust_rejection_rate;

        CD_general_size_iid(i,j) = ...
            result.cd_general_rejection_rate;

        % Store mean of CD statistics.

        CD_mean_iid(i,j) = ...
            result.mean_cd;

        CD_robust_mean_iid(i,j) = ...
            result.mean_cd_robust;

        CD_general_mean_iid(i,j) = ...
            result.mean_cd_general;

        % Store variance of CD statistics.

        CD_var_iid(i,j) = ...
            result.var_cd;

        CD_robust_var_iid(i,j) = ...
            result.var_cd_robust;

        CD_general_var_iid(i,j) = ...
            result.var_cd_general;

        % Pairwise diagnostics

        % result.pearson_stats is a reps x 1 cell array.
        %
        % Each cell contains all pairwise Pearson statistics from one
        % Monte Carlo replication.
        %
        % Convert this into a matrix:
        % rows    = Monte Carlo replications
        % columns = pairs (i,j)

        pearson_matrix = ...
            cell2mat(result.pearson_stats');

        giraitis_matrix = ...
            cell2mat(result.giraitis_stats');

        % Monte Carlo mean of each pair.

        pearson_pair_mean_iid{i,j} = ...
            mean(pearson_matrix, 1);

        giraitis_pair_mean_iid{i,j} = ...
            mean(giraitis_matrix, 1);

        % Monte Carlo variance of each pair.

        pearson_pair_var_iid{i,j} = ...
            var(pearson_matrix, 0, 1);

        giraitis_pair_var_iid{i,j} = ...
            var(giraitis_matrix, 0, 1);

        %% DETERMINISTIC SCALE

        % Panel regression:
        % y_it = 0.2 + 0.5 x_it + u_it
        % u_it = h_t epsilon_it
        %
        % where:
        % h_t = (3/T) floor(t/10)

        result = ...
            robust_scale_power_general_regress_mc(N, T, reps, 'deterministic', 'size', [0 0], alpha);

        % Store empirical sizes.

        CD_size_deterministic(i,j) = ...
            result.cd_rejection_rate;

        CD_robust_size_deterministic(i,j) = ...
            result.cd_robust_rejection_rate;

        CD_general_size_deterministic(i,j) = ...
            result.cd_general_rejection_rate;

        % Store mean of CD statistics.

        CD_mean_deterministic(i,j) = ...
            result.mean_cd;

        CD_robust_mean_deterministic(i,j) = ...
            result.mean_cd_robust;

        CD_general_mean_deterministic(i,j) = ...
            result.mean_cd_general;

        % Store variance of CD statistics.

        CD_var_deterministic(i,j) = ...
            result.var_cd;

        CD_robust_var_deterministic(i,j) = ...
            result.var_cd_robust;

        CD_general_var_deterministic(i,j) = ...
            result.var_cd_general;

        % Pairwise diagnostics

        pearson_matrix = ...
            cell2mat(result.pearson_stats');

        giraitis_matrix = ...
            cell2mat(result.giraitis_stats');

        % Monte Carlo mean of each pair.

        pearson_pair_mean_deterministic{i,j} = ...
            mean(pearson_matrix, 1);

        giraitis_pair_mean_deterministic{i,j} = ...
            mean(giraitis_matrix, 1);

        % Monte Carlo variance of each pair.

        pearson_pair_var_deterministic{i,j} = ...
            var(pearson_matrix, 0, 1);

        giraitis_pair_var_deterministic{i,j} = ...
            var(giraitis_matrix, 0, 1);

        %% STOCHASTIC SCALE

        % Panel regression:
        % y_it = 0.2 + 0.5 x_it + u_it
        % u_it = h_t epsilon_it
        %
        % where:
        % h_t = |T^(-1/2) sum_{s=1}^{t} zeta_s|

        result = ...
            robust_scale_power_general_regress_mc(N, T, reps, 'stochastic', 'size', [0 0], alpha);

        % Store empirical sizes.

        CD_size_stochastic(i,j) = ...
            result.cd_rejection_rate;

        CD_robust_size_stochastic(i,j) = ...
            result.cd_robust_rejection_rate;

        CD_general_size_stochastic(i,j) = ...
            result.cd_general_rejection_rate;

        % Store mean of CD statistics.

        CD_mean_stochastic(i,j) = ...
            result.mean_cd;

        CD_robust_mean_stochastic(i,j) = ...
            result.mean_cd_robust;

        CD_general_mean_stochastic(i,j) = ...
            result.mean_cd_general;

        % Store variance of CD statistics.

        CD_var_stochastic(i,j) = ...
            result.var_cd;

        CD_robust_var_stochastic(i,j) = ...
            result.var_cd_robust;

        CD_general_var_stochastic(i,j) = ...
            result.var_cd_general;

        % Pairwise diagnostics

        pearson_matrix = ...
            cell2mat(result.pearson_stats');

        giraitis_matrix = ...
            cell2mat(result.giraitis_stats');

        % Monte Carlo mean of each pair.

        pearson_pair_mean_stochastic{i,j} = ...
            mean(pearson_matrix, 1);

        giraitis_pair_mean_stochastic{i,j} = ...
            mean(giraitis_matrix, 1);

        % Monte Carlo variance of each pair.

        pearson_pair_var_stochastic{i,j} = ...
            var(pearson_matrix, 0, 1);

        giraitis_pair_var_stochastic{i,j} = ...
            var(giraitis_matrix, 0, 1);

    end

end

%% PREALLOCATE POWER STORAGE: WEAK FACTOR

% Standard CD.

CD_power_weak_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_power_weak_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_power_weak_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% Robust CD.

CD_robust_power_weak_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_power_weak_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_robust_power_weak_stochastic = ...
    zeros(length(T_vals), length(N_vals));

% General CD.

CD_general_power_weak_iid = ...
    zeros(length(T_vals), length(N_vals));

CD_general_power_weak_deterministic = ...
    zeros(length(T_vals), length(N_vals));

CD_general_power_weak_stochastic = ...
    zeros(length(T_vals), length(N_vals));

%% POWER EXPERIMENT: gamma_i ~ U(0.1,0.3)

for i = 1:length(T_vals)

    for j = 1:length(N_vals)

        T = T_vals(i);
        N = N_vals(j);

        %% UNIT SCALE

        result = ...
            robust_scale_power_general_regress_mc(N, T, reps, 'iid', 'power', gamma_weak, alpha);

        % Store empirical powers.

        CD_power_weak_iid(i,j) = ...
            result.cd_rejection_rate;

        CD_robust_power_weak_iid(i,j) = ...
            result.cd_robust_rejection_rate;

        CD_general_power_weak_iid(i,j) = ...
            result.cd_general_rejection_rate;

        %% DETERMINISTIC SCALE

        result = ...
            robust_scale_power_general_regress_mc(N, T, reps, 'deterministic', 'power', gamma_weak, alpha);

        % Store empirical powers.

        CD_power_weak_deterministic(i,j) = ...
            result.cd_rejection_rate;

        CD_robust_power_weak_deterministic(i,j) = ...
            result.cd_robust_rejection_rate;

        CD_general_power_weak_deterministic(i,j) = ...
            result.cd_general_rejection_rate;

        %% STOCHASTIC SCALE

        result = ...
            robust_scale_power_general_regress_mc(N, T, reps, 'stochastic', 'power', gamma_weak, alpha);

        % Store empirical powers.

        CD_power_weak_stochastic(i,j) = ...
            result.cd_rejection_rate;

        CD_robust_power_weak_stochastic(i,j) = ...
            result.cd_robust_rejection_rate;

        CD_general_power_weak_stochastic(i,j) = ...
            result.cd_general_rejection_rate;

    end

end

% %% PREALLOCATE POWER STORAGE: STRONG FACTOR
% 
% % Standard CD.
% 
% CD_power_strong_iid = ...
%     zeros(length(T_vals), length(N_vals));
% 
% CD_power_strong_deterministic = ...
%     zeros(length(T_vals), length(N_vals));
% 
% CD_power_strong_stochastic = ...
%     zeros(length(T_vals), length(N_vals));
% 
% % Robust CD.
% 
% CD_robust_power_strong_iid = ...
%     zeros(length(T_vals), length(N_vals));
% 
% CD_robust_power_strong_deterministic = ...
%     zeros(length(T_vals), length(N_vals));
% 
% CD_robust_power_strong_stochastic = ...
%     zeros(length(T_vals), length(N_vals));
% 
% % General CD.
% 
% CD_general_power_strong_iid = ...
%     zeros(length(T_vals), length(N_vals));
% 
% CD_general_power_strong_deterministic = ...
%     zeros(length(T_vals), length(N_vals));
% 
% CD_general_power_strong_stochastic = ...
%     zeros(length(T_vals), length(N_vals));
% 
% %% POWER EXPERIMENT: gamma_i ~ U(0.3,0.5)
% 
% for i = 1:length(T_vals)
% 
%     for j = 1:length(N_vals)
% 
%         T = T_vals(i);
%         N = N_vals(j);
% 
%         %% UNIT SCALE
% 
%         result = ...
%             robust_scale_power_general_regress_mc(N, T, reps, 'iid', 'power', gamma_strong, alpha);
% 
%         % Store empirical powers.
% 
%         CD_power_strong_iid(i,j) = ...
%             result.cd_rejection_rate;
% 
%         CD_robust_power_strong_iid(i,j) = ...
%             result.cd_robust_rejection_rate;
% 
%         CD_general_power_strong_iid(i,j) = ...
%             result.cd_general_rejection_rate;
% 
%         %% DETERMINISTIC SCALE
% 
%         result = ...
%             robust_scale_power_general_regress_mc(N, T, reps, 'deterministic', 'power', gamma_strong, alpha);
% 
%         % Store empirical powers.
% 
%         CD_power_strong_deterministic(i,j) = ...
%             result.cd_rejection_rate;
% 
%         CD_robust_power_strong_deterministic(i,j) = ...
%             result.cd_robust_rejection_rate;
% 
%         CD_general_power_strong_deterministic(i,j) = ...
%             result.cd_general_rejection_rate;
% 
%         %% STOCHASTIC SCALE
% 
%         result = ...
%             robust_scale_power_general_regress_mc(N, T, reps, 'stochastic', 'power', gamma_strong, alpha);
% 
%         % Store empirical powers.
% 
%         CD_power_strong_stochastic(i,j) = ...
%             result.cd_rejection_rate;
% 
%         CD_robust_power_strong_stochastic(i,j) = ...
%             result.cd_robust_rejection_rate;
% 
%         CD_general_power_strong_stochastic(i,j) = ...
%             result.cd_general_rejection_rate;
% 
%     end
% 
% end

%% EMPIRICAL SIZE TABLES

CD_size_iid_table = array2table( ...
    CD_size_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_size_iid_table = array2table( ...
    CD_robust_size_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_size_iid_table = array2table( ...
    CD_general_size_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_size_deterministic_table = array2table( ...
    CD_size_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_size_deterministic_table = array2table( ...
    CD_robust_size_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_size_deterministic_table = array2table( ...
    CD_general_size_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_size_stochastic_table = array2table( ...
    CD_size_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_size_stochastic_table = array2table( ...
    CD_robust_size_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_size_stochastic_table = array2table( ...
    CD_general_size_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

%% EMPIRICAL POWER TABLES: WEAK FACTOR

CD_power_weak_iid_table = array2table( ...
    CD_power_weak_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_power_weak_iid_table = array2table( ...
    CD_robust_power_weak_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_power_weak_iid_table = array2table( ...
    CD_general_power_weak_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_power_weak_deterministic_table = array2table( ...
    CD_power_weak_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_power_weak_deterministic_table = array2table( ...
    CD_robust_power_weak_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_power_weak_deterministic_table = array2table( ...
    CD_general_power_weak_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_power_weak_stochastic_table = array2table( ...
    CD_power_weak_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_power_weak_stochastic_table = array2table( ...
    CD_robust_power_weak_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_power_weak_stochastic_table = array2table( ...
    CD_general_power_weak_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

% %% EMPIRICAL POWER TABLES: STRONG FACTOR
% 
% CD_power_strong_iid_table = array2table( ...
%     CD_power_strong_iid, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_robust_power_strong_iid_table = array2table( ...
%     CD_robust_power_strong_iid, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_general_power_strong_iid_table = array2table( ...
%     CD_general_power_strong_iid, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_power_strong_deterministic_table = array2table( ...
%     CD_power_strong_deterministic, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_robust_power_strong_deterministic_table = array2table( ...
%     CD_robust_power_strong_deterministic, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_general_power_strong_deterministic_table = array2table( ...
%     CD_general_power_strong_deterministic, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_power_strong_stochastic_table = array2table( ...
%     CD_power_strong_stochastic, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_robust_power_strong_stochastic_table = array2table( ...
%     CD_robust_power_strong_stochastic, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));
% 
% CD_general_power_strong_stochastic_table = array2table( ...
%     CD_general_power_strong_stochastic, ...
%     'VariableNames', strcat("N", string(N_vals)), ...
%     'RowNames', strcat("T", string(T_vals)));

%% DISPLAY SIZE RESULTS

disp(' ')
disp('STANDARD CD SIZE: I.I.D.')
disp(CD_size_iid_table)

disp(' ')
disp('ROBUST CD SIZE: I.I.D.')
disp(CD_robust_size_iid_table)

disp(' ')
disp('GENERAL CD SIZE: I.I.D.')
disp(CD_general_size_iid_table)


disp(' ')
disp('STANDARD CD SIZE: DETERMINISTIC')
disp(CD_size_deterministic_table)

disp(' ')
disp('ROBUST CD SIZE: DETERMINISTIC')
disp(CD_robust_size_deterministic_table)

disp(' ')
disp('GENERAL CD SIZE: DETERMINISTIC')
disp(CD_general_size_deterministic_table)


disp(' ')
disp('STANDARD CD SIZE: STOCHASTIC')
disp(CD_size_stochastic_table)

disp(' ')
disp('ROBUST CD SIZE: STOCHASTIC')
disp(CD_robust_size_stochastic_table)

disp(' ')
disp('GENERAL CD SIZE: STOCHASTIC')
disp(CD_general_size_stochastic_table)

%% DISPLAY POWER RESULTS: WEAK FACTOR

disp(' ')
disp('STANDARD CD POWER: U(0.1,0.3), I.I.D.')
disp(CD_power_weak_iid_table)

disp(' ')
disp('ROBUST CD POWER: U(0.1,0.3), I.I.D.')
disp(CD_robust_power_weak_iid_table)

disp(' ')
disp('GENERAL CD POWER: U(0.1,0.3), I.I.D.')
disp(CD_general_power_weak_iid_table)


disp(' ')
disp('STANDARD CD POWER: U(0.1,0.3), DETERMINISTIC')
disp(CD_power_weak_deterministic_table)

disp(' ')
disp('ROBUST CD POWER: U(0.1,0.3), DETERMINISTIC')
disp(CD_robust_power_weak_deterministic_table)

disp(' ')
disp('GENERAL CD POWER: U(0.1,0.3), DETERMINISTIC')
disp(CD_general_power_weak_deterministic_table)


disp(' ')
disp('STANDARD CD POWER: U(0.1,0.3), STOCHASTIC')
disp(CD_power_weak_stochastic_table)

disp(' ')
disp('ROBUST CD POWER: U(0.1,0.3), STOCHASTIC')
disp(CD_robust_power_weak_stochastic_table)

disp(' ')
disp('GENERAL CD POWER: U(0.1,0.3), STOCHASTIC')
disp(CD_general_power_weak_stochastic_table)

% %% DISPLAY POWER RESULTS: STRONG FACTOR
% 
% disp(' ')
% disp('STANDARD CD POWER: U(0.3,0.5), I.I.D.')
% disp(CD_power_strong_iid_table)
% 
% disp(' ')
% disp('ROBUST CD POWER: U(0.3,0.5), I.I.D.')
% disp(CD_robust_power_strong_iid_table)
% 
% disp(' ')
% disp('GENERAL CD POWER: U(0.3,0.5), I.I.D.')
% disp(CD_general_power_strong_iid_table)
% 
% 
% disp(' ')
% disp('STANDARD CD POWER: U(0.3,0.5), DETERMINISTIC')
% disp(CD_power_strong_deterministic_table)
% 
% disp(' ')
% disp('ROBUST CD POWER: U(0.3,0.5), DETERMINISTIC')
% disp(CD_robust_power_strong_deterministic_table)
% 
% disp(' ')
% disp('GENERAL CD POWER: U(0.3,0.5), DETERMINISTIC')
% disp(CD_general_power_strong_deterministic_table)
% 
% 
% disp(' ')
% disp('STANDARD CD POWER: U(0.3,0.5), STOCHASTIC')
% disp(CD_power_strong_stochastic_table)
% 
% disp(' ')
% disp('ROBUST CD POWER: U(0.3,0.5), STOCHASTIC')
% disp(CD_robust_power_strong_stochastic_table)
% 
% disp(' ')
% disp('GENERAL CD POWER: U(0.3,0.5), STOCHASTIC')
% disp(CD_general_power_strong_stochastic_table)

%% CD MEAN TABLES

% Standard CD.

CD_mean_iid_table = array2table( ...
    CD_mean_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_mean_deterministic_table = array2table( ...
    CD_mean_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_mean_stochastic_table = array2table( ...
    CD_mean_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

% Robust CD.

CD_robust_mean_iid_table = array2table( ...
    CD_robust_mean_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_mean_deterministic_table = array2table( ...
    CD_robust_mean_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_mean_stochastic_table = array2table( ...
    CD_robust_mean_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

% General CD.

CD_general_mean_iid_table = array2table( ...
    CD_general_mean_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_mean_deterministic_table = array2table( ...
    CD_general_mean_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_mean_stochastic_table = array2table( ...
    CD_general_mean_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

%% CD VARIANCE TABLES

% Standard CD.

CD_var_iid_table = array2table( ...
    CD_var_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_var_deterministic_table = array2table( ...
    CD_var_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_var_stochastic_table = array2table( ...
    CD_var_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

% Robust CD.

CD_robust_var_iid_table = array2table( ...
    CD_robust_var_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_var_deterministic_table = array2table( ...
    CD_robust_var_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_robust_var_stochastic_table = array2table( ...
    CD_robust_var_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

% General CD.

CD_general_var_iid_table = array2table( ...
    CD_general_var_iid, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_var_deterministic_table = array2table( ...
    CD_general_var_deterministic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

CD_general_var_stochastic_table = array2table( ...
    CD_general_var_stochastic, ...
    'VariableNames', strcat("N", string(N_vals)), ...
    'RowNames', strcat("T", string(T_vals)));

%% DISPLAY CD MEANS

disp(' ')
disp('MEAN OF STANDARD CD: UNIT SCALE')
disp(CD_mean_iid_table)

disp(' ')
disp('MEAN OF ROBUST CD: UNIT SCALE')
disp(CD_robust_mean_iid_table)

disp(' ')
disp('MEAN OF GENERAL CD: UNIT SCALE')
disp(CD_general_mean_iid_table)


disp(' ')
disp('MEAN OF STANDARD CD: DETERMINISTIC SCALE')
disp(CD_mean_deterministic_table)

disp(' ')
disp('MEAN OF ROBUST CD: DETERMINISTIC SCALE')
disp(CD_robust_mean_deterministic_table)

disp(' ')
disp('MEAN OF GENERAL CD: DETERMINISTIC SCALE')
disp(CD_general_mean_deterministic_table)


disp(' ')
disp('MEAN OF STANDARD CD: STOCHASTIC SCALE')
disp(CD_mean_stochastic_table)

disp(' ')
disp('MEAN OF ROBUST CD: STOCHASTIC SCALE')
disp(CD_robust_mean_stochastic_table)

disp(' ')
disp('MEAN OF GENERAL CD: STOCHASTIC SCALE')
disp(CD_general_mean_stochastic_table)

%% DISPLAY CD VARIANCES

disp(' ')
disp('VARIANCE OF STANDARD CD: UNIT SCALE')
disp(CD_var_iid_table)

disp(' ')
disp('VARIANCE OF ROBUST CD: UNIT SCALE')
disp(CD_robust_var_iid_table)

disp(' ')
disp('VARIANCE OF GENERAL CD: UNIT SCALE')
disp(CD_general_var_iid_table)


disp(' ')
disp('VARIANCE OF STANDARD CD: DETERMINISTIC SCALE')
disp(CD_var_deterministic_table)

disp(' ')
disp('VARIANCE OF ROBUST CD: DETERMINISTIC SCALE')
disp(CD_robust_var_deterministic_table)

disp(' ')
disp('VARIANCE OF GENERAL CD: DETERMINISTIC SCALE')
disp(CD_general_var_deterministic_table)


disp(' ')
disp('VARIANCE OF STANDARD CD: STOCHASTIC SCALE')
disp(CD_var_stochastic_table)

disp(' ')
disp('VARIANCE OF ROBUST CD: STOCHASTIC SCALE')
disp(CD_robust_var_stochastic_table)

disp(' ')
disp('VARIANCE OF GENERAL CD: STOCHASTIC SCALE')
disp(CD_general_var_stochastic_table)

%% INSPECT PAIRWISE DIAGNOSTICS

% Example:
% T = 100 and N = 50.
%
% T = 100 is position 4 in T_vals.
% N = 50 is position 3 in N_vals.

% i_example = 4;
% j_example = 3;

% %% UNIT SCALE
% 
% disp(' ')
% disp('PAIRWISE PEARSON MEANS: I.I.D., T=100, N=50')
% disp(pearson_pair_mean_iid{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE PEARSON VARIANCES: I.I.D., T=100, N=50')
% disp(pearson_pair_var_iid{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE GIRAITIS MEANS: I.I.D., T=100, N=50')
% disp(giraitis_pair_mean_iid{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE GIRAITIS VARIANCES: I.I.D., T=100, N=50')
% disp(giraitis_pair_var_iid{i_example,j_example}')

%% DETERMINISTIC SCALE

% disp(' ')
% disp('PAIRWISE PEARSON MEANS: DETERMINISTIC, T=100, N=50')
% disp(pearson_pair_mean_deterministic{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE PEARSON VARIANCES: DETERMINISTIC, T=100, N=50')
% disp(pearson_pair_var_deterministic{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE GIRAITIS MEANS: DETERMINISTIC, T=100, N=50')
% disp(giraitis_pair_mean_deterministic{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE GIRAITIS VARIANCES: DETERMINISTIC, T=100, N=50')
% disp(giraitis_pair_var_deterministic{i_example,j_example}')

% %% STOCHASTIC SCALE
% 
% disp(' ')
% disp('PAIRWISE PEARSON MEANS: STOCHASTIC, T=100, N=50')
% disp(pearson_pair_mean_stochastic{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE PEARSON VARIANCES: STOCHASTIC, T=100, N=50')
% disp(pearson_pair_var_stochastic{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE GIRAITIS MEANS: STOCHASTIC, T=100, N=50')
% disp(giraitis_pair_mean_stochastic{i_example,j_example}')
% 
% disp(' ')
% disp('PAIRWISE GIRAITIS VARIANCES: STOCHASTIC, T=100, N=50')
% disp(giraitis_pair_var_stochastic{i_example,j_example}')