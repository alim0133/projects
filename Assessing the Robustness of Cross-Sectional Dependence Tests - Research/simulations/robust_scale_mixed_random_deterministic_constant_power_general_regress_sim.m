
clear all
clc

% Set seed.

rng(36006963);

%% MONTE CARLO SETTINGS

% Time dimensions.

T_vals = [10 20 50 100];

% Cross-sectional dimensions.

N_vals = [10 20 50 100];

% Pervasiveness parameter 0 <= a <= 1.
%
% Number of deterministic-scale individuals:
% N_deterministic = floor(N^a)

a_vals = [0.6 0.75 0.9];

% Number of Monte Carlo replications.

reps = 5000;

% Significance level.

alpha = 0.05;

% Factor loading range under power:
% gamma_i ~ U(0.1,0.3)

gamma_weak = [0.1 0.3];

%% PREALLOCATE SIZE STORAGE

% Dimensions: T x N x a

CD_size = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

CD_robust_size = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

CD_general_size = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

%% PREALLOCATE POWER STORAGE

CD_power = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

CD_robust_power = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

CD_general_power = ...
    zeros(length(T_vals), length(N_vals), length(a_vals));

%% PREALLOCATE NUMBER OF DETERMINISTIC UNITS

N_deterministic_storage = ...
    zeros(length(N_vals), length(a_vals));

%% SIZE EXPERIMENT

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

%% POWER EXPERIMENT

for k = 1:length(a_vals)

    a = a_vals(k);

    for i = 1:length(T_vals)

        for j = 1:length(N_vals)

            T = T_vals(i);
            N = N_vals(j);

            % Power DGP (gamma_i ~ U(0.1,0.3)):
            % y_it = beta_0 + beta_1 x_it + u_it
            % u_it = gamma_i f_t + h_it epsilon_it
            %
            % Deterministic-scale units:
            % h_it = (3/T) floor(t/10)
            %
            % Remaining units:
            % h_it = 1

            result = ...
                robust_scale_mixed_random_deterministic_constant_power_general_regress_mc( ...
                    N, T, reps, a, 'power', gamma_weak, alpha);

            % Store empirical powers.

            CD_power(i,j,k) = ...
                result.cd_rejection_rate;

            CD_robust_power(i,j,k) = ...
                result.cd_robust_rejection_rate;

            CD_general_power(i,j,k) = ...
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

%% DISPLAY SIZE RESULTS

for k = 1:length(a_vals)

    a = a_vals(k);

    CD_size_table = array2table( ...
        CD_size(:,:,k), ...
        'VariableNames', strcat("N", string(N_vals)), ...
        'RowNames', strcat("T", string(T_vals)));

    CD_robust_size_table = array2table( ...
        CD_robust_size(:,:,k), ...
        'VariableNames', strcat("N", string(N_vals)), ...
        'RowNames', strcat("T", string(T_vals)));

    CD_general_size_table = array2table( ...
        CD_general_size(:,:,k), ...
        'VariableNames', strcat("N", string(N_vals)), ...
        'RowNames', strcat("T", string(T_vals)));

    disp(' ')
    disp(['STANDARD CD EMPIRICAL SIZE: a = ', num2str(a)])
    disp(CD_size_table)

    disp(' ')
    disp(['ROBUST CD EMPIRICAL SIZE: a = ', num2str(a)])
    disp(CD_robust_size_table)

    disp(' ')
    disp(['GENERAL CD EMPIRICAL SIZE: a = ', num2str(a)])
    disp(CD_general_size_table)

end

%% DISPLAY POWER RESULTS

for k = 1:length(a_vals)

    a = a_vals(k);

    CD_power_table = array2table( ...
        CD_power(:,:,k), ...
        'VariableNames', strcat("N", string(N_vals)), ...
        'RowNames', strcat("T", string(T_vals)));

    CD_robust_power_table = array2table( ...
        CD_robust_power(:,:,k), ...
        'VariableNames', strcat("N", string(N_vals)), ...
        'RowNames', strcat("T", string(T_vals)));

    CD_general_power_table = array2table( ...
        CD_general_power(:,:,k), ...
        'VariableNames', strcat("N", string(N_vals)), ...
        'RowNames', strcat("T", string(T_vals)));

    disp(' ')
    disp(['STANDARD CD POWER: a = ', num2str(a), ...
          ', gamma_i ~ U(0.1,0.3)'])
    disp(CD_power_table)

    disp(' ')
    disp(['ROBUST CD POWER: a = ', num2str(a), ...
          ', gamma_i ~ U(0.1,0.3)'])
    disp(CD_robust_power_table)

    disp(' ')
    disp(['GENERAL CD POWER: a = ', num2str(a), ...
          ', gamma_i ~ U(0.1,0.3)'])
    disp(CD_general_power_table)

end
