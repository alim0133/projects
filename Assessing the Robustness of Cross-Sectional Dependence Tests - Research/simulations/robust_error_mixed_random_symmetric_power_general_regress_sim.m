
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
% Number of product-error individuals:
% N_product = floor(N^a)

a_vals = [1];

% Number of Monte Carlo replications.

reps = 5000;

% Significance level.

alpha = 0.05;

% Factor loading range for power.

gamma_weak = [0.1 0.3];

%% PREALLOCATE SIZE STORAGE

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

%% PREALLOCATE NUMBER OF PRODUCT-ERROR UNITS

N_product_storage = ...
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
            % u_it = gamma_i f_t + v_it
            %
            % Product-error individuals:
            % v_it = epsilon_it epsilon_i+1,t
            %
            % Remaining individuals:
            % v_it = epsilon_it

            result = ...
                robust_error_mixed_random_symmetric_power_general_regress_mc( ...
                    N, T, reps, a, 'size', [0 0], alpha);

            % Store number of product-error individuals.

            N_product_storage(j,k) = ...
                result.N_product;

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
            % u_it = gamma_i f_t + v_it
            %
            % Product-error individuals:
            % v_it = epsilon_it epsilon_i+1,t
            %
            % Remaining individuals:
            % v_it = epsilon_it

            result = ...
                robust_error_mixed_random_symmetric_power_general_regress_mc( ...
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

%% DISPLAY NUMBER OF PRODUCT-ERROR UNITS

N_product_table = array2table( ...
    N_product_storage, ...
    'VariableNames', strcat("a", string(a_vals)), ...
    'RowNames', strcat("N", string(N_vals)));

disp(' ')
disp('NUMBER OF PRODUCT-ERROR UNITS')
disp(N_product_table)

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
