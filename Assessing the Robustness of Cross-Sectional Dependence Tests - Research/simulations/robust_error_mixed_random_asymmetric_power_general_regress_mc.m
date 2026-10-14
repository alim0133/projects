
function result = ...
    robust_error_mixed_random_asymmetric_power_general_regress_mc(N, T, reps, a, scenario, gamma_range, alpha)

    %% INPUT CHECK

    % 0 <= a <= 1 controls number of individuals with asymmetric
    % product errors.

    if a < 0 || a > 1

        error('a must satisfy 0 <= a <= 1.')

    end

    %% GROUP SIZES

    % Number of asymmetric product-error individuals.

    N_product = floor(N^a);

    % Ensure number remains between 1 and N.

    N_product = min(N, max(1, N_product));

    % Remaining individuals have iid normal innovations.

    N_iid = N - N_product;

    %% REGRESSION COEFFICIENTS

    % Panel regression:
    % y_it = beta_0 + beta_1 x_it + u_it
    %
    % where:
    % beta_0 = 0.2
    % beta_1 = 0.5

    beta_0 = 0.2;
    beta_1 = 0.5;

    %% PREALLOCATE STORAGE

    % Three statistics:
    % 1. Standard CD
    % 2. Robust CD based on Giraitis pair statistics
    % 3. Pesaran general CD

    cd_stats = zeros(reps,1);
    cd_robust_stats = zeros(reps,1);
    cd_general_stats = zeros(reps,1);

    % Rejection decisions.

    cd_rejections = zeros(reps,1);
    cd_robust_rejections = zeros(reps,1);
    cd_general_rejections = zeros(reps,1);

    % Store pairwise statistics.

    all_pearson_stats = cell(reps,1);
    all_giraitis_stats = cell(reps,1);

    % Store product-error units from each replication.

    all_product_units = cell(reps,1);

    %% FACTOR LOADINGS

    % Size:
    % gamma_i = 0
    %
    % Power:
    % gamma_i ~ U(gamma_lower,gamma_upper)
    %
    % gamma_i is drawn once and held fixed across all Monte Carlo replications.

    switch lower(scenario)

        case 'size'

            gamma = zeros(N,1);

        case 'power'

            gamma_lower = gamma_range(1);

            gamma_upper = gamma_range(2);

            gamma = gamma_lower + (gamma_upper - gamma_lower) * rand(N,1);

        otherwise

            error('Unknown scenario.')

    end

    %% MONTE CARLO REPLICATIONS

    for r = 1:reps

        %% RANDOMLY ASSIGN ERROR TYPES

        % Randomly order the N individuals.

        random_order = randperm(N);

        % First N_product individuals receive asymmetric product errors.

        product_units = random_order(1:N_product);

        % Remaining individuals receive iid errors.

        iid_units = random_order(N_product+1:N);

        % Store selected product-error units.

        all_product_units{r} = product_units;

        %% DRAW OBSERVED REGRESSOR

        % x_it ~ iid N(0,1).

        X = randn(N,T);

        %% DRAW PRIMITIVE I.I.D. NORMAL INNOVATIONS

        % epsilon_it ~ iid N(0,1).

        epsilon = randn(N,T);

        %% CONSTRUCT ASYMMETRIC STANDARDIZED INNOVATIONS

        % z_it = (epsilon_it^2 - 1) / sqrt(2)
        %
        % E(z_it)   = 0
        % Var(z_it) = 1.

        z = (epsilon.^2 - 1) / sqrt(2);

        %% CONSTRUCT CROSS-SECTIONAL LEAD

        % z_lead(i,t) = z_i+1,t
        %
        % with circular wrap-around:
        % z_lead(N,t) = z_1t.

        z_lead = circshift(z, -1, 1);

        %% CONSTRUCT MIXED IDIOSYNCRATIC INNOVATIONS

        % Start by assigning iid normal errors to every individual:
        % v_it = epsilon_it.

        V = epsilon;

        % Replace selected rows with asymmetric product errors:
        % v_it = z_it z_i+1,t.

        V(product_units,:) = z(product_units,:) .* z_lead(product_units,:);

        %% DRAW COMMON FACTOR

        % f_t ~ iid N(0,1).

        f = randn(1,T);

        %% CONSTRUCT FULL ERROR PROCESS

        % u_it = gamma_i f_t + v_it.

        U = gamma * f + V;

        %% GENERATE PANEL REGRESSION

        % y_it = beta_0 + beta_1 x_it + u_it.

        Y = beta_0 + beta_1 * X + U;

        %% ESTIMATE PANEL REGRESSIONS

        % Preallocate OLS residual matrix.

        E = zeros(N,T);

        for i = 1:N

            % Add an intercept.

            X_i = [ones(T,1), X(i,:)'];

            % Dependent variable for individual i.

            y_i = Y(i,:)';

            % OLS estimator:
            % beta_hat_i = (X_i'X_i)^(-1) X_i'y_i
            %
            % MATLAB backslash computes the OLS solution.

            beta_hat_i = X_i \ y_i;

            % OLS residuals:
            % e_hat_i = y_i - X_i beta_hat_i.

            E(i,:) = (y_i - X_i * beta_hat_i)';

        end

        %% NUMBER OF UNIQUE PAIRS

        % Number of pairs (i,j), where i < j.

        number_pairs = N*(N-1)/2;

        % Preallocate pairwise standard Pearson and Giraitis statistics.

        pearson_stats = zeros(number_pairs,1);
        giraitis_stats = zeros(number_pairs,1);

        %% CONSTRUCT SCALED RESIDUALS FOR GENERAL CD

        % For each individual:
        % sum_squares(i) = sum_t e_hat_it^2.

        sum_squares = sum(E.^2, 2);

        % Protect against a zero-variance residual series.

        if any(sum_squares <= eps)

            error('An OLS residual series has numerically zero variance.')

        end

        % Scaled residual:
        % xi_it = e_hat_it / sqrt(sum_t e_hat_it^2).

        Xi = E ./ sqrt(sum_squares);

        %% COMPUTE PAIRWISE STATISTICS

        % Counter identifies position of each pair.

        pair_counter = 0;

        % Preallocate general CD numerator and squared denominator.

        general_numerator = 0;
        general_denominator_squared = 0;

        for i = 1:N-1

            % OLS residuals for individual i.

            e_i = E(i,:);

            % Scaled OLS residuals for individual i.

            xi_i = Xi(i,:);

            for j = i+1:N

                pair_counter = pair_counter + 1;

                % OLS residuals for individual j.

                e_j = E(j,:);

                % Scaled OLS residuals for individual j.

                xi_j = Xi(j,:);

                %% PAIRWISE PEARSON STATISTIC

                % Pearson numerator:
                % sum_t e_hat_it e_hat_jt.

                pearson_numerator = sum(e_i .* e_j);

                % Pearson denominator:
                % sqrt(sum_t e_hat_it^2 * sum_t e_hat_jt^2).

                pearson_denominator = sqrt(sum(e_i.^2) * sum(e_j.^2));

                % Protect against a zero denominator.

                if pearson_denominator <= eps

                    error(['The Pearson denominator is numerically ', ...
                           'zero for pair (%d,%d).'], i, j)

                end

                % Pairwise Pearson statistic:
                %
                % t_ij = sqrt(T) * rho_hat_ij
                %      = sqrt(T) * sum_t e_hat_it e_hat_jt /
                %        sqrt(sum_t e_hat_it^2 * sum_t e_hat_jt^2).

                pearson_stats(pair_counter) = sqrt(T) * pearson_numerator / pearson_denominator;

                %% PAIRWISE GIRAITIS STATISTIC

                % Cross-products:
                % e_hat_it e_hat_jt.

                cross_products = e_i .* e_j;

                % Giraitis numerator:
                % sum_t e_hat_it e_hat_jt.

                giraitis_numerator = sum(cross_products);

                % Giraitis denominator:
                % sqrt(sum_t e_hat_it^2 e_hat_jt^2).

                giraitis_denominator = sqrt(sum(cross_products.^2));

                % Protect against a zero denominator.

                if giraitis_denominator <= eps

                    error(['The Giraitis denominator is numerically ', ...
                           'zero for pair (%d,%d).'], i, j)

                end

                % Pairwise Giraitis statistic:
                %
                % t_tilde_ij = sum_t e_hat_it e_hat_jt /
                %              sqrt(sum_t e_hat_it^2 e_hat_jt^2).

                giraitis_stats(pair_counter) = giraitis_numerator / giraitis_denominator;

                %% GENERAL CD COMPONENTS

                % Scaled cross-products:
                % xi_it xi_jt.

                scaled_cross_products = xi_i .* xi_j;

                % General CD numerator:
                % sum_{i<j} sum_t xi_it xi_jt.
                %
                % Add sum_t xi_it xi_jt to the overall general CD numerator.

                general_numerator = general_numerator + sum(scaled_cross_products);

                % General CD squared denominator:
                % sum_{i<j} sum_t xi_it^2 xi_jt^2.
                %
                % Add sum_t xi_it^2 xi_jt^2 to the squared general CD denominator.

                general_denominator_squared = general_denominator_squared + sum(scaled_cross_products.^2);

            end

        end

        %% STORE PAIRWISE STATISTICS

        % Store all pairwise Pearson statistics from replication r.

        all_pearson_stats{r} = pearson_stats;

        % Store all pairwise Giraitis statistics from replication r.

        all_giraitis_stats{r} = giraitis_stats;

        %% STANDARD CD STATISTIC

        % CD = sqrt(2/N(N-1)) sum_{i<j} t_ij.

        cd_stats(r) = sqrt(2 / (N*(N-1))) * sum(pearson_stats);

        %% ROBUST CD STATISTIC

        % CD_robust = sqrt(2/N(N-1)) sum_{i<j} t_tilde_ij.

        cd_robust_stats(r) = sqrt(2 / (N*(N-1))) * sum(giraitis_stats);

        %% GENERAL CD STATISTIC

        % CD_general = sum_{i<j} sum_t xi_it xi_jt / sqrt(sum_{i<j} sum_t xi_it^2 xi_jt^2)

        if general_denominator_squared <= eps

            error('The general CD denominator is numerically zero.')

        end

        cd_general_stats(r) = general_numerator / sqrt(general_denominator_squared);

        %% TWO-SIDED REJECTION RULES

        % Two-sided standard-normal critical value.

        critical_value = norminv(1 - alpha/2);

        % H0:
        % rho_ij = 0 for every i ~= j.
        %
        % Reject whenever:
        % |statistic| > z_(1-alpha/2).

        cd_rejections(r) = abs(cd_stats(r)) > critical_value;
        cd_robust_rejections(r) = abs(cd_robust_stats(r)) > critical_value;
        cd_general_rejections(r) = abs(cd_general_stats(r)) > critical_value;

    end

    %% MONTE CARLO SUMMARY

    % Under scenario = 'size':
    % gamma_i = 0
    %
    % so the rejection frequency is empirical size.
    %
    % Under scenario = 'power':
    % gamma_i ~= 0
    %
    % so the rejection frequency is empirical power.

    cd_rejection_rate = mean(cd_rejections);
    cd_robust_rejection_rate = mean(cd_robust_rejections);
    cd_general_rejection_rate = mean(cd_general_rejections);

    %% STORE RESULTS

    % Empirical size or power.

    result.cd_rejection_rate = cd_rejection_rate;
    result.cd_robust_rejection_rate = cd_robust_rejection_rate;
    result.cd_general_rejection_rate = cd_general_rejection_rate;

    % Pervasiveness parameter.

    result.a = a;

    % Number of iid and asymmetric product-error units.

    result.N_iid = N_iid;
    result.N_product = N_product;

    % Means.

    result.mean_cd = mean(cd_stats);
    result.mean_cd_robust = mean(cd_robust_stats);
    result.mean_cd_general = mean(cd_general_stats);

    % Variances.

    result.var_cd = var(cd_stats);
    result.var_cd_robust = var(cd_robust_stats);
    result.var_cd_general = var(cd_general_stats);

    % Store pairwise statistics.

    result.pearson_stats = all_pearson_stats;
    result.giraitis_stats = all_giraitis_stats;

    % Store CD statistics.

    result.cd_stats = cd_stats;
    result.cd_robust_stats = cd_robust_stats;
    result.cd_general_stats = cd_general_stats;

    % Store rejection decisions.

    result.cd_rejections = cd_rejections;
    result.cd_robust_rejections = cd_robust_rejections;
    result.cd_general_rejections = cd_general_rejections;

    % Store randomly selected asymmetric product-error units
    % from every replication.

    result.product_units = all_product_units;

    % Store factor loadings.

    result.gamma = gamma;

    % Store true regression coefficients.

    result.beta_0 = beta_0;
    result.beta_1 = beta_1;

end
