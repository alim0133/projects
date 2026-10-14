
function result = main_deseasonalised(abp_file, start_quarter, end_quarter, alpha)

    %% READ ABP DATA

    % ABP data contain:
    % quarter = quarter, e.g. "1975q1"
    % msacode = MSA code
    % cpi     = Consumer Price Index
    % hp      = Freddie Mac House Price Index

    data = readtable(abp_file, 'TextType', 'string');

    %% CONSTRUCT QUARTER INDEX

    % ABP quarter labels are stored as strings such as "1975q1".

    quarter_string = lower(string(data.quarter));

    year_number = str2double(extractBefore(quarter_string, "q"));
    quarter_number = str2double(extractAfter(quarter_string, "q"));

    % Numeric ABP quarter index used for sorting and date selection.

    data.quarter_index = 4*year_number + quarter_number;

    %% PARSE REQUESTED DATE RANGE

    start_string = lower(string(start_quarter));
    end_string = lower(string(end_quarter));

    start_y = str2double(extractBefore(start_string, "q"));
    start_q = str2double(extractAfter(start_string, "q"));

    end_y = str2double(extractBefore(end_string, "q"));
    end_q = str2double(extractAfter(end_string, "q"));

    if isnan(start_y) || isnan(start_q) || start_q < 1 || start_q > 4
        error('start_quarter must have the form "YYYYQ1", ..., "YYYYQ4".')
    end

    if isnan(end_y) || isnan(end_q) || end_q < 1 || end_q > 4
        error('end_quarter must have the form "YYYYQ1", ..., "YYYYQ4".')
    end

    start_index = 4*start_y + start_q;
    end_index = 4*end_y + end_q;

    if start_index > end_index
        error('start_quarter must not be after end_quarter.')
    end

    %% IDENTIFY MSAs AND QUARTERS

    % Unique MSA codes in chronological order.

    msa_codes = unique(data.msacode, 'sorted');

    % Unique ABP quarters in chronological order.

    quarter_indices = unique(data.quarter_index, 'sorted');

    % Construct readable quarter labels.

    year_numbers = floor((quarter_indices - 1)/4);
    quarter_numbers = mod(quarter_indices - 1, 4) + 1;

    quarter_labels = string(year_numbers) + "Q" + string(quarter_numbers);

    N = length(msa_codes);
    T_level = length(quarter_indices);

    % Check requested range against available ABP data.

    if start_index < quarter_indices(1) || end_index > quarter_indices(end)
        error('Requested date range lies outside the available ABP data.')
    end

    %% CONSTRUCT N x T ABP DATA MATRICES

    % HP(i,t)  contains the nominal house-price index for MSA i at time t.
    % CPI(i,t) contains the CPI matched to MSA i at time t.

    HP = NaN(N, T_level);
    CPI = NaN(N, T_level);

    % Map every ABP observation to its MSA and quarter.

    [~, msa_index] = ismember(data.msacode, msa_codes);
    [~, time_index] = ismember(data.quarter_index, quarter_indices);

    for row = 1:height(data)

        i = msa_index(row);
        t = time_index(row);

        HP(i,t) = data.hp(row);
        CPI(i,t) = data.cpi(row);

    end

    % Check that the house-price and CPI panels are balanced.

    if any(isnan(HP), 'all') || any(isnan(CPI), 'all')
        error('The ABP house-price/CPI panel contains missing observations.')
    end

    % Logarithms require strictly positive values.

    if any(HP <= 0, 'all') || any(CPI <= 0, 'all')
        error('House-price and CPI observations must be strictly positive.')
    end

    %% CONSTRUCT QUARTERLY REAL HOUSE-PRICE CHANGES

    % Real house-price change:
    % Pi_it = 100 * [log(HP_it/CPI_it) - log(HP_i,t-1/CPI_i,t-1)]

    log_real_price = log(HP) - log(CPI);

    Pi_full = 100 * diff(log_real_price, 1, 2);

    % Each change dated quarter t uses levels from t and t-1.

    change_indices = quarter_indices(2:end);
    change_labels = quarter_labels(2:end);

    %% SELECT DATE RANGE

    % Select the requested range after constructing changes so that, for
    % example, a change starting in 2000Q1 can use the 1999Q4 level.

    if start_index < change_indices(1)
        error(['The first available real house-price change is 1975Q2. ', ...
               'Choose start_quarter >= "1975Q2".'])
    end

    if end_index > change_indices(end)
        error(['The final available house-price change is 2014Q4. ', ...
               'Choose end_quarter <= "2014Q4".'])
    end

    selected = change_indices >= start_index & ...
               change_indices <= end_index;

    if ~any(selected)
        error('No real house-price changes are available in the requested range.')
    end

    Y = Pi_full(:, selected);

    selected_indices = change_indices(selected);
    selected_quarters = change_labels(selected);

    % Effective time dimension after first differencing and date selection.

    T = size(Y, 2);

    % Quarter number for every selected observation.

    selected_quarter_numbers = mod(selected_indices - 1, 4) + 1;

    %% CONSTRUCT QUARTERLY SEASONAL DUMMIES

    % Q1 is the omitted reference quarter.
    %
    % D2 = 1 for Q2
    % D3 = 1 for Q3
    % D4 = 1 for Q4

    D2 = selected_quarter_numbers == 2;
    D3 = selected_quarter_numbers == 3;
    D4 = selected_quarter_numbers == 4;

    %% DE-SEASONALISE REAL HOUSE-PRICE CHANGES

    % For each MSA i estimate:
    %
    % Pi_it = a_i ...
    %       + delta_i2 D2_t ...
    %       + delta_i3 D3_t ...
    %       + delta_i4 D4_t ...
    %       + e_it
    %
    % e_it is the de-seasonalised house-price change.

    residuals = zeros(N, T);

    intercepts = zeros(N, 1);
    seasonal_coefficients = zeros(N, 3);

    % Common regressors for every MSA.

    X_seasonal = [ones(T,1), D2, D3, D4]; % T x 4

    for i = 1:N

        % Real house-price changes for MSA i.

        y_i = Y(i,:)'; % T x 1

        % OLS coefficients:
        % beta_i = [a_i; delta_i2; delta_i3; delta_i4; gamma_i]

        beta_i = X_seasonal \ y_i; % 4 x 1

        intercepts(i) = beta_i(1); % N x 1
        seasonal_coefficients(i,:) = beta_i(2:4)'; % N x 3

        % De-seasonalised and de-factored residuals.

        residuals(i,:) = (y_i - X_seasonal*beta_i)'; % 1 x T

    end

    %% PREPARE RESIDUALS FOR CD STATISTICS

    % The regressions contain an intercept, so the residual means are
    % already zero up to numerical precision.
    %
    % Demean again explicitly so the implementation matches the pairwise
    % definitions used in the CD statistics.

    E = residuals - mean(residuals, 2);

    %% NUMBER OF UNIQUE PAIRS

    % Number of pairs (i,j), where i < j.
    
    number_pairs = N*(N-1)/2;

    pearson_stats = zeros(number_pairs, 1);
    giraitis_stats = zeros(number_pairs, 1);

    %% CONSTRUCT SCALED RESIDUALS FOR GENERAL CD

    % For each MSA:
    % sum_squares(i) = sum_t e_it^2

    sum_squares = sum(E.^2, 2);

    if any(sum_squares <= eps)
        error('A residual series has numerically zero variance.')
    end

    % Scaled residual:
    % xi_it = e_it / sqrt(sum_t e_it^2)

    Xi = E ./ sqrt(sum_squares);

    %% COMPUTE PAIRWISE STATISTICS

    pair_counter = 0;

    % General CD numerator and squared denominator.

    general_numerator = 0;
    general_denominator_squared = 0;

    for i = 1:N-1

        e_i = E(i,:);
        xi_i = Xi(i,:);

        for j = i+1:N

            pair_counter = pair_counter + 1;

            e_j = E(j,:);
            xi_j = Xi(j,:);

            %% PAIRWISE PEARSON STATISTIC

            pearson_numerator = sum(e_i .* e_j);
            pearson_denominator = sqrt(sum(e_i.^2) * sum(e_j.^2));

            if pearson_denominator <= eps
                error(['The Pearson denominator is numerically ', ...
                       'zero for pair (%d,%d).'], i, j)
            end
            
            % t_ij = sqrt(T) * rho_hat_ij
            %      = sqrt(T) * sum_t e_it e_jt / sqrt(sum_t e_it^2 * sum_t e_jt^2)

            pearson_stats(pair_counter) = sqrt(T) * pearson_numerator / pearson_denominator;

            %% PAIRWISE GIRAITIS STATISTIC

            cross_products = e_i .* e_j;

            giraitis_numerator = sum(cross_products);
            giraitis_denominator = sqrt(sum(cross_products.^2));

            if giraitis_denominator <= eps
                error(['The Giraitis denominator is numerically ', ...
                       'zero for pair (%d,%d).'], i, j)
            end
            
            % ttilde_ij = sum_t e_it e_jt / sqrt(sum_t e_it^2 e_jt^2)

            giraitis_stats(pair_counter) = giraitis_numerator / giraitis_denominator;

            %% GENERAL CD COMPONENTS

            % Scaled cross-products:
            % xi_it xi_jt.

            scaled_cross_products = xi_i .* xi_j;

            % Add sum_t xi_it xi_jt to the overall numerator.

            general_numerator = general_numerator + sum(scaled_cross_products);

            % Add sum_t xi_it^2 xi_jt^2 to the squared denominator.

            general_denominator_squared = general_denominator_squared + sum(scaled_cross_products.^2);

        end

    end

    %% STANDARD CD STATISTIC

    % CD = sqrt(2/N(N-1)) * sum_{i<j} t_ij

    cd_stat = sqrt(2 / (N*(N-1))) * sum(pearson_stats);

    %% ROBUST CD STATISTIC

    % CD_robust = sqrt(2/N(N-1)) * sum_{i<j} ttilde_ij

    cd_robust_stat = sqrt(2 / (N*(N-1))) * sum(giraitis_stats);

    %% GENERAL CD STATISTIC

    % CD_general = sum_{i<j} sum_t xi_it xi_jt / sqrt[sum_{i<j} sum_t xi_it^2 xi_jt^2]

    if general_denominator_squared <= eps
        error('The general CD denominator is numerically zero.')
    end

    cd_general_stat = general_numerator / sqrt(general_denominator_squared);

    %% REJECTION DECISIONS

    % Two-sided standard-normal critical value.

    critical_value = norminv(1 - alpha/2);

    % Reject H0 if:
    % |statistic| > z_(1-alpha/2)

    cd_rejection = abs(cd_stat) > critical_value;
    cd_robust_rejection = abs(cd_robust_stat) > critical_value;
    cd_general_rejection = abs(cd_general_stat) > critical_value;

    %% P-VALUES

    cd_pvalue = 2 * (1 - normcdf(abs(cd_stat)));
    cd_robust_pvalue = 2 * (1 - normcdf(abs(cd_robust_stat)));
    cd_general_pvalue = 2 * (1 - normcdf(abs(cd_general_stat)));

    %% DISPLAY RESULTS

    fprintf('\nABP EMPIRICAL APPLICATION\n')
    fprintf('-------------------------\n')

    fprintf('Selected range: %s to %s\n', ...
        char(selected_quarters(1)), ...
        char(selected_quarters(end)))

    fprintf('Number of MSAs: %d\n', N)
    fprintf('Number of quarterly changes: %d\n', T)
    fprintf('Number of MSA pairs: %d\n\n', number_pairs)

    fprintf('Standard CD statistic: %.4f\n', cd_stat)
    fprintf('Standard CD p-value: %.6f\n', cd_pvalue)
    fprintf('Standard CD rejection: %d\n\n', cd_rejection)

    fprintf('Robust CD statistic: %.4f\n', cd_robust_stat)
    fprintf('Robust CD p-value: %.6f\n', cd_robust_pvalue)
    fprintf('Robust CD rejection: %d\n\n', cd_robust_rejection)

    fprintf('General CD statistic: %.4f\n', cd_general_stat)
    fprintf('General CD p-value: %.6f\n', cd_general_pvalue)
    fprintf('General CD rejection: %d\n', cd_general_rejection)


    %% STORE RESULTS

    result.N = N;
    result.T = T;
    result.number_pairs = number_pairs;

    result.msa_codes = msa_codes;
    result.quarters = selected_quarters;

    % Raw quarterly real house-price changes.

    result.real_house_price_changes = Y;

    % Seasonal-filtering coefficients and de-seasonalised series.

    result.intercepts = intercepts;
    result.seasonal_coefficients = seasonal_coefficients;
    result.residuals = residuals;

    % Pairwise quantities used in the CD statistics.

    result.pearson_stats = pearson_stats;
    result.giraitis_stats = giraitis_stats;
    result.scaled_residuals = Xi;

    % CD statistics.

    result.cd_stat = cd_stat;
    result.cd_robust_stat = cd_robust_stat;
    result.cd_general_stat = cd_general_stat;

    result.cd_pvalue = cd_pvalue;
    result.cd_robust_pvalue = cd_robust_pvalue;
    result.cd_general_pvalue = cd_general_pvalue;

    result.cd_rejection = cd_rejection;
    result.cd_robust_rejection = cd_robust_rejection;
    result.cd_general_rejection = cd_general_rejection;

end
