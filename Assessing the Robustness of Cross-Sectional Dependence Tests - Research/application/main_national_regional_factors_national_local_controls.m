
function result = main_national_regional_factors_national_local_controls(abp_file, gdpc1_file, start_quarter, end_quarter, alpha)

    %% READ DATA

    % ABP data contain:
    % quarter    = quarter, e.g. "1975q1"
    % msacode    = MSA code
    % regioncode = BEA region
    % cpi        = Consumer Price Index
    % hp         = Freddie Mac House Price Index
    % pop        = population
    % pipc       = nominal income per capita

    data = readtable(abp_file, 'TextType', 'string');

    % FRED GDPC1 data contain:
    % observation_date = first day of each quarter
    % GDPC1            = U.S. real GDP level

    gdp_data = readtable(gdpc1_file, 'TextType', 'string');

    %% CONSTRUCT QUARTER INDEX

    % ABP quarter labels are stored as strings such as "1975q1".

    quarter_string = lower(string(data.quarter));

    year_number = str2double(extractBefore(quarter_string, "q"));
    quarter_number = str2double(extractAfter(quarter_string, "q"));

    % Numeric ABP quarter index used for sorting and date selection.

    data.quarter_index = 4*year_number + quarter_number;

    % FRED GDPC1 dates are stored as dates such as "1975-01-01".

    gdp_dates = datetime( ...
        string(gdp_data.observation_date), 'InputFormat', 'yyyy-MM-dd');

    gdp_year_number = year(gdp_dates);
    gdp_quarter_number = floor((month(gdp_dates) - 1)/3) + 1;

    % Numeric GDP quarter index.

    gdp_data.quarter_index = 4*gdp_year_number + gdp_quarter_number;

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

    % Construct readable ABP quarter labels.

    year_numbers = floor((quarter_indices - 1)/4);
    quarter_numbers = mod(quarter_indices - 1, 4) + 1;

    quarter_labels = string(year_numbers) + "Q" + string(quarter_numbers);

    N = length(msa_codes);
    T_level = length(quarter_indices);

    % Check requested range against available ABP data.

    if start_index < quarter_indices(1) || end_index > quarter_indices(end)
        error('Requested date range lies outside the available ABP data.')
    end

    % Unique GDP quarters in chronological order.
    
    [gdp_quarter_indices, gdp_sort_order] = unique(gdp_data.quarter_index, 'sorted');

    gdp_levels = gdp_data.GDPC1(gdp_sort_order);

    % Real GDP levels must be observed and strictly positive before taking logs.

    if any(isnan(gdp_levels)) || any(gdp_levels <= 0)
        error('GDPC1 contains missing or non-positive observations.')
    end

    %% CONSTRUCT N x T ABP DATA MATRICES

    % HP(i,t) contains the nominal house-price index for MSA i at time t.
    % CPI(i,t) contains the CPI matched to MSA i at time t.
    % POP(i,t)  contains population for MSA i at time t.
    % PIPC(i,t) contains nominal income per capita for MSA i at time t.

    HP = NaN(N, T_level);
    CPI = NaN(N, T_level);
    POP = NaN(N, T_level);
    PIPC = NaN(N, T_level);

    % Map every ABP observation to its MSA and quarter.

    [~, msa_index] = ismember(data.msacode, msa_codes);
    [~, time_index] = ismember(data.quarter_index, quarter_indices);

    for row = 1:height(data)

        i = msa_index(row);
        t = time_index(row);

        HP(i,t) = data.hp(row);
        CPI(i,t) = data.cpi(row);
        POP(i,t) = data.pop(row);
        PIPC(i,t) = data.pipc(row);

    end

    % Check that the four ABP panels are balanced.

    if any(isnan(HP), 'all') || any(isnan(CPI), 'all') || ...
       any(isnan(POP), 'all') || any(isnan(PIPC), 'all')

        error('The ABP panel contains missing observations.')
    end

    % Logarithms require strictly positive values.

    if any(HP <= 0, 'all') || any(CPI <= 0, 'all') || ...
       any(POP <= 0, 'all') || any(PIPC <= 0, 'all')

        error('HP, CPI, population and income observations must be strictly positive.')
    end

    %% IDENTIFY BEA REGION FOR EACH MSA

    % The region code should be constant over time for a given MSA.
    
    msa_region_codes = NaN(N, 1);
    
    for i = 1:N
    
        region_code = unique(data.regioncode(msa_index == i));
    
        if numel(region_code) ~= 1 || isnan(region_code(1))
            error('Each MSA must be associated with exactly one valid region code.')
        end
    
        msa_region_codes(i) = region_code(1);
    
    end
    
    % Unique BEA region codes represented in the sample.
    
    region_codes = unique(msa_region_codes, 'sorted');
    
    R = length(region_codes);

    %% MATCH U.S. REAL GDP LEVELS TO ABP QUARTERS

    % Match each GDPC1 quarter to each ABP quarter.

    [is_match, gdp_index] = ismember(quarter_indices, gdp_quarter_indices);

    if any(~is_match)
        error('GDPC1 is missing one or more quarters required by the ABP sample.')
    end

    % Store the matched GDPC1 levels.

    GDP = gdp_levels(gdp_index)'; % 1 x T_level

    %% CONSTRUCT QUARTERLY GROWTH RATES

    % Real house-price change:
    % Pi_it = 100 * [log(HP_it/CPI_it) - log(HP_i,t-1/CPI_i,t-1)]

    log_real_price = log(HP) - log(CPI);

    Pi_full = 100 * diff(log_real_price, 1, 2);

    % Population growth:
    % GPOP_it = 100 * [log(POP_it) - log(POP_i,t-1)]

    GPOP_full = 100 * diff(log(POP), 1, 2);

    % Real per-capita income growth:
    % GINC_it = 100 * [log(PIPC_it/CPI_it) - log(PIPC_i,t-1/CPI_i,t-1)]

    log_real_income = log(PIPC) - log(CPI);

    GINC_full = 100 * diff(log_real_income, 1, 2);

    % U.S. real GDP growth:
    % GGDP_t = 100 * [log(RGDP_t) - log(RGDP_t-1)]

    GGDP_full = 100 * diff(log(GDP));

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
    GPOP = GPOP_full(:, selected);
    GINC = GINC_full(:, selected);
    GGDP = GGDP_full(selected);

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

    %% ESTIMATE NATIONAL AND REGIONAL FACTORS FOR EACH VARIABLE

    % Each variable is de-factored using its OWN national and regional
    % cross-sectional averages.
    %
    % House prices:
    %   national Pi factor  = average Pi across all MSAs
    %   regional Pi factor  = average Pi across MSAs in region r
    %
    % Population:
    %   national GPOP factor = average GPOP across all MSAs
    %   regional GPOP factor = average GPOP across MSAs in region r
    %
    % Income:
    %   national GINC factor = average GINC across all MSAs
    %   regional GINC factor = average GINC across MSAs in region r

    pi_national_factor = mean(Y, 1);
    gpop_national_factor = mean(GPOP, 1);
    ginc_national_factor = mean(GINC, 1);

    pi_regional_averages = zeros(R, T);
    gpop_regional_averages = zeros(R, T);
    ginc_regional_averages = zeros(R, T);

    pi_regional_factor = zeros(N, T);
    gpop_regional_factor = zeros(N, T);
    ginc_regional_factor = zeros(N, T);

    for r = 1:R

        % Identify MSAs belonging to region r.

        region_members = msa_region_codes == region_codes(r);
        N_r = sum(region_members);

        if N_r == 0
            error('A region contains no MSAs.')
        end

        % Regional cross-sectional averages.

        pi_regional_averages(r,:) = mean(Y(region_members,:), 1);
        gpop_regional_averages(r,:) = mean(GPOP(region_members,:), 1);
        ginc_regional_averages(r,:) = mean(GINC(region_members,:), 1);

        % Give every MSA in region r its corresponding regional factor.

        pi_regional_factor(region_members,:) = repmat(pi_regional_averages(r,:), N_r, 1);
        gpop_regional_factor(region_members,:) = repmat(gpop_regional_averages(r,:), N_r, 1);
        ginc_regional_factor(region_members,:) = repmat(ginc_regional_averages(r,:), N_r, 1);

    end

    %% DE-SEASONALISE AND DE-FACTOR REAL HOUSE-PRICE CHANGES

    % For MSA i belonging to region r, estimate:
    %
    % Pi_it = a_i ...
    %       + delta_i2 D2_t ...
    %       + delta_i3 D3_t ...
    %       + delta_i4 D4_t ...
    %       + gamma_N,i Pibar_t ...
    %       + gamma_R,i Pibar_rt ...
    %       + pi_it
    %
    % pi_it is the de-seasonalised and de-factored house-price change.

    pi = zeros(N, T);

    pi_intercepts = zeros(N, 1);
    pi_seasonal_coefficients = zeros(N, 3);
    pi_national_factor_loadings = zeros(N, 1);
    pi_regional_factor_loadings = zeros(N, 1);

    for i = 1:N

        % Real house-price changes for MSA i.

        y_i = Y(i,:)'; % T x 1

        % MSA-specific design matrix.
    
        X_pi_i = [ones(T,1), D2, D3, D4, pi_national_factor', pi_regional_factor(i,:)']; % T x 6

        % OLS coefficients:
        % beta_i = [a_i; delta_i2; delta_i3; delta_i4; gamma_N,i; gamma_R,i]

        beta_pi_i = X_pi_i \ y_i; % 6 x 1

        pi_intercepts(i) = beta_pi_i(1); % N x 1
        pi_seasonal_coefficients(i,:) = beta_pi_i(2:4)'; % N x 3
        pi_national_factor_loadings(i) = beta_pi_i(5); % N x 1
        pi_regional_factor_loadings(i) = beta_pi_i(6);

        % De-seasonalised and de-factored residuals.

        pi(i,:) = (y_i - X_pi_i*beta_pi_i)'; % 1 x T

    end

    %% DE-SEASONALISE AND DE-FACTOR POPULATION GROWTH

    % For MSA i belonging to region r, estimate:
    %
    % GPOP_it = a_pop,i ...
    %         + delta_pop,i2 D2_t ...
    %         + delta_pop,i3 D3_t ...
    %         + delta_pop,i4 D4_t ...
    %         + gamma_pop,N,i GPOPbar_t ...
    %         + gamma_pop,R,i GPOPbar_rt ...
    %         + gpop_it
    %
    % gpop_it is the de-seasonalised and de-factored population growth.

    gpop = zeros(N, T);

    gpop_intercepts = zeros(N, 1);
    gpop_seasonal_coefficients = zeros(N, 3);
    gpop_national_factor_loadings = zeros(N, 1);
    gpop_regional_factor_loadings = zeros(N, 1);

    for i = 1:N

        gpop_i = GPOP(i,:)';

        X_gpop_i = [ones(T,1), D2, D3, D4, gpop_national_factor', gpop_regional_factor(i,:)'];

        beta_gpop_i = X_gpop_i \ gpop_i;

        gpop_intercepts(i) = beta_gpop_i(1);
        gpop_seasonal_coefficients(i,:) = beta_gpop_i(2:4)';
        gpop_national_factor_loadings(i) = beta_gpop_i(5);
        gpop_regional_factor_loadings(i) = beta_gpop_i(6);

        gpop(i,:) = (gpop_i - X_gpop_i*beta_gpop_i)';

    end

    %% DE-SEASONALISE AND DE-FACTOR REAL PER-CAPITA INCOME GROWTH

    % For MSA i belonging to region r, estimate:
    %
    % GINC_it = a_inc,i ...
    %         + delta_inc,i2 D2_t ...
    %         + delta_inc,i3 D3_t ...
    %         + delta_inc,i4 D4_t ...
    %         + gamma_inc,N,i GINCbar_t ...
    %         + gamma_inc,R,i GINCbar_rt ...
    %         + ginc_it
    %
    % ginc_it is the de-seasonalised and de-factored real per-capita income growth.

    ginc = zeros(N, T);

    ginc_intercepts = zeros(N, 1);
    ginc_seasonal_coefficients = zeros(N, 3);
    ginc_national_factor_loadings = zeros(N, 1);
    ginc_regional_factor_loadings = zeros(N, 1);

    for i = 1:N

        ginc_i = GINC(i,:)';

        X_ginc_i = [ones(T,1), D2, D3, D4, ginc_national_factor', ginc_regional_factor(i,:)'];

        beta_ginc_i = X_ginc_i \ ginc_i;

        ginc_intercepts(i) = beta_ginc_i(1);
        ginc_seasonal_coefficients(i,:) = beta_ginc_i(2:4)';
        ginc_national_factor_loadings(i) = beta_ginc_i(5);
        ginc_regional_factor_loadings(i) = beta_ginc_i(6);

        ginc(i,:) = (ginc_i - X_ginc_i*beta_ginc_i)';

    end
    
    %% REGRESS DE-FACTORED HOUSE-PRICE CHANGES ON DE-FACTORED REGRESSORS

    % For each MSA i estimate:
    %
    % pi_it = c_i ...
    %       + beta_pop,i * gpop_it ...
    %       + beta_inc,i * ginc_it ...
    %       + beta_gdp,i * GGDP_t ...
    %       + u_it
    %
    % GGDP_t is common across all MSAs, while its coefficient is allowed
    % to differ across MSAs.
    %
    % The residual e_it is used in the three CD statistics.

    residuals = zeros(N, T);

    regression_intercepts = zeros(N, 1);
    population_coefficients = zeros(N, 1);
    income_coefficients = zeros(N, 1);
    gdp_coefficients = zeros(N, 1);

    for i = 1:N

        pi_i = pi(i,:)'; % T x 1

        Z_i = [ones(T,1), gpop(i,:)', ginc(i,:)', GGDP']; % T x 4

        theta_i = Z_i \ pi_i; % 4 x 1

        regression_intercepts(i) = theta_i(1);
        population_coefficients(i) = theta_i(2);
        income_coefficients(i) = theta_i(3);
        gdp_coefficients(i) = theta_i(4);

        residuals(i,:) = (pi_i - Z_i*theta_i)'; % 1 x T

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
    fprintf('Number of BEA regions: %d\n', R)
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
    result.R = R;
    result.T = T;
    result.number_pairs = number_pairs;

    result.msa_codes = msa_codes;
    result.region_codes = region_codes;
    result.msa_region_codes = msa_region_codes;
    result.quarters = selected_quarters;

    % Raw quarterly growth-rate panels used in the selected sample.

    result.real_house_price_changes = Y;
    result.population_growth = GPOP;
    result.real_income_growth = GINC;
    result.real_gdp_growth = GGDP;

    % National factors.

    result.pi_national_factor = pi_national_factor;
    result.gpop_national_factor = gpop_national_factor;
    result.ginc_national_factor = ginc_national_factor;

    % Regional factors and regional averages.

    result.pi_regional_factor = pi_regional_factor;
    result.gpop_regional_factor = gpop_regional_factor;
    result.ginc_regional_factor = ginc_regional_factor;

    result.pi_regional_averages = pi_regional_averages;
    result.gpop_regional_averages = gpop_regional_averages;
    result.ginc_regional_averages = ginc_regional_averages;

    % De-seasonalised and de-factored series.

    result.pi = pi;
    result.gpop = gpop;
    result.ginc = ginc;

    % House-price filtering coefficients.

    result.pi_intercepts = pi_intercepts;
    result.pi_seasonal_coefficients = pi_seasonal_coefficients;
    result.pi_national_factor_loadings = pi_national_factor_loadings;
    result.pi_regional_factor_loadings = pi_regional_factor_loadings;

    % Population-growth filtering coefficients.

    result.gpop_intercepts = gpop_intercepts;
    result.gpop_seasonal_coefficients = gpop_seasonal_coefficients;
    result.gpop_national_factor_loadings = gpop_national_factor_loadings;
    result.gpop_regional_factor_loadings = gpop_regional_factor_loadings;

    % Real-income-growth filtering coefficients.

    result.ginc_intercepts = ginc_intercepts;
    result.ginc_seasonal_coefficients = ginc_seasonal_coefficients;
    result.ginc_national_factor_loadings = ginc_national_factor_loadings;
    result.ginc_regional_factor_loadings = ginc_regional_factor_loadings;

    % Final regression coefficients and residuals.

    result.regression_intercepts = regression_intercepts;
    result.population_coefficients = population_coefficients;
    result.income_coefficients = income_coefficients;
    result.gdp_coefficients = gdp_coefficients;
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