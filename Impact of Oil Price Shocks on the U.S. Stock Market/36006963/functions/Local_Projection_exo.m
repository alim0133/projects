function [IRF,se] = Local_Projection_exo(y,epsilon,w1,w2,p,hor,varargin)
% Benjamin Wong
% Monash University
% Edited June 2020
% Estimates a Local Projection assuming we have an exogenous variable

%% INPUTS
%
%y                        Independent Variable (LHS)
%epsilon                  the shock of interest
%w1                       Set of controls which include contemporaneous
%                         term (i.e. t)
%w2                       Set of controls which only has lags
%hor                      The number of horizons to estimate the Local
%                         Projectin
%differenced              if LHS variable is differenced to get back the
%                         level
%varargin                 enter "Differenced" if variable is in change on
%                         LHS so IRF is in the level

%% OUTPUTS
%IRF                      The estimated impulse response function
%se                       Vector of standard errors constructed using
%                           Newey-West standard errors

%% Preliminaries

%Ensure controls, shocks, and LHS are of the same size

if sum(strcmp(varargin,'Differenced'))>0
    if size(y,1)-size(w2,1) ~= 1
        error('Check Dimension of Matrices')
    end
else
    if size(y,1) ~= size(w2,1)
        error('Check Dimension of Matrices')
    end
end

if size(w2,1) ~= size(epsilon,1)
    error('Check Dimension of Matrices')
end


if size(w1,1) ~= 0      %in case the set in w1 is empty
    if size(w2,1) ~= size(w1,1)
        error('Check Dimension of Matrices')
    end
end



%% Preliminaries

T = size(w2,1)-p;

Y = y(p+1:end,:);     %Cut Away first p lags



%Create Lag Matrices

X = [];     %create lags of w1 and w2

if size(w1,1) ~=0
    X = w1(p+1:end,:);
end

for i = 1:p
    if size(w1,1) ~=0
        Z = w1(p+1-i:end-i,:);
        X = [X Z];
    end
    
    Z = w2(p+1-i:end-i,:);
    X = [X Z];
end

%Put a constant and the exogenous variable at the front of the matrix X
X = [ones(T,1) epsilon(p+1:end,:) X];



%%
IRF = NaN(1,hor+1);
se = NaN(1,hor+1);

%% Run Sequence of regressions

%Coefficient of interest is the second variable in X

for jj = 0:hor
    
    if sum(strcmp(varargin,'Differenced'))>0
        Yd = Y(jj+2:end) - Y(1:end-jj-1);
        
        B = inv(X(1:end-jj,:)'*X(1:end-jj,:))*(X(1:end-jj,:)'*Yd);
        IRF(1,jj+1) = B(2); %Second coefficient is of interest (the "epsilon")
        
        % Residual
        e = Yd - X(1:end-jj,:)*B;
    else
        B = inv(X(1:end-jj,:)'*X(1:end-jj,:))*(X(1:end-jj,:)'*Y(jj+1:end,:));
        IRF(1,jj+1) = B(2); %Second coefficient is of interest (the "epsilon")
        
        % Residual
        e = Y(jj+1:end,:) - X(1:end-jj,:)*B;
    end
    
    nobs = T-jj;
    Z = X(1:end-jj,:);
    
    %Newey West Standard errors
    Q = 0;
    for l = 0:hor   %This portion of the loop is to calculate the autocorrelation for the Newey-West std error
        w_l = 1-l/(hor+1);
        for t = l+1:nobs
            if (l==0)   % This calculates the S_0 portion
                Q = Q  + e(t) ^2 * Z(t, :)' * Z(t,:);
            else        % This calculates the off-diagonal terms
                Q = Q + w_l * e(t) * e(t-l)* ...
                    (Z(t, :)' * Z(t-l,:) + Z(t-l, :)' * Z(t,:));
            end
        end
    end
    Q = (1/(nobs-size(X,2))) .*Q;
    covar_mat = nobs.*((X(1:end-jj,:)'*X(1:end-jj,:))\Q/(X(1:end-jj,:)'*X(1:end-jj,:)));
    se(1,jj+1) = sqrt(covar_mat(2,2));
    
end
end