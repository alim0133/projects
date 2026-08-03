clear all
clc

addpath('functions')
addpath('figures')
addpath('data')

% Load data
rawdata = xlsread('data.xlsx','oil_market_data','B2:D398');
CPI = xlsread('data.xlsx','US_monthly_data','C318:C714');
SP500 = xlsread('sp500.xlsx','sp500','B2:B398');
yearlab =(1974+1/12:1/12:2006+12/12)'; % timeline

% Transform data
data = [100*diff(log(rawdata(:,1))) rawdata(2:end,3) 100*log(100*(rawdata(2:end,2)./CPI(2:end,1))) 100*(diff(log(SP500(:,1)))-diff(log(CPI(:,1))))];
data(:,1) = -1*data(:,1); % normalize to represent negative supply shock
variables = {'Oil Production (%)','Real Activity','Real Oil Price (%)', 'Cumulative Real Stock Returns (%)'};
shocks = {'Oil Supply Shock','Aggregate Demand Shock','Oil-Specific Demand Shock','Other Shocks to Stock Returns'};

p = 24; % lags
hor = 15; % IRF and FEVD horizon
CI = [0.05 0.16 0.84 0.95]'; % to build IRF confidence interval
n = size(data,2);

% Estimate VAR(24)
[phi,sigma,x,e,phi_boot,sigma_boot] = olsvar(data,p,'double_bootstrap');

% IRF and FEVD
A0 = chol(sigma,'lower'); % cholesky identification
[IRF,FEVD] = calculate_IRF_FEVD(phi,A0,hor,'FEVD');
IRF(1:n:end,:) = cumsum(IRF(1:n:end,:),2); % cumulate because oil production and real stock returns are first differenced
IRF(4:n:end,:) = cumsum(IRF(4:n:end,:),2);

% IRF confidence interval
IRF_boot = NaN(n^2,hor+1,1000);

for jj = 1:1000 % bootstrap after bootstrap
    A0_boot = chol(sigma_boot(:,:,jj),'lower');
    IRF_boot(:,:,jj) = calculate_IRF_FEVD(phi_boot(:,:,jj),A0_boot,hor);
    IRF_boot(1:n:end,:,jj) = cumsum(IRF_boot(1:n:end,:,jj),2);
    IRF_boot(4:n:end,:,jj) = cumsum(IRF_boot(4:n:end,:,jj),2);
end

IRF_quantiles = quantile(IRF_boot,CI,3);

% Expand y-axis to fit IRF
for jj = 1:n
    y_axis(jj,1) = 1.1*min(min(squeeze(IRF_quantiles(jj:n:end,:,1))));
    y_axis(jj,2) = 1.1*max(max(squeeze(IRF_quantiles(jj:n:end,:,4))));
end

% Graph IRF

%figure
%fig = gcf; fig.PaperOrientation = 'landscape';
%set(fig,'PaperUnits','normalized');
%set(fig,'PaperPosition', [0 0 1 1]);
%for jj = 1:n
%    for ii = 1:n
%        subplot(n,n,(jj-1)*n+ii)
%        plot(0:hor,IRF((jj-1)*n+ii,:),'-k','LineWidth',3)
%        hold on
%        plot([0 hor],zeros(2,1),'-k','LineWidth',1)
%        hold on
%        plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+ii,:,2:3)),'--b','LineWidth',1.5)
%        hold on
%        plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+ii,:,[1 4])),'--r','LineWidth',1.5)
%        ylim(y_axis(ii,:))
%        ylabel(variables{ii})
%        title(shocks{jj})
%        set(gca,'FontSize',11)
%        if jj == n
%            xlabel('Months')
%        end
%    end
%end
%print( [ pwd '\figures\SP500_IRF' ] , '-dpdf' );

figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
for jj = 1:n-1
    subplot(1,n-1,jj)
    plot(0:hor,IRF((jj-1)*n+n-1,:),'-k','LineWidth',3)
    hold on
    plot([0 hor],zeros(2,1),'-k','LineWidth',1)
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n-1,:,2:3)),'--b','LineWidth',1.5)
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n-1,:,[1 4])),'--r','LineWidth',1.5)
    ylim(y_axis(n-1,:))
    ylabel(variables{n-1})
    title(shocks{jj})
    set(gca,'FontSize',11)
    xlabel('Months')
end
print( [ pwd '\figures\SP500_IRF1' ] , '-dpdf' );

figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
for jj = 1:n-1
    subplot(1,n-1,jj)
    plot(0:hor,IRF((jj-1)*n+n,:),'-k','LineWidth',3)
    hold on
    plot([0 hor],zeros(2,1),'-k','LineWidth',1)
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n,:,2:3)),'--b','LineWidth',1.5)
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n,:,[1 4])),'--r','LineWidth',1.5)
    ylim(y_axis(n,:))
    ylabel(variables{n})
    title(shocks{jj})
    set(gca,'FontSize',11)
    xlabel('Months')
end
print( [ pwd '\figures\SP500_IRF2' ] , '-dpdf' );

% Historical decomposition
HD = historical_decomposition(phi,e,A0);

% Graph historical decomposition

%figure
%fig = gcf; fig.PaperOrientation = 'landscape';
%set(fig,'PaperUnits','normalized');
%set(fig,'PaperPosition', [0 0 1 1]);
%subplot(4,1,1)
%plot(yearlab,data(:,3),'-r','LineWidth',3)
%hold on
%plot([yearlab(1) yearlab(end)],zeros(2,1),'-k')
%ylim([0.95*min(data(:,3)) 1.05*max(data(:,3))])
%xlim([yearlab(1) yearlab(end)])
%set(gca,'FontSize',11)
%title('100 times log Real Price of Oil in 1983 US$')
%for jj = 1:3
%    subplot(4,1,jj+1)
%    plot(yearlab,HD{3,1}(:,jj),'-r','LineWidth',3)
%    hold on
%    plot([yearlab(1) yearlab(end)],zeros(2,1),'-k')
%    ylim([1.05*min(min(HD{3,1})) 1.05*max(max(HD{3,1}))])
%    xlim([yearlab(1) yearlab(end)])
%    set(gca,'FontSize',11)
%    title(['Contribution of ',shocks{jj},' to the Real Price of Oil'])
%end
%print( [ pwd '\figures\SP500_HD' ] , '-dpdf' );

figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
for jj = 1:3
    subplot(3,1,jj)
    plot(yearlab,HD{3,1}(:,jj),'-r','LineWidth',3)
    hold on
    plot([yearlab(1) yearlab(end)],zeros(2,1),'-k')
    ylim([1.05*min(min(HD{3,1})) 1.05*max(max(HD{3,1}))])
    xlim([yearlab(1) yearlab(end)])
    set(gca,'FontSize',11)
    title(['Cumulative Effect of ',shocks{jj},' on Real Price of Oil'])
end
print( [ pwd '\figures\SP500_HD1' ] , '-dpdf' );

% Report FEVD for horizons 1,2,3,12
disp('FEVD for horizons 1, 2, 3, 12 is:');
disp([FEVD(4,:,2);FEVD(4,:,3);FEVD(4,:,4);FEVD(4,:,13)]);
%FEVD(4,:,501) % for horizon 500 (set hor=500 first)

% Graph FEVD
figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
bar(squeeze(FEVD(4,:,:))','stacked')
xlabel('Horizon')
ylim([0 100])
legend(shocks)
title('Decomposition of Real Oil Price Shocks to Variability in SP500 Real Returns')
set(gca,'FontSize',12)
print( [ pwd '\figures\SP500_FEVD' ] , '-dpdf' );

% Local Projection IRF
for ii = 1:n % shock
    if ii == 1 % qt is empty
        for  jj = 1:n % variable
            if jj == 1 % first variable is differenced
    
                [IRF_LP((ii-1)*n+jj,:),se((ii-1)*n+jj,:)] = Local_Projection_exo(100*log(rawdata(:,1)),data(:,ii),[],data,p,hor,'Differenced');
            
            elseif jj == 4 % last variable is differenced
    
                [IRF_LP((ii-1)*n+jj,:),se((ii-1)*n+jj,:)] = Local_Projection_exo(100*(log(SP500(:,1))-log(CPI(:,1))),data(:,ii),[],data,p,hor,'Differenced');
    
            else
                
                [IRF_LP((ii-1)*n+jj,:),se((ii-1)*n+jj,:)] = Local_Projection_exo(data(:,jj),data(:,ii),[],data,p,hor);   
            end
        end
    else % qt contains variables that come before shocked variable
        for  jj = 1:n % variable
            if jj == 1 % first variable is differenced
    
                [IRF_LP((ii-1)*n+jj,:),se((ii-1)*n+jj,:)] = Local_Projection_exo(100*log(rawdata(:,1)),data(:,ii),data(:,1:ii-1),data(:,ii:end),p,hor,'Differenced');
                % w contains lags of data(:,ii:end) instead of data(:,1:end) to prevent multicollinearity issues with qt
            elseif jj == 4 % last variable is differenced
    
                [IRF_LP((ii-1)*n+jj,:),se((ii-1)*n+jj,:)] = Local_Projection_exo(100*(log(SP500(:,1))-log(CPI(:,1))),data(:,ii),data(:,1:ii-1),data(:,ii:end),p,hor,'Differenced');
    
            else
                
                [IRF_LP((ii-1)*n+jj,:),se((ii-1)*n+jj,:)] = Local_Projection_exo(data(:,jj),data(:,ii),data(:,1:ii-1),data(:,ii:end),p,hor);
            end
        end
    end
end

% Normalize LP IRF with VAR IRF
IRF_LP(1:4,:) = IRF_LP(1:4,:)*IRF(1,1);
se(1:4,:) = se(1:4,:)*IRF(1,1);

IRF_LP(5:8,:) = IRF_LP(5:8,:)*IRF(6,1);
se(5:8,:) = se(5:8,:)*IRF(6,1);

IRF_LP(9:12,:) = IRF_LP(9:12,:)*IRF(11,1);
se(9:12,:) = se(9:12,:)*IRF(11,1);

IRF_LP(13:16,:) = IRF_LP(13:16,:)*IRF(16,1);
se(13:16,:) = se(13:16,:)*IRF(16,1);

% Expand y-axis to fit IRF
y_axis = [];
for jj = 1:n
    y_axis(jj,1) = 1.8*min(min([IRF(jj:n:end,:);IRF_LP(jj:n:end,:)]));
    y_axis(jj,2) = 1.2*max(max([IRF(jj:n:end,:);IRF_LP(jj:n:end,:)]));
end

% Graph LP IRF

%figure; fig = gcf; fig.PaperOrientation = 'landscape';
%set(fig,'PaperUnits','normalized');
%set(fig,'PaperPosition', [0 0 1 1]);
%for jj = 1:n
%    for ii = 1:n
%        subplot(n,n,(jj-1)*n+ii)
%        h1=plot(0:hor,IRF_LP((jj-1)*n+ii,:),'-k','LineWidth',3);
%        hold on
%        plot(0:hor,[IRF_LP((jj-1)*n+ii,:)+se((jj-1)*n+ii,:);IRF_LP((jj-1)*n+ii,:)-se((jj-1)*n+ii,:)],'--k','LineWidth',1.5)
%        hold on
%        h2=plot(0:hor,IRF((jj-1)*n+ii,:),'-r','LineWidth',2);
%        hold on
%        plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+ii,:,2:3)),'--r','LineWidth',1.5)
%        hold on
%        plot([0 hor],zeros(2,1),'-k','LineWidth',1)
%        ylim(y_axis(ii,:))
%        ylabel(variables{ii})
%        title(shocks{jj})
%        set(gca,'FontSize',11)
%        if jj == n
%            xlabel('Months')
%        end
%    end
%end
%legend([h1 h2],{'Local Projection','VAR'},'Orientation','horizontal')
%print( [ pwd '\figures\LP_IRF' ] , '-dpdf' );

figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
for jj = 1:n-1
    subplot(1,n-1,jj)
    h1=plot(0:hor,IRF_LP((jj-1)*n+n-1,:),'-k','LineWidth',3);
    hold on
    plot(0:hor,[IRF_LP((jj-1)*n+n-1,:)+se((jj-1)*n+n-1,:);IRF_LP((jj-1)*n+n-1,:)-se((jj-1)*n+n-1,:)],'--k','LineWidth',1.5)
    hold on
    h2=plot(0:hor,IRF((jj-1)*n+n-1,:),'-r','LineWidth',3);
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n-1,:,2:3)),'--r','LineWidth',1.5)
    hold on
    plot([0 hor],zeros(2,1),'-k','LineWidth',1)
    ylim(y_axis(n-1,:))
    ylabel(variables{n-1})
    title(shocks{jj})
    set(gca,'FontSize',11)
    xlabel('Months')
end
legend([h1 h2],{'Local Projection','VAR'},'Orientation','horizontal')
print( [ pwd '\figures\LP_IRF1' ] , '-dpdf' );

figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
for jj = 1:n-1
    subplot(1,n-1,jj)
    h1=plot(0:hor,IRF_LP((jj-1)*n+n,:),'-k','LineWidth',3);
    hold on
    plot(0:hor,[IRF_LP((jj-1)*n+n,:)+se((jj-1)*n+n,:);IRF_LP((jj-1)*n+n,:)-se((jj-1)*n+n,:)],'--k','LineWidth',1.5)
    hold on
    h2=plot(0:hor,IRF((jj-1)*n+n,:),'-r','LineWidth',3);
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n,:,2:3)),'--r','LineWidth',1.5)
    hold on
    plot([0 hor],zeros(2,1),'-k','LineWidth',1)
    ylim(y_axis(n,:))
    ylabel(variables{n})
    title(shocks{jj})
    set(gca,'FontSize',11)
    xlabel('Months')
end
legend([h1 h2],{'Local Projection','VAR'},'Orientation','horizontal')
print( [ pwd '\figures\LP_IRF2' ] , '-dpdf' );

%% Repeat without oil production and real activity

% Transform data
data = [100*log(100*(rawdata(2:end,2)./CPI(2:end,1))) 100*(diff(log(SP500(:,1)))-diff(log(CPI(:,1))))];
variables = {'Real Oil Price (%)', 'Cumulative Real Stock Returns (%)'};
shocks = {'Real Oil Price Shock','Other Shocks to Stock Returns'};

p = 24; % lags
hor = 15; % IRF and FEVD horizon
CI = [0.05 0.16 0.84 0.95]'; % to build IRF confidence interval
n = size(data,2);

% Estimate VAR(24)
[phi,sigma,x,e,phi_boot,sigma_boot] = olsvar(data,p,'double_bootstrap');

% IRF and FEVD
A0 = chol(sigma,'lower'); % cholesky identification
[IRF,FEVD] = calculate_IRF_FEVD(phi,A0,hor,'FEVD');
IRF(2:n:end,:) = cumsum(IRF(2:n:end,:),2); % cumulate because real stock returns is first differenced

% IRF confidence interval
IRF_boot = NaN(n^2,hor+1,1000);

for jj = 1:1000 % bootstrap after bootstrap
    A0_boot = chol(sigma_boot(:,:,jj),'lower');
    IRF_boot(:,:,jj) = calculate_IRF_FEVD(phi_boot(:,:,jj),A0_boot,hor);
    IRF_boot(2:n:end,:,jj) = cumsum(IRF_boot(2:n:end,:,jj),2);
end

IRF_quantiles = quantile(IRF_boot,CI,3);

% To input IRF scale
for jj = 1:n
    y_axis(jj,1) = 1.1*min(min(squeeze(IRF_quantiles(jj:n:end,:,1))));
    y_axis(jj,2) = 1.1*max(max(squeeze(IRF_quantiles(jj:n:end,:,4))));
end

% Graph IRF

%figure
%fig = gcf; fig.PaperOrientation = 'landscape';
%set(fig,'PaperUnits','normalized');
%set(fig,'PaperPosition', [0 0 1 1]);
%for jj = 1:n
%    for ii = 1:n
%        subplot(n,n,(jj-1)*n+ii)
%        plot(0:hor,IRF((jj-1)*n+ii,:),'-k','LineWidth',3)
%        hold on
%        plot([0 hor],zeros(2,1),'-k','LineWidth',1)
%        hold on
%        plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+ii,:,2:3)),'--b','LineWidth',1.5)
%        hold on
%        plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+ii,:,[1 4])),'--r','LineWidth',1.5)
%        ylim(y_axis(ii,:))
%        ylabel(variables{ii})
%        title(shocks{jj})
%        set(gca,'FontSize',11)
%        if jj == n
%            xlabel('Months')
%        end
%    end    
%end
%print( [ pwd '\figures\bivariate_IRF' ] , '-dpdf' );

figure
fig = gcf; fig.PaperOrientation = 'landscape';
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
for jj = 1:n-1
    subplot(1,n-1,jj)
    plot(0:hor,IRF((jj-1)*n+n,:),'-k','LineWidth',3)
    hold on
    plot([0 hor],zeros(2,1),'-k','LineWidth',1)
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n,:,2:3)),'--b','LineWidth',1.5)
    hold on
    plot(0:hor,squeeze(IRF_quantiles((jj-1)*n+n,:,[1 4])),'--r','LineWidth',1.5)
    ylim(y_axis(n,:))
    ylabel(variables{n})
    title(shocks{jj})
    set(gca,'FontSize',11)
    xlabel('Months')
end
print( [ pwd '\figures\bivariate_IRF1' ] , '-dpdf' );