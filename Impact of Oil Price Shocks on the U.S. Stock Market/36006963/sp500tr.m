clear all
clc

addpath('functions')
addpath('figures')
addpath('data')

% Load data
rawdata = xlsread('data.xlsx','oil_market_data','B177:D398');
CPI = xlsread('data.xlsx','US_monthly_data','C493:C714');
SP500TR = xlsread('sp500total.xlsx','sp500','B2:B223');
yearlab =(1988+8/12:1/12:2006+12/12)'; % timeline

% Transform data
data = [100*diff(log(rawdata(:,1))) rawdata(2:end,3) 100*log(100*(rawdata(2:end,2)./CPI(2:end,1))) 100*(diff(log(SP500TR(:,1)))-diff(log(CPI(:,1))))];
data(:,1) = -1*data(:,1); % normalize to represent negative supply shock
variables = {'Oil Production (%)','Real Activity','Real Oil Price (%)', 'Cumulative Real Stock Returns (%)'};
shocks = {'Oil Supply Shock','Aggregate Demand Shock','Oil-Specific Demand Shock','Other Shocks to Stock Returns'};

p = 18; % lags
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
%print( [ pwd '\figures\extended_IRF' ] , '-dpdf' );

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
print( [ pwd '\figures\SP500TR_IRF1' ] , '-dpdf' );

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
print( [ pwd '\figures\SP500TR_IRF2' ] , '-dpdf' );

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
%print( [ pwd '\figures\extended_HD' ] , '-dpdf' );

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
print( [ pwd '\figures\SP500TR_HD1' ] , '-dpdf' );

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
set(gca,'FontSize',14)
print( [ pwd '\figures\SP500TR_FEVD' ] , '-dpdf' );