% read in stock data and convert into matlab matrix for computation
matrix20102015 = table2array(readtable('output2010-2015.csv'));
matrix20152020 = table2array(readtable('output2015-2020.csv'));

% get return from raw ending prices 
return_matrix_20102015 = price2ret(matrix20102015);
return_matrix_20152020 = price2ret(matrix20152020);

% get average return from historical return
avg_return = mean(return_matrix_20102015,1);

% get the cov matrix from 2010 - 2015 data
cov20102015 = cov(return_matrix_20102015);

% get the RSI from starting from 2015.1.1 and break them into 30 days interval
RSI20152020 = table2array(readtable('stock_RSI_without_date.csv')); % 1298 lines of data here

% get 3-month risk-free from CSV table
rf = table2array(readtable('rf2517.csv'));

% get risk-free that we are going to use
rf20152020 = rf(1220:1220+1297);

% 30 * 40 = 1200 << 1297 total rows in RSI form 2015 - 2020
% In our simulation, we use 40 periods, each with 30 trading days of data.

% initializing variables
% days that have passed since first trading day of 2015
index = 1;
% rolling return starting from 2015.1.1 (starting value is 1)
overall_return = zeros(40,1);
% monthly return for each period
monthly_return = zeros(40,1);
% proportion of assest stored in bank for each period
bank_overall = zeros(40,1);
% the number of times the i-th stock is selected to construct a portfolio during 2015-2020
selected = zeros(95,1);
sharpe = zeros(40,1);
all_risky = zeros(40,95);

for iter = 1:40
    disp(iter)
    first_ind = index;     % starting day
    last_ind = index + 30;      % ending day
    index = index + 30;

    temp_return = avg_return;   % avg return array that we will use for markowiz
    temp_cov = cov20102015;     % covariance matrix that we will use for markowiz
    long = zeros(95,1);         % stocks that we will put into portfolio depending on RSI
    short = zeros(95,1);        % stocks that we will put into portfolio depending on RSI
    count = 0;  % the number of stocks chosen from the base 95 ones to construct portfolio
    
    
%     for l = 1:95
%         if long(l) == 0 && short(l)==0
%             temp_return(:,l) = 0;
%             temp_cov(:,l) = 0;
%             temp_cov(l,:) = 0;
%         end 
%     end 
    
    mu0_3month = rf20152020(first_ind)/100;     % 3-month rf rate
    mu0 = (1 + mu0_3month)^(1/66) - 1;          % daily rf rate: converting 3-month risk-free to daily risk-free
    mu = temp_return';                          % risky assets: average historical daily returns for
    V = temp_cov;                               % risky assets: covariance matrix
    sigma = .05;                                % a user chosen upper bound on portfolio risk
    xx0 = 1/96;                                 % initial bank proportion: will be optimized
    xx = 95/96;                                 % initial risky assest proportion: will be optimized
    trans_cost = 0;                             % the transaction cost per wealth unit bought or sold
    
    [x0,x] = cvx_markowitz1_2(mu0,mu,V,sigma,xx0,xx,trans_cost);%,long,short);   % x0 is bank, x is risky
    bank = x0;                                  % current period bank proportion
    risky_assets = x;                          % current period risky proportion
    all_risky(iter,:) = risky_assets';
    % current period return for the entire month
    curr_mu = (matrix20152020(last_ind,:) - matrix20152020(first_ind,:))./matrix20152020(first_ind,:);
    
    % monthly return
    % monthly = mu0_3month * x0 + curr_mu * x;
    monthly = mu0_3month * x0 + curr_mu * x;
    % annualizing the monthly return
    annualized = (1+monthly)^12-1;
    
    % calculate sharpe ratio for the current period
    
    % testing iter=32 because it's very large
%     if iter==32
%         % disp("mu "+mu)
%         % disp("xx: "+xx)
%         % disp("mu0: "+ mu0)
%         disp("chosen "+find(long))
%         
%         disp("annualize "+annualized)
%         disp("monthly "+monthly)
%         
%         disp("bank: "+ bank)
%         disp("risky_assets: "+ round(risky_assets,3))
%         disp("curr_mu': "+ curr_mu')
%     end
    
    overall_return(iter) = annualized;
    monthly_return(iter) = monthly;
    bank_overall(iter) = bank;
    sharpe(iter) = (monthly - mu0_3month) / (risky_assets'* temp_cov * risky_assets);
    
end

figure
plot(monthly_return)
title('monthly returns for i-th period')
figure()
plot(cumprod(monthly_return+1))
title('cumulative returns')
% figure
% plot(bank_overall)
% title('proportion in bank for i-th period')
% figure
% plot(sharpe)
% title('sharpe ratio for i-th period')
% figure
% bar(selected)
% title('the number of times the i-th stock is selected')

temp = (1 + mean(rf20152020)/100)^(250/66) - 1;
cumulative_return = cumprod(monthly_return+1);
annual_return = cumulative_return(end:end)^(1/5)-1;
overall_sharpe = (mean(annual_return) - temp) / (std(monthly_return)*sqrt(12))

