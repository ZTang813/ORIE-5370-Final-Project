% 0. regular RSI            -> GJL write about intuition about why percentage is low
% 1. quantile trading: short 0.95 percentile RSI stocks, long 0.05
% percentile RSI stocks
% 2. differenced RSI (quantiled) trading: take difference of RSIs between this
% month and last month, short >.95 percentile change stocks, long <.05
% percentile RSI stocks
% 3. quantile + RSI


% --------- Below is the code for Differenced RSI + Quantile Trading stradegy ------------

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

BAND_WIDTH = [10 20 30];
RSI_WIDTH = [20 30 40];
[m1,n1] = size(BAND_WIDTH);
[m2,n2] = size(RSI_WIDTH);

MONTHLY_ALL = zeros(40,n1*n2);
CUM_ALL = zeros(40,n1*n2);
SHARPE_ALL = zeros(1,n1*n2);

for sim1 = 1:n1
for sim2 = 1:n2
    band = BAND_WIDTH(sim1);
    rsi = RSI_WIDTH(sim2);
    disp(band)
    disp(rsi)
    index = 1;
    
    monthly_return = zeros(40,1);                       % monthly return for each period
    bank_overall = zeros(40,1);                         % proportion of assest stored in bank for each period
    all_risky = zeros(40,95);                           % the weight of each asset at time period i
    
    for iter = 1:40
        % disp(iter)
        first_ind = index;                              % starting day
        last_ind = index + 30;                          % ending day
        index = index + 30;
        if iter == 1 
            continue                                    % because we start from 2015.01.01
        end
        
        selected = zeros(95,1);                         % whether the stock is selected in the current period
        temp_return = avg_return;                       % avg return array that we will use for markowiz
        temp_cov = cov20102015;                         % covariance matrix that we will use for markowiz
        long = zeros(95,1);                             % stocks that we will put into portfolio depending on RSI
        short = zeros(95,1);                            % stocks that we will put into portfolio depending on RSI
        count = 0;                                      % the number of stocks chosen from the base 95 ones to construct portfolio

        % if RSI on the starting date is between 35 and 65, we discard it
        for stock = 1:95            % 95 stocks in total
            % choose and assign long or short

            diff = RSI20152020(first_ind,:) - RSI20152020(first_ind - 30,:);

            upper = prctile(diff,100-band);
            lower = prctile(diff,band);

            upperRSI = prctile(RSI20152020(first_ind,:),100-rsi);
            lowerRSI = prctile(RSI20152020(first_ind,:),rsi);

            if diff(stock)>=upper && RSI20152020(first_ind,stock)<=lowerRSI
                long(stock) = 1;
                count=count+1;
                selected(stock) = selected(stock)+1;
            end
            if diff(stock) <= lower && RSI20152020(first_ind,stock)>=upperRSI
                short(stock) = 1;
                count=count+1;
                selected(stock) = selected(stock)+1;
            end
        end
        
        if count==0                                         % debugging step
            % disp("iter " + iter + " has no valid stocks")
            continue
        end
        
        % now we have the stocks that we want to long, short or do nothing and
        % the stock's information is stored in the variable 'long'
        
        % disp("long "+find(long))
        % disp("short "+find(short))
        
        % Generating valid covariance matrixes and return arrays for Markowitz optimization
        long_short = long+short;
        long(find(1-long_short)) = [];
        short(find(1-long_short)) = [];
        
        temp_return(find(1-selected)) = [];
        temp_cov(find(1-selected),:) = [];
        temp_cov(:,find(1-selected)) = [];   
        
        % -------------------------- Begin Markowitz Optimization ------------------------
        mu0_3month = rf20152020(first_ind)/100;     % 3-month rf rate
        mu0 = (1 + mu0_3month)^(1/66) - 1;          % daily rf rate: converting 3-month risk-free to daily risk-free
        mu = temp_return';                          % risky assets: average historical daily returns for
        V = temp_cov;                               % risky assets: covariance matrix
        sigma = .05;                                % a user chosen upper bound on portfolio risk
        xx0 = 1/(count+1);                          % initial bank proportion: will be optimized
        xx = 1/(count+1) * ones(count,1);           % initial risky assest proportion: will be optimized
        trans_cost = 0;                             % the transaction cost per wealth unit bought or sold
        
        [x0,x] = cvx_markowitz1_2(mu0,mu,V,sigma,xx0,xx,trans_cost,long,short);   % x0 is bank, x is risky
        bank = x0;                                  % current period bank proportion
        risky_assets = x;                           % current period risky proportion
         
        % current period return for the entire month
        curr_mu = (matrix20152020(last_ind,:) - matrix20152020(first_ind,:))./matrix20152020(first_ind,:);   
        curr_mu(find(1-selected)) = [];             % remove the un-selected stocks from the initial pool of 95 stocks

        monthly = mu0_3month * x0 + curr_mu * x;    % monthly return
        annualized = (1+monthly)^12-1;              % annualizing the monthly return (not used)
        monthly_return(iter) = monthly;             % saving data in global variable
        bank_overall(iter) = bank;                  % saving data in global variable

    end
%     disp(monthly);
    temp = (1 + mean(rf20152020)/100)^(250/66) - 1;                                     
    cumulative_return = cumprod(monthly_return+1);                                      % calculating cumulative return
    annual_return = cumulative_return(end:end)^(1/5)-1;                                 % calculating annual return (not used)
    overall_sharpe = (mean(annual_return) - temp) / (std(monthly_return)*sqrt(12));     % calculating overall sharpe ratio for this simulation
    %disp((sim1-1)*n2 + sim2)
    
    MONTHLY_ALL(:,(sim1-1)*n2 + sim2) = monthly_return;
    CUM_ALL(:,(sim1-1)*n2 + sim2) = cumulative_return;
    SHARPE_ALL((sim1-1)*n2 + sim2) = overall_sharpe;
    
end
end

% -------------------------- Plotting monthly returns ----------------------------
figure
for i = 1:n1
for j = 1:n2
   plot(MONTHLY_ALL(:,(i-1)*n2 + j),'LineWidth',1.5,'DisplayName','Quantile '+string(BAND_WIDTH(i)) + '% and RSI% ' + string(RSI_WIDTH(j))+'%');
   hold on
end
end
yline(0,'-.r','DisplayName','zero');
hold on
xlabel('period')
ylabel('return')
title('Differenced RSI + Quantile Trading - monthly returns for i-th period')
hold off
legend show

% -------------------------- Plotting cumulative returns -----------------------------
figure()
for i = 1:n1
for j = 1:n2
   plot(CUM_ALL(:,(i-1)*n2 + j),'LineWidth',1.5,'DisplayName','Quantile '+string(BAND_WIDTH(i)) + '% and RSI% ' + string(RSI_WIDTH(j))+'%');
   hold on
end
end
yline(1,'-.r','DisplayName','one');
hold on
title('Differenced RSI + Quantile Trading - cumulative returns until i-th period')
xlabel('period')
ylabel('cumulative return')
hold off
legend show


