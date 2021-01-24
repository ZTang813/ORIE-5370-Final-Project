% 0. regular RSI            
% 1. quantile trading: short 0.95 percentile RSI stocks, long 0.05
% percentile RSI stocks
% 2. differenced RSI (quantiled) trading: take difference of RSIs between this
% month and last month, short >.95 percentile change stocks, long <.05
% percentile RSI stocks
% 3. quantile + RSI


% ------------------ Below is the code for Quantile Trading stradegy ---------------------

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

BAND_WIDTH = [1 3 5 7];
[m,n] = size(BAND_WIDTH);
MONTHLY_ALL = zeros(40,n);
CUM_ALL = zeros(40,n);
SHARPE_ALL = zeros(1,n);
% RISKY_ALL = zeros(n,40,95);

for sim = 1:n 
    var = BAND_WIDTH(sim);
    disp(var)
    index = 1;
    
    monthly_return = zeros(40,1);                           % monthly return for each period
    bank_overall = zeros(40,1);                             % proportion of assest stored in bank for each period
    selected = zeros(95,1);                                 % the number of times the i-th stock is selected to construct a portfolio during 2015-2020
    all_risky = zeros(40,95);                               % the weight of each asset at time period i

    for iter = 1:40
        % disp(iter)
        first_ind = index;                                  % starting day
        last_ind = index + 30;                              % ending day
        index = index + 30;                
        selected = zeros(95,1);  
        temp_return = avg_return;                           % avg return array that we will use for markowiz
        temp_cov = cov20102015;                             % covariance matrix that we will use for markowiz
        long = zeros(95,1);                                 % stocks that we will put into portfolio depending on RSI
        short = zeros(95,1);                                % stocks that we will put into portfolio depending on RSI
        count = 0;                                          % the number of stocks chosen from the base 95 ones to construct portfolio

        % if RSI on the starting date is between 35 and 65, we discard it
        for stock = 1:95                                    % 95 stocks in total
            % choose and assign long or short
            upper = prctile(RSI20152020(first_ind,:),100-var);
            lower = prctile(RSI20152020(first_ind,:),var);

            if RSI20152020(first_ind,stock)>=upper
                short(stock) = 1;
                count=count+1;
                selected(stock) = selected(stock)+1;
            end

            if RSI20152020(first_ind,stock) <= lower
                long(stock) = 1;
                count=count+1;
                selected(stock) = selected(stock)+1;
            end
        end

        if count==0
            % disp("iter " + iter + " has no valid stocks")
            continue
        end
        long_short = long+short;
        long(find(1-long_short)) = [];
        short(find(1-long_short)) = [];
        
        temp_return(find(1-selected)) = [];
        temp_cov(find(1-selected),:) = [];
        temp_cov(:,find(1-selected)) = [];
        
%         disp(temp_cov)
%         disp(temp_return)
        
        % now we have the stocks that we want to long, short or do nothing and
        % the stock's information is stored in the variable 'long'

        % disp("long "+find(long))
        % disp("short "+find(short))

        % -------------------------- Begin Markowitz Optimization ------------------------
        mu0_3month = rf20152020(first_ind)/100;         % 3-month rf rate
        mu0 = (1 + mu0_3month)^(1/66) - 1;              % daily rf rate: converting 3-month risk-free to daily risk-free
        mu = temp_return';                              % risky assets: average historical daily returns for
        V = temp_cov;                                   % risky assets: covariance matrix
        sigma = .05;                                    % a user chosen upper bound on portfolio risk
        xx0 = 1/(count+1);                              % initial bank proportion: will be optimized
        xx = 1/(count+1) * ones(count,1);               % initial risky assest proportion: will be optimized
        trans_cost = 0; %0.0005                         % the transaction cost per wealth unit bought or sold
        
        [x0,x] = cvx_markowitz1_2(mu0,mu,V,sigma,xx0,xx,trans_cost,long,short);   % x0 is bank, x is risky
        bank = x0;                                      % current period bank proportion
        risky_assets = x;                               % current period risky proportion
        % disp(risky_assets)
        
        % current period return for the entire month
        curr_mu = (matrix20152020(last_ind,:) - matrix20152020(first_ind,:))./matrix20152020(first_ind,:);   
        curr_mu(find(1-selected)) = [];
        
        monthly = mu0_3month * x0 + curr_mu * x;        % monthly return
        
        annualized = (1+monthly)^12-1;                  % annualizing the monthly return (not used)
        monthly_return(iter) = monthly;                 % saving data in global variable
        bank_overall(iter) = bank;                      % saving data in global variable
        % RISKY_ALL = zeros(n,40,95);
        
        
    end

    temp = (1 + mean(rf20152020)/100)^(250/66) - 1;                                     
    cumulative_return = cumprod(monthly_return+1);                                      % calculating cumulative return
    annual_return = cumulative_return(end:end)^(1/5)-1;                                 % calculating annual return (not used)
    overall_sharpe = (mean(annual_return) - temp) / (std(monthly_return)*sqrt(12));     % calculating overall sharpe ratio for this simulation
    
    MONTHLY_ALL(:,sim) = monthly_return;
    CUM_ALL(:,sim) = cumulative_return;
    SHARPE_ALL(sim) = overall_sharpe;
    
end

% -------------------------- Plotting monthly returns ----------------------------
figure
for i = 1:n
   plot(MONTHLY_ALL(:,i),'LineWidth',1.5,'DisplayName','Quantile <'+string(BAND_WIDTH(i)) + '% and >' + string(100-BAND_WIDTH(i))+"%");
   hold on
end
yline(0,'-.r','DisplayName','zero');
hold on
xlabel('period')
ylabel('return')
title('Quantile Trading - monthly returns for i-th period')
hold off
legend show

% -------------------------- Plotting cumulative returns -----------------------------
figure()
for i = 1:n
   plot(CUM_ALL(:,i),'LineWidth',1.5,'DisplayName','Quantile <'+string(BAND_WIDTH(i)) + '% and >' + string(100-BAND_WIDTH(i))+"%");
   hold on
end
yline(1,'-.r','DisplayName','one');
hold on
title('Quantile Trading - cumulative returns until i-th period')
xlabel('period')
ylabel('cumulative return')
hold off
legend show


