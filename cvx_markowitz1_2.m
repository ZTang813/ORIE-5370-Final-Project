% In addition to the regular markowitz optimization, we introduce two new
% variables and restictions long and short to help us determine whether or
% not to long or short the corresponding stock in the portfolio

function [x0,x] =  cvx_markowitz1_2(mu0,mu,V,sigma,xx0,xx,trans_cost,long,short);
    
    n = length(mu);
    cvx_begin quiet
    
        variables x0 x(n) y(n) total_trans_cost
        
        maximize((mu0+1) * x0 + x' * (mu+1)); % maximizing wealth
        
        subject to
                x'*V*x <= sigma^2;
                x0 + sum(x) + total_trans_cost == 1;
                x == xx + y;
                trans_cost * sum(abs(y)) <= total_trans_cost;
                
                % restrictions:
                x0 >= -1;
                x0 <= 1;
                abs(x) <= 1;
                
                x(find(long)) >=0;
                x(find(short))<=0;
                x(find(1-(long+short))) == 0;
    cvx_end
    
    (mu0+1) * x0 + x' * (mu+1);
end

