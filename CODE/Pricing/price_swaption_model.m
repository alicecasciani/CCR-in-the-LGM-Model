function price = price_swaption_model(dates, discounts, dates_yearly, today, mkt_vol, strike, model)
% price_swaption_model Computes swaption price under specified model
%   price = price_swaption_model(dates, discounts, dates_yearly, today,
%           mkt_vol, strike, model)
%
%   models: 'BLACK' or 'BACHELIER'
%
%   INPUTS and OUTPUT same as underlying pricing functions

    switch upper(model)
        case 'BLACK'
            mkt_vol = mkt_vol / 100;
            price = BlackPrice(dates, discounts, dates_yearly, today, mkt_vol, strike);

        case 'BACHELIER'
            mkt_vol = mkt_vol / 10000;
            price = BachelierPrice(dates, discounts, dates_yearly, today, mkt_vol, strike);

        otherwise
            error('Modello non riconosciuto. Usa "BLACK" o "BACHELIER".');
    end
end
