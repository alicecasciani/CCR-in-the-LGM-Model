function stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365,indice3m,indice1y)
% Calcola il discount factor stocastico stepwise nel modello LGM (Bachelier), 
% utilizzando una volatilità piecewise costante definita tra le date in date_sigma.

% Determinazione dell'indice della volatilità da utilizzare a seconda dell'orizzonte temporale
if i <= indice3m
    indice = 1;                            % Fino a 3 mesi si usa la prima volatilità
elseif i > indice3m && i <= indice1y
    indice = 2;                            % Tra 3 mesi e 1 anno si usa la seconda
else
    indice = floor(i/indice1y) + 2;        % Oltre 1 anno si usa un'approssimazione più grossolana
end
if indice > 11
   indice = 11;                            % Cap a 11 per non andare fuori dal vettore sigma
end

int = 0;                                   % Inizializzazione dell'integrale

% Calcolo dell'integrale in intervalli multipli (fino all'intervallo corrente)
if indice > 1
    for ids = 1:indice-1
        % Calcolo degli year fraction rispetto a today
        yf_lunga = yearfrac(today, simulations_grid(i), ACT_365);
        yf_corta = yearfrac(today, simulations_grid(i-1), ACT_365);

        % Definizione dell'integranda relativa all'intervallo [date_sigma(ids), date_sigma(ids+1)]
        integranda = @(u) sigmaLGM_Bachelier(ids)^2 .* ((1 - exp(-a * (yf_lunga - u))) / a).^2 ...
                        - sigmaLGM_Bachelier(ids)^2 .* ((1 - exp(-a * (yf_corta - u))) / a).^2;

        % Integrazione su quell'intervallo
        int = int + quadgk(integranda, yearfrac(today, date_sigma(ids), ACT_365), yearfrac(today, date_sigma(ids+1), ACT_365));
    end

    % Ultimo intervallo incompleto, da date_sigma(indice) a t_{i-1}
    integranda = @(u) sigmaLGM_Bachelier(indice)^2 .* ((1 - exp(-a * (yf_lunga - u))) / a).^2 ...
                    - sigmaLGM_Bachelier(indice)^2 .* ((1 - exp(-a * (yf_corta - u))) / a).^2;
    int = int + quadgk(integranda, yearfrac(today, date_sigma(indice), ACT_365), yearfrac(today, simulations_grid(i-1), ACT_365));
else
    % Se siamo ancora nel primo intervallo di volatilità, si integra direttamente da 0 a t_{i-1}
    yf_lunga = yearfrac(today, simulations_grid(i), ACT_365);
    yf_corta = yearfrac(today, simulations_grid(i-1), ACT_365);
    integranda = @(u) sigmaLGM_Bachelier(1)^2 .* ((1 - exp(-a * (yf_lunga - u))) / a).^2 ...
                    - sigmaLGM_Bachelier(1)^2 .* ((1 - exp(-a * (yf_corta - u))) / a).^2;
    int = quadgk(integranda, 0, yearfrac(today, simulations_grid(i-1), ACT_365));
end

% Calcolo del coefficiente C che moltiplica il fattore stocastico
C = (1 - exp(-a * yearfrac(simulations_grid(i-1), simulations_grid(i), ACT_365))) / a;

% Parte deterministica A, include il fwd discount e il fattore di varianza
A = fwd_disc(i) * exp(-0.5 * int);

% Calcolo del discount factor stocastico
stoch_disc = A .* exp(-xtsLGM(:,i) .* C);
end
