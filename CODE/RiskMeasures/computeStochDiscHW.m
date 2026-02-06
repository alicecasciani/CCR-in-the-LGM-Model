function stoch_disc=computeStochDiscHW(today,simulations_grid, i, ACT_365,a,sigmaBach,xts,fwd_disc)
% Calcola il discount factor stocastico nel modello Hull-White (versione Bachelier),
% tra i tempi simulations_grid(i-1) e simulations_grid(i), in base al valore del processo xts

% Calcolo degli year fraction dal giorno odierno alle date di griglia i e i-1
yf_lunga=yearfrac(today,simulations_grid(i),ACT_365);          % year fraction fino a t_i
yf_corta=yearfrac(today,simulations_grid(i-1),ACT_365);        % year fraction fino a t_{i-1}

% Definizione della funzione integranda per il calcolo della parte deterministica della volatilità
% (differenza tra i due contributi integrati nel tempo per ottenere la varianza)
integranda = @(u) sigmaBach^2.*((1-exp(-a*(yf_lunga-u)))./a).^2 ...
                   - sigmaBach^2.*((1-exp(-a*(yf_corta-u)))./a).^2;

% Integrazione numerica dell'integranda tra 0 e t_{i-1} per ottenere la varianza incrementale
int=quadgk(integranda, 0, yearfrac(today ,simulations_grid(i-1),ACT_365));

% Calcolo del coefficiente C che moltiplica il valore del processo xts nel discount factor
C=(1-exp(-a*yearfrac(simulations_grid(i-1),simulations_grid(i),ACT_365)))/a;

% Parte deterministica A del discount factor, che include il fwd discount e il termine di varianza
A=fwd_disc(i)*exp(-0.5*int);

% Calcolo finale del discount factor stocastico: parte deterministica * parte stocastica
stoch_disc = A.*exp(-xts(:,i).*C);
