function [collateral_Delta, collateral_Epsilon] = Weekly_LGM_Collateral(Portfolio, jamshResults, dates, discounts, simulations_num, aLGM, sigmaBach)
% Imposta il seed per la generazione casuale per garantire riproducibilità
rng(42);

% Prende la data di partenza dalla prima data fornita
today=dates(1);

% Costanti per il calcolo degli year fractions
EU_30_360=6;
ACT_365=3; 

% Parametro di mean reversion del modello LGM
a=aLGM;

% Definisce una griglia temporale iniziale con oggi e oggi + 3 mesi
time_grid = [today, today + calmonths(3)];

% Aggiunge alla griglia temporale le date di business day offset di oggi + anni da 1 a 10
time_grid = businessdayoffset([time_grid, today + calyears(1):calyears(1):today + calyears(10)]);

% Estrae il tasso fisso dello swap dal portafoglio
swapRate = Portfolio.FixedLegRate;

% Estrae la maturità dello swap
maturity = Portfolio.Maturity; 

% Estrae la tipologia di fixed leg (Receive o Pay)
fixedLeg = Portfolio.FixedLeg;

% Inizializza un vettore di zeri per classificare il tipo di swap
swapType = zeros(size(fixedLeg));             

% Imposta +1 per Receive
swapType(strcmpi(fixedLeg,'Receive')) =  1;   

% Imposta -1 per Pay
swapType(strcmpi(fixedLeg,'Pay'))     = -1;   

% Assegna le date per la volatilità al time_grid
date_sigma=time_grid;

% Calcola il notional moltiplicato per 1 milione (scala grandezza)
notional = 10^6*Portfolio.NotionalAmount; 

% Prende la volatilità Bachelier calcolata da jamshResults
sigmaLGM_Bachelier = jamshResults.BACHELIER;

% Genera le date dei pagamenti della fixed leg a intervalli settimanali fino a 10 anni
fixed_leg_payments_dates = businessdayoffset(today:calweeks(1):today+calyears(10));

% Calcola i fattori di sconto interpolati per le date dei pagamenti
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);

% Calcola i delta temporali (year fraction) fra le date dei pagamenti secondo convenzione EU_30_360
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

% Definisce la griglia delle simulazioni (escludendo la prima data di pagamento)
simulations_grid = fixed_leg_payments_dates(2:end);

% Numero di step di simulazione
N=length(simulations_grid);

% Prende il nome della controparte dal portafoglio
counterparty = Portfolio.CounterpartyName;

% Conta quante controparti sono di tipo 'Epsilon'
nEpsilon = sum(strcmp(counterparty,'Epsilon'));   

% Conta quante controparti sono di tipo 'Delta'
nDelta = sum(strcmp(counterparty,'Delta'));

% Calcola il flusso di cassa fisso moltiplicando il tasso per i deltas
fixed_cash_flow = swapRate.*deltas;

% Genera una matrice di numeri casuali standard normali per le simulazioni
z=randn(simulations_num, N);

% Definisce le dimensioni per le matrici di output
r = simulations_num; c = N;      

% Inizializza matrici 3D per simulazioni LGM separate per controparti Epsilon e Delta
TLGM_Epsilon = zeros(r, c, nEpsilon);
TLGM_Delta = zeros(r, c, nDelta);

% Inizializza la matrice degli stati x_t per tutte le simulazioni e tutti i tempi
xtsLGM = zeros(simulations_num, N);

% Contatori per la gestione delle controparti Delta ed Epsilon
countDelta = 1; 
countEpsilon = 1;

% Loop sulle date della griglia di simulazione
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx); % data di simulazione corrente
    
    if idx == 1
        prev_sim_date = today; % alla prima iterazione la data precedente è today
        x_t           = zeros(1, simulations_num); % stato iniziale LGM a zero
    else
        prev_sim_date = simulations_grid(idx-1); % data di simulazione precedente
    end
    
    % Simula lo stato x_t al tempo corrente con volatilità variabile a seconda dell'indice
    if idx < 12
        % primi 11 step con volatilità sigmaLGM_Bachelier(1)
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(1), x_t, prev_sim_date, sim_date);
    elseif idx >= 12 && idx < 52
        % da 12 a 51 step con volatilità sigmaLGM_Bachelier(2)
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(2), x_t, prev_sim_date, sim_date);
    elseif idx < N
        % dopo il primo anno si prende la volatilità dalla posizione indice, max 11
        indice = floor(idx/52) + 2;
        if indice>11
            indice=11;
        end
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(indice), x_t, prev_sim_date, sim_date);
    elseif idx == N
        % all'ultimo step si usa sigmaBach passato come parametro
        x_t = simulateLGM(z(:, idx), aLGM, sigmaBach, x_t, prev_sim_date, sim_date);
    end
    
    % Salva il valore simulato x_t per tutte le simulazioni alla data idx
    xtsLGM(:, idx) = x_t;
end

for j=1:length(maturity)
    
    % Inizializza le matrici MtM (Mark-to-Market) per controparti Epsilon e Delta per ogni simulazione e step temporale
    MtMs_Epsilon = zeros(simulations_num, N);
    MtMs_Delta = zeros(simulations_num, N);

    % Costruisce la griglia temporale delle simulazioni settimanali in base alla maturità corrente
    if maturity >= 1
        % Se maturità è almeno 1 anno, crea griglia con step settimanali fino a 'maturity(j)' anni
        simulations_grid = businessdayoffset(today:calweeks(1):today+calyears(maturity(j)));
    else
        % Se maturità è meno di 1 anno, crea griglia con step settimanali fino a 'maturity(j)' mesi
        simulations_grid = businessdayoffset(today:calweeks(1):today+calmonths(maturity(j)*12));
    end

    % Rimuove la prima data (oggi) dalla griglia delle simulazioni
    simulations_grid = simulations_grid(2:end);
    Nstep = numel(simulations_grid);

    % Prima data di simulazione
    sim_date = simulations_grid(1);

    % Calcola i coefficienti A e C tramite la funzione affine_trick_LGM,
    % che dipendono dalla data di simulazione e da parametri LGM, curve di sconto e griglia temporale
    [A, C] = affine_trick_LGM(sim_date, simulations_grid, aLGM, sigmaLGM_Bachelier, discounts_timegrid, date_sigma, dates, discounts);

    for idx = 1:length(simulations_grid)
        
        % Data corrente della simulazione
        sim_date = simulations_grid(idx);
        % Stato del processo LGM (x_t) per tutti gli scenari, trasposto a vettore riga
        x_t = xtsLGM(:, idx)';
        
        % Caso maturità >= 1 anno e data simulazione prima della maturità
        if (maturity(j) >= 1) && sim_date < businessdayoffset(today + calyears(maturity(j)))

            % Se la controparte è Epsilon
            if counterparty(j) == "Epsilon"
                % Calcola i flussi attualizzati usando la formula affine, per tutti gli scenari
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);  
                % Somma dei flussi anticipati
                sum_term = sum(term1, 1);
                % Termine terminale, rappresenta il valore finale della componente esponenziale
                terminal = A(end) * exp(-C(end) .* x_t);
                % Correzione per gli step di simulazione intermedi (adjust=1 se idx < Nstep)
                adjust = (idx < Nstep);
                % Calcolo MtM moltiplicando per notional e tipo di swap (pay/receive)
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                % Assegna MtM alla matrice di controparti Epsilon
                MtMs_Epsilon(:, idx) = MtM;
                % Shift dei vettori A e C per il prossimo step
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];

            % Se la controparte è Delta
            elseif counterparty(j) == "Delta"
                % Calcoli analoghi per la controparte Delta
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
            end

        % Caso maturità inferiore a 1 anno e sim_date prima della maturità in mesi
        elseif sim_date < businessdayoffset(today + calmonths(maturity(j)*12))

            % Logica identica al blocco precedente, ma con scadenza in mesi
            if counterparty(j) == "Epsilon"
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                MtMs_Epsilon(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];

            elseif counterparty(j) == "Delta"
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
            end
        end

    end

    % Assegna i risultati MtM calcolati ai tensori TLGM in base alla controparte
    if counterparty(j) == "Epsilon"
        TLGM_Epsilon(:,:,countEpsilon) = MtMs_Epsilon;
        countEpsilon = countEpsilon + 1;
    else
        TLGM_Delta(:,:,countDelta) = MtMs_Delta;
        countDelta = countDelta + 1;
    end

end

% Somma lungo la terza dimensione (numero contratti) per ottenere MtM aggregati per ogni controparte
MtM_Epsilon = sum(TLGM_Epsilon, 3);
MtM_Delta = sum(TLGM_Delta, 3);

% Inizializza le matrici di collateral (collateralizzazione)
collateral_Epsilon = zeros(simulations_num, N);
collateral_Delta = zeros(simulations_num, N);

% Copie di MtM per calcolo collateralizzato
MtM_coll_Epsilon = MtM_Epsilon;
MtM_coll_Delta = MtM_Delta;

% Griglia temporale dei pagamenti fissi (escludendo la prima data)
simulations_grid = fixed_leg_payments_dates(2:end);

% Calcola i fattori di sconto tra step nella griglia delle simulazioni
disc_simgrid = interpolation(discounts, dates, simulations_grid);
fwd_disc = disc_simgrid(2:end) ./ disc_simgrid(1:end-1);

% Indici di riferimento per trimestri e anni (12 settimane = 3 mesi, 52 settimane = 1 anno)
indice3m = 12;
indice1y = 52;

% Loop per aggiornare le esposizioni collateralizzate della controparte Epsilon
for i=1:N-1
    if i == 1
        % Primo step: non si aggiorna collateral, si lascia MtM invariato
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i);
    else
        % Calcola fattore di sconto stocastico per il passo i
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        % Aggiunge l’effetto del collateral precedente attualizzato
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i-1) ./ stoch_disc;
        % Aggiorna collateral come negativo dell’esposizione MtM collateralizzata
        collateral_Epsilon(:, i) = -MtM_coll_Epsilon(:, i);
        % Resetta esposizione corrente a zero per evitare doppio conteggio
        MtM_coll_Epsilon(:, i) = 0;
    end
end

% Stessa logica per controparti Delta
for i=1:N-1   
    if i == 1
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i);
    else
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i-1) ./ stoch_disc;
        collateral_Delta(:, i) = -MtM_coll_Delta(:, i);
        MtM_coll_Delta(:, i) = 0;
    end
end
