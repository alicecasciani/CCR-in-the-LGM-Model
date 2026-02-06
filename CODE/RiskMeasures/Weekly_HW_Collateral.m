function [MtM_coll_Delta, MtM_coll_Epsilon, collateral_Delta, collateral_Epsilon] = Weekly_HW_Collateral(Portfolio, dates, discounts, simulations_num, aBach, sigmaBach)
% Imposta il seed del generatore di numeri casuali per riproducibilità
rng(42);

% Definisce la data di oggi come il primo elemento dell'array 'dates'
today = dates(1);

% Definizioni convenzioni day count
EU_30_360 = 6;
ACT_365 = 3;

% Parametro di mean reversion per il modello LGM
a = aBach;

% Estrazione dati dal portafoglio
swapRate = Portfolio.FixedLegRate;      % Tasso fisso dello swap
maturity = Portfolio.Maturity;           % Scadenza dello swap
fixedLeg = Portfolio.FixedLeg;           % Tipo di flusso fisso (Pay o Receive)

% Codifica del tipo di swap: +1 per Receive, -1 per Pay, 0 altrimenti
swapType = zeros(size(fixedLeg));             
swapType(strcmpi(fixedLeg,'Receive')) =  1;   
swapType(strcmpi(fixedLeg,'Pay'))     = -1;   

% Calcolo del nozionale totale (scala 10^6)
notional = 10^6 * Portfolio.NotionalAmount; 

% Costruzione del calendario dei pagamenti fissi con step settimanale
fixed_leg_payments_dates = businessdayoffset(today:calweeks(1):today+calyears(10));

% Interpolazione dei fattori di sconto per le date di pagamento fisse
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);

% Calcolo delle frazioni di anno fra i pagamenti, secondo convenzione EU_30_360
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

% Definisce la griglia temporale delle simulazioni (date di pagamento escluse la prima)
simulations_grid = fixed_leg_payments_dates(2:end);
N = length(simulations_grid);   % Numero di step temporali

% Estrazione nome delle controparti e conteggio di quante sono Delta ed Epsilon
counterparty = Portfolio.CounterpartyName;
nEpsilon = sum(strcmp(counterparty, 'Epsilon'));
nDelta = sum(strcmp(counterparty, 'Delta'));

% Calcolo dei flussi fissi moltiplicando il tasso per i delta temporali
fixed_cash_flow = swapRate .* deltas;

% Generazione matrice di variabili gaussiane standard per la simulazione Monte Carlo
z = randn(simulations_num, N);

% Inizializzazione matrici per risultati LGM per entrambe le controparti
r = simulations_num; c = N;      % Dimensioni di ciascuna matrice di simulazione
TLGM_Epsilon = zeros(r, c, nEpsilon);
TLGM_Delta = zeros(r, c, nDelta);

% Matrice per salvare i valori simulati x_t del modello Hull-White (LGM)
xtsHW = zeros(simulations_num, N);

% Contatori per le controparti Delta ed Epsilon (per eventuale popolamento dati)
countDelta = 1; 
countEpsilon = 1;

% Loop temporale sulle date di simulazione
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    if idx == 1
        % Per la prima data di simulazione, il tempo precedente è "today"
        prev_sim_date = today;
        % Valore iniziale della simulazione LGM a zero (condizione iniziale)
        x_t = zeros(1, simulations_num);
    else
        % Altrimenti il tempo precedente è la data di simulazione precedente
        prev_sim_date = simulations_grid(idx - 1);
    end
    % Simulazione del valore LGM al tempo sim_date usando il passo di discretizzazione
    x_t = simulateLGM(z(:, idx), aBach, sigmaBach, x_t, prev_sim_date, sim_date);
    % Salvataggio dei valori simulati nella matrice xtsHW
    xtsHW(:, idx) = x_t;
end
for j = 1:length(maturity)
    % Inizializza le matrici MtM per le due controparti ad ogni iterazione
    MtMs_Epsilon = zeros(simulations_num, N);
    MtMs_Delta = zeros(simulations_num, N);
    
    % Costruzione della griglia temporale di simulazione in base alla maturity
    if maturity >= 1
        % Se la maturity è almeno 1 anno, griglia settimanale fino a "maturity" anni
        simulations_grid = businessdayoffset(today : calweeks(1) : today + calyears(maturity(j)));
    else
        % Se la maturity è inferiore a 1 anno, griglia settimanale fino a "maturity" mesi
        simulations_grid = businessdayoffset(today : calweeks(1) : today + calmonths(maturity(j)*12));
    end
    
    % Le prime date vengono escluse dalla griglia (inizia dal secondo elemento)
    simulations_grid = simulations_grid(2:end);
    Nstep = numel(simulations_grid);  % Numero di step nella simulazione
    
    % Calcolo dei coefficienti affine A e C per ogni data di simulazione
    sim_date = simulations_grid(1);
    [A, C] = affine_trick(sim_date, simulations_grid, aBach, sigmaBach, discounts_timegrid, dates, discounts);
    
    % Loop per calcolare i MtM a ogni data di simulazione
    for idx = 1:length(simulations_grid)
        sim_date = simulations_grid(idx);
        
        % Estrazione del valore simulato x_t del modello LGM per la data idx
        x_t = xtsHW(:, idx)';
        
        % Condizione se maturity >= 1 anno e sim_date prima della scadenza
        if (maturity(j) >= 1) && sim_date < businessdayoffset(today + calyears(maturity(j)))
            
            if counterparty(j) == "Epsilon"
                % Calcolo del termine esponenziale per i flussi fissi ponderati dai coefficienti affine
                term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);                  % Somma dei termini (1×N scenari)
                terminal = A(end) * exp(-C(end) * x_t);    % Termine terminale (1×N scenari)
                adjust = (idx < Nstep);                     % Aggiustamento booleano per ultimo step
                MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
                MtMs_Epsilon(:, idx) = MtM;                 % Salvataggio dei MtM per Epsilon
                
                % Shift di A e C per il prossimo step (aggiunta zero in testa)
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
                
            elseif counterparty(j) == "Delta"
                % Stesse operazioni per controparti Delta
                term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) * x_t);
                adjust = (idx < Nstep);
                MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
            end
        
        % Condizione se maturity < 1 anno e sim_date prima della scadenza in mesi
        elseif sim_date < businessdayoffset(today + calmonths(maturity(j) * 12))
            
            if counterparty(j) == "Epsilon"
                % Calcolo identico al caso precedente
                term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) * x_t);
                adjust = (idx < Nstep);
                MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
                MtMs_Epsilon(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
                
            elseif counterparty(j) == "Delta"
                % Stesse operazioni per controparti Delta
                term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) * x_t);
                adjust = (idx < Nstep);
                MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
            end
        end
    end
    
    % Salvataggio dei risultati nelle matrici TLGM a seconda della controparte
    if counterparty(j) == "Epsilon"
        TLGM_Epsilon(:, :, countEpsilon) = MtMs_Epsilon;
        countEpsilon = countEpsilon + 1;
    else
        TLGM_Delta(:, :, countDelta) = MtMs_Delta;
        countDelta = countDelta + 1;
    end

end

% Somma lungo la terza dimensione per ottenere i valori totali MtM
MtM_Epsilon = sum(TLGM_Epsilon, 3);
MtM_Delta = sum(TLGM_Delta, 3);

% Inizializzazione delle matrici per il collateral per entrambe le controparti
collateral_Epsilon = zeros(simulations_num, N);
collateral_Delta = zeros(simulations_num, N);

% Copia dei MtM totali per il calcolo del collateral
MtM_coll_Epsilon = MtM_Epsilon;
MtM_coll_Delta = MtM_Delta;

% Ricalcolo della griglia temporale e fattori di sconto forward
simulations_grid = fixed_leg_payments_dates(2:end);
disc_simgrid = interpolation(discounts, dates, simulations_grid);
fwd_disc = disc_simgrid(2:end) ./ disc_simgrid(1:end-1);

% Calcolo del collateral e aggiustamento dei MtM per la controparte Epsilon
for i = 1:N-1
    if i == 1
        % Per il primo step si mantiene il valore MtM senza modifiche
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i);
    else
        % Calcolo del fattore di sconto stocastico per il passo i
        stoch_disc = computeStochDiscHW(today, simulations_grid, i, ACT_365, a, sigmaBach, xtsHW, fwd_disc);
        % Aggiornamento del MtM includendo il collateral precedente attualizzato
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i - 1) ./ stoch_disc;
        % Calcolo del collateral corrente come opposto del MtM aggiustato
        collateral_Epsilon(:, i) = -MtM_coll_Epsilon(:, i);
        % Reset del MtM a zero dopo aggiustamento collateral
        MtM_coll_Epsilon(:, i) = 0;
    end
end

% Calcolo del collateral e aggiustamento dei MtM per la controparte Delta (analogo a Epsilon)
for i = 1:N-1
    if i == 1
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i);
    else
        stoch_disc = computeStochDiscHW(today, simulations_grid, i, ACT_365, a, sigmaBach, xtsHW, fwd_disc);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i - 1) ./ stoch_disc;
        collateral_Delta(:, i) = -MtM_coll_Delta(:, i);
        MtM_coll_Delta(:, i) = 0;
    end
end

