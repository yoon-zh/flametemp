%% Calculate adiabatic flame temperature given mix of fuels
% "yoon-zh"
% v2.preview

% Changes:
% Added input parser
% Added function tutorial (flametemp_tutorial.mlx)
% Defined functions and commented process
% Fixed problems in chemical balancing
% Vectorized multiple calculations for efficiency
% Fixed problems in interpolation function

% To do:
%   Create file "equations.pdf" and write down equations properly
%   Finish defining functions (precise, moles, etc)
%   Review precision of calculations (compare to cantera)
%   Expand data in table for available fuels
%   Consider using external database over writing all data from the book
%   Add support for fuel combustions including sulfur (ex. coal)
%   Review logic of repetitive tasks and improve if possible (ex. eqstr)

function t_flame = flametemp(fuels, varargin)
    %% Process arguments
    % Create input parser instance
    p = inputParser;
    
    % Add required parameter: fuels
    addRequired(p, 'fuels', @(x) validateFuelInput(x));
    
    % Add optional parameters with defaults
    addParameter(p, 'fuel_state', [], @(x) validateStateInput(x)); % Default handled in post-processing
    addParameter(p, 'fractions', [], @(x) validateattributes(x, {'numeric'}, {'vector', 'nonnegative', '<=', 1})); % Default handled in post-processing
    addParameter(p, 'real_air', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'phi', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
    addParameter(p, 'air_temp', 25, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -273.15}));
    addParameter(p, 'pressure', 101.325, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'xCO2', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
    addParameter(p, 'xH2O', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
    addParameter(p, 'num_comb', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
    addParameter(p, 'eff', [], @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1})); % Optional, no default
    addParameter(p, 'Wnet', [], @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'})); % Optional, no default
    addParameter(p, 'Qnet', [], @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'})); % Optional, no default
    
    % Add debugging parameters
    addParameter(p, 'printmeall', false, @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
    addParameter(p, 'precise', false, @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
    
    % Parse inputs
    parse(p, fuels, varargin{:});
    
    % Cross-validation checks
    if ~isempty(p.Results.fuel_state) && (numel(p.Results.fuels) ~= numel(p.Results.fuel_state))
        error('flametemp:InputSizeMismatch', 'fuels and fuel_state must have the same number of elements.');
    end
    if ~isempty(p.Results.fractions) && (numel(p.Results.fuels) ~= numel(p.Results.fractions))
        error('flametemp:InputSizeMismatch', 'fractions must have the same number of elements as fuels.');
    end
    if ~isempty(p.Results.fractions) && (abs(sum(p.Results.fractions) - 1) > 1e-6)
        error('flametemp:InvalidFractions', 'Fuel fractions must sum to 1.');
    end
    
    % Validate Qnet, Wnet, and eff rules
    Qnet = p.Results.Qnet;
    Wnet = p.Results.Wnet;
    eff = p.Results.eff;
    num_comb = p.Results.num_comb;
    
    % Check if any of Qnet, Wnet, or eff is provided
    if ~isempty(Qnet) || ~isempty(Wnet) || ~isempty(eff)
        % At least two of them must be provided
        if (isempty(Qnet) + isempty(Wnet) + isempty(eff)) > 1
            error('flametemp:InvalidCombination', ...
                'If any of Qnet, Wnet, or eff is provided, at least two must be provided.');
        end
        
        % Calculate missing values based on provided inputs
        if ~isempty(Qnet) && ~isempty(Wnet)
            % Calculate efficiency
            eff = Wnet / (num_comb * Qnet);
        elseif ~isempty(Wnet) && ~isempty(eff)
            % Calculate Qnet
            Qnet = Wnet / (num_comb * eff);
        elseif ~isempty(Qnet) && ~isempty(eff)
            % Calculate Wnet
            Wnet = Qnet * num_comb * eff;
        end
        
        % If all three are provided, ensure consistency
        if ~isempty(Qnet) && ~isempty(Wnet) && ~isempty(eff)
            if abs(Qnet - (Wnet / (num_comb * eff))) > 1e-6
                error('flametemp:InconsistentValues', ...
                    'Qnet, Wnet, and eff must satisfy Qnet == Wnet / (num_comb * eff).');
            end
        end
    end
    
    % Debugging outputs
    if p.Results.printmeall
        disp('All parsed inputs:');
        disp(p.Results);
    end
    if p.Results.precise
        disp('Precise calculations enabled.');
    end
    % C_TO_KELVIN_ADD = 273.15; % ADD TO CONVERT °C TO K % Unused for now
    KELVIN_TO_C_ADD = -273.15; % ADD TO CONVERT K TO °C

    % Move results to modifiable struct with shorter name "t"
    t = p.Results;
    
    % Post-process defaults for fuel_state and fractions
    if isempty(t.fuel_state)
        % Default fuel_state to '(g)' for all fuels
        t.fuel_state = repmat("(g)", size(t.fuels));
    end
    if isempty(t.fractions)
        % Default fractions to equal proportions
        t.fractions = ones(size(t.fuels)) / numel(t.fuels);
    end

    %% Calculations
    % Initialize fuel compositions map
    fuel_comp = initFuelCompositions();
    
    % Obtain moles of air
    n_air = molesAir(t.fuels, t.fractions, t.real_air, fuel_comp);

    % Obtain moles of H2O from humid air
    nH2Or = molesH2OinAir(t.phi, t.air_temp, t.pressure, n_air);

    % Array of moles of reactants and products
    n_r = molesReactants(t.fractions, n_air, nH2Or);
    n_p = molesProducts(t.fuels, t.fractions, t.real_air, t.xCO2, t.xH2O, fuel_comp, n_r);

    % Print chemical formula
    dispChemEq(t.fuels, t.fuel_state, n_r, n_p);

    % Obtain enthalpy of reactants
    h_r = enthalpyReactants(t.fuels, t.fuel_state, t.air_temp);
    
    % Obtain base enthalpy of products, excluding h @ AFT
    h_pBase = enthalpyProductsBase();

    % Calculate the enthalpy @ AFT assuming all elements are N2
    h_flame = enthalpyAFT(n_r, n_p, h_r, h_pBase);

    % Interpolate for t_flame (IN KELVIN)
    t_flame = AFTintp(n_p, h_flame);

    % Print AFT in K and °C
    dispAFT(t_flame, KELVIN_TO_C_ADD);

    %% TO IMPLEMENT
    % If precise value required, calculate with 100 extra iterations
    % if t.precise
    %     t_flame = preciseAFT(t_flame);
    % end
    
    % if ~isempty(Qnet) || ~isempty(Wnet) || ~isempty(eff)
    %     molar_fractions = molarFractions();
    %     disp(molar_fractions)
    % end

    % Debug prints
    % if t.printmeall
    %     disp("Moles of air: " + n_air)
    % end
end

%% Validate fuels input
function validateFuelInput(x)
    if ~(isstring(x) || iscellstr(x)) || isempty(x)
        error('fuels must be a non-empty string array or cell array of characters.');
    end
end

%% Validate fuel_state input
function validateStateInput(x)
    validStates = {'(l)', '(g)'};
    if ~isempty(x) && ~(isstring(x) || iscellstr(x)) || ~all(ismember(x, validStates))
        error('fuel_state must contain only "(l)" or "(g)" for each fuel.');
    end
end

%% Define fuel compositions (element count)
function fuel_comp = initFuelCompositions() % [C H O N]
    fuel_comp = containers.Map(...
    {'C2H2','C2H5OH','CH3OH','C4H10','C12H26','C2H6','CH4','C8H18','C3H8','NH3', 'C3H6'},...
    {[2,2,0,0],[2,6,1,0],[1,4,1,0],[4,10,0,0],[12,26,0,0],[2,6,0,0],[1,4,0,0],[8,18,0,0],[3,8,0,0],[0,3,0,1],[3,6,0,0]});
end

%% Obtain moles of air
function n_air = molesAir(fuels, fractions, real_air, fuel_comp)
    % In equations.pdf see equation 1
    comp = cell2mat(cellfun(@(f) fuel_comp(f), cellstr(fuels), 'UniformOutput', false).');
    
    % Calculate terms for each fuel
    terms = comp(:, 1) + comp(:, 2)/4 - comp(:, 3)/2;
    
    % Compute ideal air
    ideal_air = sum(fractions(:) .* terms);
    
    % Calculate total moles of air
    n_air = ideal_air * real_air;
end

%% Calculate moles of water in reactant from humid air
function nH2Or = molesH2OinAir(phi, air_temp, Pt, n_air)
    data = readmatrix('TABLES.xlsx', 'Sheet', 'Tsat');
    if ismember(air_temp, data(:,1)) % If exact tsat exists
        indx = air_temp == data(:,1);
        psat = data(indx, 2); % Return exact psat
    else % If not, interpolate
        T_indices = find(data(:,1) <= air_temp);
        to_indx = T_indices(end);
        tf_indx = find(data(:,1) > air_temp, 1);
        t_intp = [data(to_indx, 1); data(tf_indx, 1)];
        p_intp = [data(to_indx, 2); data(tf_indx, 2)];
        psat = LagrangeIntp(t_intp, p_intp, air_temp);
    end
    % In equations.pdf see equation 2
    Pv = psat*phi;
    nH2Or = (Pv/Pt) * n_air/(1-(Pv/Pt));
end

%% Calculate moles of reactants
function n_r = molesReactants(fractions, n_air, nH2Or)
    % In equations.pdf see equation 3
    nO2r = n_air;
    nN2r = 3.76 * n_air;
    n_r = [fractions, nO2r, nN2r, nH2Or]; % Store moles of reactants
end

%% Balance chemical equation
function n_p = molesProducts(fuels, fractions, real_air, xCO2, xH2O, fuel_comp, n_r)
    nO2r = n_r(numel(fuels) + 1);
    nN2r = n_r(numel(fuels) + 2);
    nH2Or = n_r(numel(fuels) + 3);
    % Retrieve element count of fuels [C H O N]
    fuel_matrix = arrayfun(@(f) fuel_comp(f), fuels, 'UniformOutput', false);
    % Convert cell to matrix, now it's Nx4
    fuel_matrix = cell2mat(fuel_matrix');
    fractions = fractions(:);  % force fractions into Nx1 column
    nC = sum(fractions .* fuel_matrix(:,1));
    nH = sum(fractions .* fuel_matrix(:,2)) + 2*nH2Or;
    nO = sum(fractions .* fuel_matrix(:,3)) + 2*nO2r + nH2Or;
    nN = sum(fractions .* fuel_matrix(:,4)) + 2*nN2r;
    % In equations.pdf see equation 4
    if real_air >= 1
        %    nCO2       nCO         nH2O        nH2     nO2
        A = [1,         1,          0,          0,      0; % C
             0,         0,          2,          2,      0; % H
             2,         1,          1,          0,      2; % O
             1-xCO2,-xCO2,          0,          0,      0; % CO2 conversion
             0,         0,     1-xH2O,      -xH2O,      0]; % H2O conversion
        B = [nC;        nH;        nO;          0;      0];
        sol = A\B;
        nO2p = sol(5);
    elseif real_air < 1
        %    nCO2       nCO         nH2O        nH2
        A = [1,         1,          0,          0; % C
             0,         0,          2,          2; % H
             2,         1,          1,          0; % O
             1-xCO2,-xCO2,          0,          0; % CO2 conversion
             0,         0,     1-xH2O,      -xH2O]; % H2O conversion
        B = [nC;        nH;        nO;          0];
        sol = A\B;
        nO2p = 0;
    end
    %     [nCO2p,    nCOp,  nH2Op,   nH2p, nO2p,   nN2p]
    n_p = [sol(1), sol(2), sol(3), sol(4), nO2p, (nN/2)];
end

%% Print chemical equation
function dispChemEq(fuels, fuel_state, n_r, n_p)
    eqstr = "";
    numFuels = numel(fuels);

    % Reactants
    for i = 1:numFuels
        if n_r(i) ~= 1
            eqstr = eqstr + n_r(i) + " ";
        end
        eqstr = eqstr + fuels(i) + fuel_state(i) + " + ";
    end
    if n_r(end) ~= 0 % nH2Or = 0 when the air is dry
        eqstr = eqstr + n_r(numFuels+1) + " (O2 + 3.76N2) + " + n_r(end) + " H2O => ";
    else
        eqstr = eqstr + n_r(numFuels+1) + " (O2 + 3.76N2) => ";
    end

    % Products
    eqstr = eqstr + n_p(1) + " CO2 + ";
    if n_p(2) ~= 0 % nCOp = 0 when xCO2 = 1
        eqstr = eqstr + n_p(2) + " CO + ";
    end

    eqstr = eqstr + n_p(3) + " H2O + ";
    if n_p(4) ~= 0 % nH2p = 0 when xH2O = 1
        eqstr = eqstr + n_p(4) + " H2 + ";
    end

    if n_p(5) ~= 0 % nO2p = 0 when real_air >= 1
        eqstr = eqstr + n_p(5)+" O2 + ";
    end
    eqstr = eqstr + n_p(6)+" N2";
    disp(eqstr)
end

%% Matrix for enthalpy of reactants
function h_r = enthalpyReactants(fuels, fuel_state, air_temp)
    % First get enthalpy of formation for fuels
    % Which react at T=25°C (h@25 - hamb = 0)
    hf_fuels = callhf(fuels, fuel_state);
    
    % Then calculate enthalpy for components in air
    air_comp = ["O2", "N2", "H2O"];
    hf_air = callhf(air_comp);
    h_air = callh(air_comp, air_temp);
    hamb_air = callhamb(air_comp);
    % In equations.pdf see equation 5
    h_total_air = hf_air + h_air - hamb_air;
    % Return matrix of enthalpies
    h_r = [hf_fuels, h_total_air];
end

%% Matrix for enthalpy of products excluding enthalpy of AFT
function h_pBase = enthalpyProductsBase()
    PRODUCTS_STR = ["CO2", "CO", "H2O", "H2", "O2", "N2"];
    hf_p = callhf(PRODUCTS_STR);
    hamb_p = callhamb(PRODUCTS_STR);
    h_pBase = hf_p - hamb_p;
end

%% Calculate enthalpy value assuming all products are equal
function h_flame = enthalpyAFT(n_r, n_p, h_r, h_pBase)
    % In equations.pdf see equation 6
    h_flame = (sum(n_r.*h_r) - sum(n_p.*h_pBase)) / sum(n_p);
end

%% Calculate AFT from enthalpy
function t_flame = AFTintp(n_p, h_flame)
    AFT_OFFSET = 0.8; % Offset to interpolate for AFT
    % First assume everything is N2 to find a temperature of reference
    t_max = tfind("N2", h_flame);
    
    % Consider t_max and a lower temperature to interpolate t_flame
    h_max_total = hflame(t_max, n_p);
    t_low = AFT_OFFSET*t_max;
    h_low_total = hflame(t_low, n_p);
    h_range = [h_low_total h_max_total];
    t_range = [t_low t_max];
    % Multiply h_flame to the sum of moles of products, since we account
    % now for their contribution (this is total enthalpy, in kJ/K)
    h_toFind = h_flame * sum(n_p);

    % In equations.pdf see equation 7
    t_flame = LagrangeIntp(h_range, t_range, h_toFind);
end

%% Print results of AFT
function dispAFT(t_flame, KELVIN_TO_C_ADD)
    disp("AFT: " + t_flame + " K")
    t_flame = t_flame + KELVIN_TO_C_ADD;
    disp("AFT: " + t_flame + " °C")
end

%% TO IMPLEMENT

%% Calculate precise t_flame
% function aft_precise = preciseAFT(t_flame)
%     t_graph = linspace(t_low, t_max);
%     h_graph = zeros(size(t_graph));
%     for i = 1:numel(t_graph)
%         h_graph(i) = hflame(t_graph(i), n_p);
%     end
%     hprec_ind = find(h_graph <= val_tofind);
%     ho_inx = hprec_ind(end);
%     hf_inx = find(h_graph > val_tofind, 1);
%     if numel(hf_inx) == 0
%         t_flame_precs = t_graph(ho_inx(end));
%     else
%         tfl_intp = [t_graph(ho_inx); t_graph(hf_inx)];
%         hfl_intp = [h_graph(ho_inx); h_graph(hf_inx)];
%         t_flame_precs = LagrangeIntp(hfl_intp, tfl_intp, val_tofind);
%     end
%     disp("Precise adiabatic flame temperature: " + t_flame_precs + " K")
%     disp("Precise T° in Celsius: " + (t_flame_precs-273.15) + " °C")
%     figure(123)
%     hold on
%     grid on
%     plot(t_graph, h_graph, "-b.")
%     plot(t_flame_precs, val_tofind, "rs",...
%         'MarkerSize',10,...
%         'MarkerEdgeColor','m',...
%         'MarkerFaceColor',[0.5,0.5,0.5])
%     xlabel("T (K)")
%     ylabel("h (kJ)")
%     title("Adiabatic flame temperature vs enthalpy")
%     legend("Enthalpies", "Precise T°")
%     hold off
%     aft_precise = t_flame;
% end

%% Calculate molar fractions
% function molar_fractions = molarFractions()
%     PRODUCTS_STR = ["CO2", "CO", "H2O", "H2", "O2", "N2"];
%     molar_fractions = n_p/sum(n_p);
%     disp("Molar fractions: ")
%     for i = 1:numel(PRODUCTS_STR)
%         disp(PRODUCTS_STR(i) + ": " + molar_fractions(i))
%     end
% end

%% Calculate temperature of exiting products
                % CO2       CO      H2O         H2      O2      N2 
% molar_weight = [44.0095, 28.0101, 18.01528, 2.01588, 31.9988, 28.0134]; % kg/kmol
% mass_p = molar_weight .* molar_fractions;

% %% Find molar flow of combustion
% Q = Wnet/(2*eff); % For 2 combusters, in MW
% disp("Q for each combustor: " + Q + " MW");
% 
% % Calculate Cp assuming ideal gas
% %        CO2       CO     H2O      H2       O2       N2
% acp = [22.26,   28.16,  32.24,   29.11,   25.48,    28.9];
% bcp = [5.981,  0.1675, 0.1923, -0.1916,    1.52, -0.1571] * (10^(-2));
% ccp = [-3.501, 0.5372,  1.055,  0.4003, -0.7155,  0.8081] * (10^(-5));
% dcp = [7.469,  -2.222, -3.595, -0.8704,   1.312,  -2.873] * (10^(-9));
% 
% Cp_p = zeros(size(products));
% for i = 1:numel(products)
%     Cp_p(i) = acp(i) + bcp(i)*t_flame + ccp(i)*(t_flame^2) + dcp(i)*(t_flame^3);
%     Cp_p(i) = Cp_p(i) * molar_fractions(i);
% end
% Cp = sum(Cp_p);
% 
% if printmeall
%     disp("Cp of each component: ")
%     disp([products', Cp_p'])
%     disp("Cp final: " + Cp + " kJ/kmol*K")
% end
% 
% texit = linspace(0.25*t_flame, 0.75*t_flame, 25);
% molar_flow = (10^3) * Q./(Cp.*(t_flame - texit));
% 
% figure()
% grid on
% plot(texit, molar_flow, "-b.")
% xlabel("Exit temperature (K)")
% ylabel("Molar flow (kmol/s)")
% title("Molar flow vs exit T°")
% if printmeall
%     disp("________________________")
%     varNames = ["Exit temperature (K)", "Molar flow (kmol/s)"];
%     disp(table(texit', molar_flow', 'VariableNames', varNames))
% end
% 
% % In Celsius
% texitC = texit - 273.15;
% figure()
% grid on
% plot(texitC, molar_flow, "-b.")
% xlabel("Exit temperature (°C)")
% ylabel("Molar flow (kmol/s)")
% title("Molar flow vs exit T°")
% if printmeall
%     disp("________________________")
%     varNames = ["Exit temperature (°C)", "Molar flow (kmol/s)"];
%     disp(table(texitC', molar_flow', 'VariableNames', varNames))
% end


%% Fuctions
% Find enthalpies
function h_r = hfind(search_var, propty, t_test)
    if propty == "hf" % Find enthalpy of formation
        if exist("t_test", "var") % User can input (l), (g) into the t_test var
            if isa(t_test, "string")
                search_var = search_var + t_test;
            end
        elseif ~contains(search_var, "(") % If user does not specify state
            search_var = search_var + "(g)"; % Default search to gas
        end
        comp = string(readcell("TABLES.xlsx", "Sheet","H°","Range","A2:A31"));
        hf = readmatrix("TABLES.xlsx", "Sheet","H°","Range","B2:B31");
        indices = comp == search_var;
        if ~any(indices)
            error('flametemp:stringNotFound', 'Enthalpy data for one or more of the fuels introduced does not exist in the database.')
        end
        h_r = hf(indices);
    
    elseif propty == "h" % Find enthalpy at given temperature
        data = readmatrix('TABLES.xlsx', 'Sheet', search_var);
        if ismember(t_test, data(:,1)) % If exact temperature exists
            indx = t_test == data(:,1);
            h_r = data(indx, 2); % Return exact enthalpy
        else % If no exact temperature exists, interpolate
            t_indices = find(data(:,1) <= t_test);
            to_indx = t_indices(end);
            tf_indx = find(data(:,1) > t_test, 1);
            t_intp = [data(to_indx, 1); data(tf_indx, 1)];
            h_intp = [data(to_indx, 2); data(tf_indx, 2)];
            h_r = LagrangeIntp(t_intp, h_intp, t_test);
        end
    end
end

% Find temperature
function t_r = tfind(search_var, h_test)
    data = readmatrix('TABLES.xlsx', 'Sheet', search_var);
    if ismember(h_test, data(:,2)) % If exact enthalpy exists
        indx = h_test == data(:,2);
        t_r = data(indx, 1); % Return exact temperature
    else % If not, interpolate
        h_indices = find(data(:,2) <= h_test);
        ho_indx = h_indices(end);
        hf_indx = find(data(:,2) > h_test, 1);
        h_intp = [data(ho_indx, 2); data(hf_indx, 2)];
        t_intp = [data(ho_indx, 1); data(hf_indx, 1)];
        t_r = LagrangeIntp(h_intp, t_intp, h_test);
    end
end

% Find enthalpy of flame given temperature
function h_fl = hflame(t_flame, n_p)
    PRODUCTS_STR = ["CO2", "CO", "H2O", "H2", "O2", "N2"];
    h_p = callh(PRODUCTS_STR, t_flame);
    h_fl = sum(n_p.*h_p);
end

% Calls enthalpies of formation
function hf = callhf(elements, elem_state)
    hf = zeros(1, numel(elements));
    for i = 1:numel(elements)
        if exist("elem_state", "var")
            hf(i) = hfind(elements(i), "hf", elem_state(i));
        else
            hf(i) = hfind(elements(i), "hf");
        end
    end
end

% Calls enthalpies at 298 K
function hamb = callhamb(elements)
    hamb = zeros(1, numel(elements));
    for i = 1:numel(elements)
        hamb(i) = hfind(elements(i), "h", 298);
    end
end

% Calls enthalpies at given temperature
function h = callh(elements, t_flame)
    h = zeros(1, numel(elements));
    for i = 1:numel(elements)
        h(i) = hfind(elements(i), "h", t_flame);
    end
end

% Lagrange interpolation from vectors x, y and test value x_intp
function y_intp = LagrangeIntp(x, y, x_intp)
    % Check input validity
    n = numel(x);
    if n < 2
        error('LagrangeIntp:InputSizeError', 'x must have at least two elements.');
    end
    if numel(y) ~= n
        error('LagrangeIntp:InputSizeError', 'x and y must be the same length.');
    end
    if numel(unique(x)) ~= n
        error('LagrangeIntp:DuplicateXValues', 'All x values must be distinct.');
    end
    
    % Ensure column vectors
    x = x(:);
    y = y(:);
    
    % Compute Lagrange interpolation
    y_intp = 0;
    for i = 1:n
        xi = x(i);
        mask = (1:n) ~= i;  % Logical mask excluding current index
        x_diff = xi - x(mask);
        terms = (x_intp - x(mask)) ./ x_diff;
        L_i = prod(terms);
        y_intp = y_intp + y(i) * L_i;
    end
end

% j.porras 
% v2.preview