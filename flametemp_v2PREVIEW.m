%% Calculate adiabatic flame temperature given mix of fuels
% "yoon-zh"
% v2.preview

% Changes:
% Added input parser
% Added function tutorial (flametemp_tutorial.mlx)
% Defined a few functions

% To do:
%   Finish defining functions
%   Cleanup functions and overall code, divide and simplify
%   Verify precision through testing and comparison with external examples
%   Export tables to a broader set of tables. Consider using external
%   database over writing all data from the book.

% precisetflame = true; % Uncomment to calculate a precise value for tflame

%% def flametemp
function tflame = flametemp(fuels, varargin)
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
    
    % Post-process defaults for fuel_state and fractions
    if isempty(p.Results.fuel_state)
        % Default fuel_state to '(g)' for all fuels
        p.Results.fuel_state = repmat("(g)", size(p.Results.fuels));
    end
    if isempty(p.Results.fractions)
        % Default fractions to equal proportions
        p.Results.fractions = ones(size(p.Results.fuels)) / numel(p.Results.fuels);
    end
    
    % Cross-validation checks
    if numel(p.Results.fuels) ~= numel(p.Results.fuel_state)
        error('flametemp:InputSizeMismatch', 'fuels and fuel_state must have the same number of elements.');
    end
    if numel(p.Results.fuels) ~= numel(p.Results.fractions)
        error('flametemp:InputSizeMismatch', 'fractions must have the same number of elements as fuels.');
    end
    if abs(sum(p.Results.fractions) - 1) > 1e-6
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
    
    % Obtain moles of air
    n_air = molesAir(fuels, real_air);

    % Obtain moles of H2O from humid air
    nH2Or = molesH2OinAir(n_air, phi, pressure);

    % Array of moles of reactants and products
    n_r = molesReactants(fractions, n_air, nH2Or);
    n_p = molesProducts(fuels, fractions, n_r, real_air);

    % Print chemical formula
    dispChemEq(fuels, fraction, fuel_state, n_r, n_p)

    
    

    % Placeholder return (replace with actual calculation)
    tflame = 1; 
    disp("Parse test")
end

% Validate fuels input
function validateFuelInput(x)
    if ~(isstring(x) || iscellstr(x)) || isempty(x)
        error('fuels must be a non-empty string array or cell array of characters.');
    end
end

% Validate fuel_state input
function validateStateInput(x)
    validStates = {'(l)', '(g)'};
    if ~isempty(x) && ~(isstring(x) || iscellstr(x)) || ~all(ismember(x, validStates))
        error('fuel_state must contain only "(l)" or "(g)" for each fuel.');
    end
end

%% Define compositions of air



%% Obtain moles of air
function n_air = molesAir(fuels, real_air)
    % Define the fuels and their compositions - struct("fuel" [C H O N]) 
    fuel_compositions = struct("C2H2", [2, 2, 0, 0], "C2H5OH", [2, 6, 1, 0],...
        "CH3OH", [1, 4, 1, 0], "C4H10", [4, 10, 0, 0], "C12H26", [12, 26, 0, 0],...
        "C2H6", [2, 6, 0, 0], "CH4", [1, 4, 0, 0], "C8H18", [8, 18, 0, 0], ...
        "C3H8", [3, 8, 0, 0], "NH3", [0, 3, 0, 1]);
    
    % Balance air
    ideal_air = 0;
    for i = 1:length(fuels)
        ideal_air = ideal_air + fractions(i) * (fuel_compositions.(fuels(i))(1) + fuel_compositions.(fuels(i))(2)/4 - fuel_compositions.(fuels(i))(3)/2);
    end
    n_air = ideal_air * real_air;
    if printmeall
        disp("Ideal air amount: " + ideal_air + " moles")
    end
end

%% Calculate moles of water in reactant from humid air
function nH2Or = molesH2OinAir(n_air, phi, Pt)
    data = readmatrix('TABLES.xlsx', 'Sheet', 'Tsat');
    if ismember(T, data(:,1)) % If exact tsat exists
        indx = T == data(:,1);
        psat = data(indx, 2); % Return exact psat
    else % If not, interpolate
        T_indices = find(data(:,1) <= T);
        to_indx = T_indices(end);
        tf_indx = find(data(:,1) > T, 1);
        t_intp = [data(to_indx, 1); data(tf_indx, 1)];
        p_intp = [data(to_indx, 2); data(tf_indx, 2)];
        psat = LagrangeIntp(t_intp, p_intp, T);
    end
    Pv = psat*phi;
    nH2Or = (Pv/Pt) * n_air/(1-(Pv/Pt));
    if printmeall
        disp("Pressure of saturation: " + double(psat) + " kPa")
    end
end

%% Calculate moles of reactants
function n_r = molesReactants(fractions, n_air)
    nO2r = n_air;
    nN2r = 3.76 * n_air;
    n_r = [fractions, nO2r, nN2r, nH2Or]; % Store moles of reactants
end

%% Balance chemical equation
function n_p = molesProducts(fuels, fractions, n_r, real_air)
    nO2r = n_r(length(fuels));
    nN2r = n_r(length(fuels) + 1);
    nH2Or = n_r(length(fuels)) + 2;
    nC = 0; nH = 0; nO = 0;
    for i = 1:length(fuels)
        nC = nC + fractions(i) * fuel_compositions.(fuels(i))(1);
        nH = nH + fractions(i) * fuel_compositions.(fuels(i))(2);
        nO = nO + fractions(i) * fuel_compositions.(fuels(i))(3);
    end
    if real_air >= 1
        A = [1,      1,      0,  0, 0; % C
            0,       0,      2,  2, 0; % H
            2,       1,      1,  0, 2; % O
            1-xCO2, -1,      0,  0, 0; % CO2 conversion
            0,       0, 1-xH2O, -1, 0]; % H2O conversion
        B = [nC; % C
        nH + 2*nH2Or; % H
        nO + 2*nO2r + nH2Or; % O
        0; 0];
        sol = A\B;
        nO2p = sol(5);
    elseif real_air < 1
        A = [1,      1,    0,    0; % C
            0,       0,    2,    2; % H
            1-xCO2, -1,    0,    0; % CO2 conversion
            0,       0, 1-xH2O, -1]; % H2O conversion
        B = [nC; % C
            nH + 2*nH2Or; % H
            0; 0];
        sol = A\B;
        nO2p = 0;
    end
    
    nCO2p = sol(1);
    nCOp = sol(2);
    nH2Op = sol(3);
    nH2p = sol(4);
    if printmeall
        disp("Chemical balance matrix: ")
        disp([A B])
    end
    nN2p = nN2r;
    n_p = [nCO2p, nCOp, nH2Op, nH2p, nO2p, nN2p]; % Store moles of products
end

function dispChemEq(fuels, fractions, fuel_state, n_r, n_p)
    eqstr = "";
    for i = 1:length(fuels)
        eqstr = eqstr + fractions(i) + fuels(i) + fuel_state(i) + " + ";
    end
    eqstr = eqstr + n_r(0) + "(O2 + 3.76N2) + " + nH2Or + "H2O => ";
    eqstr = eqstr + n_p(0)+"CO2 + "+n_p(1)+"CO + "+n_p(2)+"H2O + "+n_p(3)+"H2 + ";
    if n_p(4) ~= 0
        eqstr = eqstr + n_p(4)+"O2 + ";
    end
    eqstr = eqstr + n_p(5)+"N2";
    disp("Real combustion process: ")
    disp(eqstr)
end

%% Determine value to find
% Reactants
hf_fuels = callhf(fuels, fuel_state);
reactants_sum = 0;
for i = 1:length(fuels) % First sum fuels, which react at T=25°C
    reactants_sum = reactants_sum + fractions(i)*hf_fuels(i);
end
hf_air = callhf(air_comp);
h_air = callh(air_comp, T+273);
hamb_air = callhamb(air_comp);
for i = 1:length(air_comp) % Then, sum air
    reactants_sum = reactants_sum + n_air(i)*(hf_air(i) + h_air(i) - hamb_air(i));
end
if printmeall
    disp("Reactants: " + reactants_sum)
end

% Products, excluding hflame
hf_p = callhf(products);
hamb_p = callhamb(products);
products_sum = 0;
for i = 1:length(products)
    products_sum = products_sum + n_p(i)*(hf_p(i) - hamb_p(i));
end
val_tofind = reactants_sum - products_sum;
if printmeall
    disp("Products, excluding h @ tflame: " + products_sum)
    disp("Value to find: " + val_tofind)
end

%% Assume all products are N2, find Tmax
h_N2ideal = val_tofind / sum(n_p);
Tmax = tfind("N2", h_N2ideal);
disp("Max T° = " + Tmax + " K")

%% Consider Tmax and a lower temperature to interpolate tflame
hmax = hflame(Tmax, products, n_p);
Tlow = 0.8*Tmax;
hlow = hflame(Tlow, products, n_p);
h_tofind = [hlow hmax];
T_tofind = [Tlow Tmax];
if printmeall
    disp("h to interpolate:")
    disp(h_tofind)
    disp("T to interpolate:")
    disp(T_tofind)
end
tflame = LagrangeIntp(h_tofind, T_tofind, val_tofind);
disp("Adiabatic flame temperature: " + tflame + " K")
tflame_C = tflame - 273.15; % C
disp("T° in Celsius: " + tflame_C + " °C")

%% Calculate precise tflame
if precise
    t_graph = linspace(Tlow, Tmax);
    h_graph = zeros(size(t_graph));
    for i = 1:numel(t_graph)
        h_graph(i) = hflame(t_graph(i), products, n_p);
    end
    hprec_ind = find(h_graph <= val_tofind);
    ho_inx = hprec_ind(end);
    hf_inx = find(h_graph > val_tofind, 1);
    if numel(hf_inx) == 0
        tflame_precs = t_graph(ho_inx(end));
    else
        tfl_intp = [t_graph(ho_inx); t_graph(hf_inx)];
        hfl_intp = [h_graph(ho_inx); h_graph(hf_inx)];
        tflame_precs = LagrangeIntp(hfl_intp, tfl_intp, val_tofind);
    end
    disp("Precise adiabatic flame temperature: " + tflame_precs + " K")
    disp("Precise T° in Celsius: " + (tflame_precs-273.15) + " °C")
    figure(123)
    hold on
    grid on
    plot(t_graph, h_graph, "-b.")
    plot(tflame_precs, val_tofind, "rs",...
        'MarkerSize',10,...
        'MarkerEdgeColor','m',...
        'MarkerFaceColor',[0.5,0.5,0.5])
    xlabel("T (K)")
    ylabel("h (kJ)")
    title("Adiabatic flame temperature vs enthalpy")
    legend("Enthalpies", "Precise T°")
    hold off
    tflame = tflame_precs;
end

%% Calculate molar fractions
molar_fractions = n_p/sum(n_p);
disp("Molar fractions: ")
for i = 1:length(products)
    disp(products(i) + ": " + molar_fractions(i))
end

%% Calculate temperature of exiting products
                % CO2       CO      H2O         H2      O2      N2 
% molar_weight = [44.0095, 28.0101, 18.01528, 2.01588, 31.9988, 28.0134]; % kg/kmol
% mass_p = molar_weight .* molar_fractions;

%% Find molar flow of combustion
Q = Wnet/(2*eff); % For 2 combusters, in MW
disp("Q for each combustor: " + Q + " MW");

% Calculate Cp assuming ideal gas
%        CO2       CO     H2O      H2       O2       N2
acp = [22.26,   28.16,  32.24,   29.11,   25.48,    28.9];
bcp = [5.981,  0.1675, 0.1923, -0.1916,    1.52, -0.1571] * (10^(-2));
ccp = [-3.501, 0.5372,  1.055,  0.4003, -0.7155,  0.8081] * (10^(-5));
dcp = [7.469,  -2.222, -3.595, -0.8704,   1.312,  -2.873] * (10^(-9));

Cp_p = zeros(size(products));
for i = 1:length(products)
    Cp_p(i) = acp(i) + bcp(i)*tflame + ccp(i)*(tflame^2) + dcp(i)*(tflame^3);
    Cp_p(i) = Cp_p(i) * molar_fractions(i);
end
Cp = sum(Cp_p);

if printmeall
    disp("Cp of each component: ")
    disp([products', Cp_p'])
    disp("Cp final: " + Cp + " kJ/kmol*K")
end

texit = linspace(0.25*tflame, 0.75*tflame, 25);
molar_flow = (10^3) * Q./(Cp.*(tflame - texit));

figure()
grid on
plot(texit, molar_flow, "-b.")
xlabel("Exit temperature (K)")
ylabel("Molar flow (kmol/s)")
title("Molar flow vs exit T°")
if printmeall
    disp("________________________")
    varNames = ["Exit temperature (K)", "Molar flow (kmol/s)"];
    disp(table(texit', molar_flow', 'VariableNames', varNames))
end

% In Celsius
texitC = texit - 273.15;
figure()
grid on
plot(texitC, molar_flow, "-b.")
xlabel("Exit temperature (°C)")
ylabel("Molar flow (kmol/s)")
title("Molar flow vs exit T°")
if printmeall
    disp("________________________")
    varNames = ["Exit temperature (°C)", "Molar flow (kmol/s)"];
    disp(table(texitC', molar_flow', 'VariableNames', varNames))
end


%% Fuctions
function h_r = hfind(search_var, propty, t_test) % Find enthalpies
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


function t_r = tfind(search_var, h_test) % Find temperature
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


function h_fl = hflame(tflame, products, n_p) % Find enthalpy of flame given temperature
    h = callh(products, tflame);
    h_fl = 0;
    for i = 1:length(products)
        h_fl = h_fl + n_p(i)*h(i);
    end
end


function hf = callhf(elements, elem_state) % Calls enthalpies of formation
hf = zeros(1, length(elements));
for i = 1:length(elements)
    if exist("elem_state", "var")
        hf(i) = hfind(elements(i), "hf", elem_state(i));
    else
        hf(i) = hfind(elements(i), "hf");
    end
end
end


function hamb = callhamb(elements) % Calls enthalpies at 298 K
hamb = zeros(1, length(elements));
for i = 1:length(elements)
    hamb(i) = hfind(elements(i), "h", 298);
end
end


function h = callh(elements, tflame) % Calls enthalpies at given temperature
h = zeros(1, length(elements));
for i = 1:length(elements)
    h(i) = hfind(elements(i), "h", tflame);
end
end


function y_intp = LagrangeIntp(x, y, x_intp) % Lagrange interpolation
    % Receives 2x1 vectors "x", "y", value to interpolate "x_intp"
    % Returns interpolation value "y_intp"
    if numel(x) < 2 || numel(y) < 2
        disp("x: ")
        disp(x)
        disp("y: ")
        disp(y)
        disp("There's a mistake somewhere. Check your code.")
    end
    if size(x) ~= [numel(x),1]
        x = x';
    end
    if size(y) ~= [numel(y),1]
        y = y';
    end
    L = 0;
    for i = 1:2
        ph = 1;
        for j = 1:2
            if i ~= j
                ph = conv(ph,poly(x(j,1))) / (x(i,1) - x(j,1));
            end
        end
        L = L + ph*y(i,1);
    end
    y_intp = polyval(L,x_intp); % Interpolates x_test
end

% j.porras 
% v2.preview