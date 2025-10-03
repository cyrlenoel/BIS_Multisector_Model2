function io = f_setup_bea(folder_path, mapping_file, source, iso2, io_raw)
% -------------------------------------------------------------------------
% This function processes US BEA IO data for the US and maps it into
% broader industries based on the provided mapping file. It performs the
% following tasks:
%
% 1. Maps the data to broader industry classifications as defined in the
%    mapping file (\input\data\IO_mapping.xlsx).
%
% 2. Saves the following variables into a structured format (`io`) for
%    further use in the model:
%    a. Country code: `iso2'
%    a. Names of industries in the matrix order: `industry_names`
%    b. Total number of sectors: `nsectors`
%    c. Matrix of intermediate inputs and factor shares: `io_inputs`
%       - The top rows of the matrix represent intermediate input shares.
%       - The bottom rows represent the shares of labour and capital as
%         factors of production.
%    d. Weights of final demand components by industry:
%       - `consumption_weights`, `investment_weights`, and `government_weights`
%    e. Shares of total final demand components in GDP:
%       - `consumption_share`, `investment_share`, and `government_share`
%    f. Share of taxes in gross output: `taxes_ishare`
%
% 4. Saves the `io` struct into a `.mat` file under the `\input` folder.
% -------------------------------------------------------------------------

%% 1 Load and setup data

t_io = io_raw;
t_io.Code = t_io.Properties.RowNames;

% 1.1 Process rows
% 1.1.2 Load broader industry codes and join with the table
row_codes = readtable(mapping_file, 'Sheet', ['IO_RowCodes_' upper(source)]);
var_codes = readtable(mapping_file, 'Sheet', ['IO_VariableCodes_' upper(source)]);
t_io_brd = innerjoin(t_io,row_codes,'Keys','Code');

% 1.1.3 Aggregate rows according to the broader industry classification
t_io_brd(:, {'Code', 'Definition'}) = [];
t_io_brd = groupsummary(t_io_brd,"Code_std_20","sum");
t_io_brd.GroupCount = [];

% 1.1.4 Sort
row_codes_brd = table(unique(row_codes.Code_std_20,'stable'));
row_codes_brd.Properties.VariableNames = {'Code_std_20'};
row_codes_brd.Code_nums = (1:size(row_codes_brd,1))';
t_io_brd = innerjoin(t_io_brd,row_codes_brd);
t_io_brd = sortrows(t_io_brd,"Code_nums");
t_io_brd.Properties.VariableNames = replace(t_io_brd.Properties.VariableNames,"sum_","");
t_io_brd.Properties.RowNames = t_io_brd.Code_std_20;
t_io_brd(:,{'Code_std_20', 'Code_nums'}) = [];

% 1.2 Process columns/variables
% 1.2.1 Transpose table t_io_sum_brd
t_io_brd_tr = rows2vars(t_io_brd);
t_io_brd_tr.Properties.VariableNames = replace(t_io_brd_tr.Properties.VariableNames,'OriginalVariableNames','Code');

% 1.2.2 Aggregate rows according to the broader industry classification
t_io_brd_tr_brd = innerjoin(t_io_brd_tr,var_codes,'Keys','Code');
t_io_brd_tr_brd(:,{'Code','Definition'}) = [];
t_io_brd_tr_brd = groupsummary(t_io_brd_tr_brd,"Code_std_20","sum");
t_io_brd_tr_brd.GroupCount = [];

% 1.2.3 Sort
var_codes_brd = table(unique(var_codes.Code_std_20,'stable'));
var_codes_brd.Properties.VariableNames = {'Code_std_20'};
var_codes_brd.Code_nums = (1:size(var_codes_brd,1))';
t_io_brd_tr_brd = innerjoin(t_io_brd_tr_brd,var_codes_brd);
t_io_brd_tr_brd = sortrows(t_io_brd_tr_brd,"Code_nums");
t_io_brd_tr_brd.Properties.VariableNames = replace(t_io_brd_tr_brd.Properties.VariableNames,"sum_","");
t_io_brd_tr_brd.Properties.RowNames = t_io_brd_tr_brd.Code_std_20;
t_io_brd_tr_brd(:,{'Code_std_20', 'Code_nums'}) = [];

% 1.3 Final clean table
t_io_final = rows2vars(t_io_brd_tr_brd); % transpose back
t_io_final.Properties.RowNames = t_io_final.OriginalVariableNames;
t_io_final.OriginalVariableNames = [];


%% 2 Format the table into a struct to be fed into the model
fprintf(['\n  Saving the mat file for ' iso2 '...\n'])

% 3.1 Drop non-comparable and scrap goods
t_io_final({'Noncomparable', 'Scrap'},:) = [];

% 3.2 Define the first and last industries
first_name= {'Agriculture'} ; % First industry
last_name = {'Government'} ;% Last industry

idx_1 = find(strcmp(t_io_final.Properties.VariableNames, first_name)); % Index of first industry
idx_2 = find(strcmp(t_io_final.Properties.VariableNames, last_name)); % Index of last industry

% Check: Row and column names for the sectors should be the same
% display([t_io_final.Properties.RowNames(1:idx_2) t_io_final.Properties.VariableNames(1:idx_2)'])
assert(sum(~strcmp(t_io_final.Properties.RowNames(1:idx_2),t_io_final.Properties.VariableNames(1:idx_2)')) == 0)

% 3.3 Extract matrix of intermediate inputs
io_inputs = t_io_final(idx_1 : idx_2, 1 : idx_2);
industry_names = io_inputs.Properties.RowNames;

% 3.3.1 Replace nans and eliminate zeros
idx = ismissing(io_inputs);
io_inputs{:,:}(idx) = 0;
io_inputs = io_inputs{:,:} + 1;

% 3.4 Load factors of production >> separate function
idx_coe = find(strcmp('Compensation',t_io_final.Properties.RowNames));
idx_gos = find(strcmp('GoS',t_io_final.Properties.RowNames));
io_factors = t_io_final{[idx_coe, idx_gos], 1 : idx_2};

% 3.5 Extract taxes
taxes_matrix = t_io_final{strcmp('Taxes',t_io_final.Properties.RowNames), 1 : idx_2};

% 3.6 Extract final uses
% 3.6.1 Consumption
consumption = t_io_final.Consumption(idx_1 : idx_2) + 100;
total_consumption = sum(consumption);
consumption_weights = consumption ./ total_consumption;

% 3.6.2 Investment
investment = t_io_final.Investment(idx_1 : idx_2) + 100;
total_investment = sum(investment);
investment_weights = investment ./ total_investment;

% 3.6.3 Government
government_demand = t_io_final.Government_Demand(idx_1 : idx_2) + 100;
total_government_demand = sum(government_demand);
government_demand_weights = government_demand ./ total_government_demand;

% 3.6.4 Total GDP and shares
total_gdp = total_consumption + total_investment + total_government_demand;
consumption_share = total_consumption / total_gdp;
investment_share = total_investment / total_gdp;
government_share = total_government_demand / total_gdp;
final_demand = consumption + investment + government_demand;
industry_ii = sum(io_inputs,2);
gross_output = final_demand + industry_ii;
taxes_ishare = taxes_matrix' ./ gross_output;

% 3.7 Save data into a struct, io
io.iso2 = iso2;
io.industry_names = industry_names;
io.nsectors = size(industry_names,1);
io.io_inputs= [io_inputs; io_factors];

io.consumption_weights = consumption_weights;
io.investment_weights = investment_weights;
io.government_weights = government_demand_weights;

io.consumption_share= consumption_share;
io.investment_share = investment_share;
io.government_share = government_share;
io.taxes_ishare = taxes_ishare;

% 3.8 Save
save(fullfile(folder_path, 'input', [iso2 '_io_results_bea.mat']),'io');
fprintf('\n  %s_io_results_bea.mat saved under: \n  %s \n\n ', iso2, fullfile(folder_path, 'input'));
end

