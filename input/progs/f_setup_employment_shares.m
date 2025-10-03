function empl = f_setup_employment_shares(folder_path, mapping_file, source, iso2, t_empl)
% -------------------------------------------------------------------------
% This function processes employment data to calculate employment shares by
% industry for a given economy.  It performs the following tasks:
%
% 1. Maps the data to broader industry classifications as defined in the
%    mapping file (\input\data\IO_mapping.xlsx).
%
% 2. Aligns the industry order with the industry classification used in
%    the corresponding input-output (IO) tables.
%
% 3. Calculates employment shares for each industry.
%
% 4. Saves the following variables into a structured format (`empl`) for
%    further use in the model:
%    a. employment_share: Employment shares by industry.
%    b. industry_names  : Names of the broader industries.
%
% 5. Saves the `empl` struct for the specified country into a `.mat` file
%    under the `\input` folder.
% -------------------------------------------------------------------------

%% 1 Load and setup data

% 1.1 Keep data for 2019
ind_year  = 'x2019';
t_empl = t_empl(:,{ind_year,'SECTORNUMBER'});
t_empl.empl = t_empl.(ind_year);

% 1.2 Keep rows containing hours worked
t_empl = t_empl(startsWith(t_empl.Properties.RowNames,'5'),:); % rows starting with 5 refers to hours worked
% '1' Domestic industry output, in millions of current dollars
% '2' Domestic industry output, in millions of chain weighted (2012) dollars
% '3' Industry output chain weighted deflator, 2012=100
% '4' Total jobs, in thousands
% '5' Total hours of all persons, in millions
% '6' Wage and salary jobs, in thousands
% '7' Wage and salary hours, in millions
% '10' Self-employed and unpaid family worker jobs, in thousands
% '11' Self-employed and unpaid family worker hours, in millions
% '12' Domestic commodity output, in millions of current dollars
% '13' Domestic commodity output, in millions of chain weighted (2012) dollars
% '14' Commodity output chain weighted deflator, 2012=100
% Please refer to the emp2022_SectorPlan30.xlsx file for the industry titles.


% 1.3 Load broader industry codes and join with the table
row_codes = readtable(mapping_file,'Sheet',['Employment_' upper(source)],'ReadVariableNames',true);
t_empl.Code = t_empl.SECTORNUMBER;
t_empl_brd = innerjoin(t_empl,row_codes,'Keys','Code');

% 1.4 Drop extra rows
t_empl_brd = t_empl_brd(startsWith(t_empl_brd.Industry,'Drop') == 0,:);

% 1.5 Aggregate according to the broader industry classification
t_empl_brd = t_empl_brd(:,{'empl','Industry'});
t_empl_brd_sum = groupsummary(t_empl_brd,'Industry','sum');

% 1.6 Get the industry name order from the relevant IO mat file and sort
% industries accordingly
io_mat_file = fullfile(folder_path, 'input',[iso2 '_io_results_' source '.mat']);

if exist(io_mat_file,'file') ~= 0

    io = load(io_mat_file);
    [~,rowidx] = ismember(io.io.industry_names,t_empl_brd_sum.Industry');
    t_empl_brd_sum_ordered = t_empl_brd_sum(rowidx,:);

    t_empl_brd_sum_ordered.Properties.RowNames = t_empl_brd_sum_ordered.Industry;
    t_empl_brd_sum_ordered = t_empl_brd_sum_ordered(:,{'sum_empl'});

    %% 2 Format the table into a struct to be fed into the model
    % 2.1 Compute shares
    employment_share = t_empl_brd_sum_ordered.sum_empl ./ sum(t_empl_brd_sum_ordered.sum_empl);

    % 2.2 Save data into a struct, empl
    empl.employment_share = employment_share;
    empl.industry_names = t_empl_brd_sum_ordered.Properties.RowNames;
    assert(sum(~strcmp(io.io.industry_names,empl.industry_names)) == 0)

    % 2.3 Save as a mat file
    save(fullfile(folder_path, 'input',[iso2 '_employment_results_' source '.mat']),'empl');
    fprintf(' %s_employment_results_%s.mat saved under: \n  %s \n\n', iso2, source, fullfile(folder_path, 'input'));

else
    empl = [];
end

end

