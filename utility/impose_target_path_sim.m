function [results] = impose_target_path_sim(i, exercise_, exercises, exercises_fields, iso2)
% Impose the targetpath
if i == 1
    for tp_ = 1 : length(targetpath1)
        eval(strcat("targetpath1(",num2str(tp_),")    = +", num2str(target_change_size(tp_)), ".*ones(1,", num2str(number_of_periods), ")")) ;%- log(y_bar(targetnum1));        % Desired change in commodity prices (in deviation from steady state)
    end
else
    if contains(eval(['exercises.', char(exercises_fields(exercise_)) '.targetvar1']),'mm')
        load(['output\temporary_shock\' iso2 '\targetpath1_eps_m.mat']);
    else
        load(['output\temporary_shock\' iso2 '\targetpath1_' eval(['exercises.', char(exercises_fields(exercise_)) '.targetvar1']), '.mat']);
    end
end
end