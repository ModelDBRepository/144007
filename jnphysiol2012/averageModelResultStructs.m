function mdl_avg = averageModelResultStructs(model_files)
% function mdl_avg = averageModelResults(model_files)
% Need to open a set of looming simulation results and average them
% The input model_files is a cell array of strings.  The function
% loads them in sequence and averages the fields.


% I looked into doing this a very general way, but it seemed a little difficult.  Plus, 
% these structures are actually really big, since there are many traces saved.  So, the 
% strategy is to cut out the 'trials' field, and array the rest of the 'stim' struct

stim_array = [];
for i = 1:length(model_files)
    load(model_files{i}); %this should load a 'stim' and 'stim_str' variable, 'stim' is a 1xN condition struct
    stim = rmfield(stim, 'trial'); %this field contains trial info and takes up much space.  Deleting for summary data.
    fn = fieldnames(stim);
    for j=1:length(fn)
        a = textscan(fn{j}, 'a_%s');
        if (~isempty(a{1}))
            is_matrix(j) = 1;
        else
            is_matrix(j) = 0;
        end
    end
    stim = rmfield(stim, {fn{logical(is_matrix)}});
    stim_array = cat(1, stim_array, stim); %cat them keeping the l/v values (trial types) separate
end

% Generic way to average the numeric fields within a structure
as = size(stim_array);
fn = fieldnames(stim_array);
mdl_avg = struct([]);
for i = 1:length(fn)
    if isnumeric([stim_array.(fn{i})])
        temp = zeros([as(1) as(2) size(stim_array(1).(fn{i}))]) .* NaN; 
        for j=1:as(1)
            for k=1:as(2)
                temp(j,k,:,:) = stim_array(j,k).(fn{i});
            end
        end
        mdl_avg(1).(fn{i}) = mean(temp);
    end
end
    



    
