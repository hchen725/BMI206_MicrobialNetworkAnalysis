% Get the cascade number for each species reactions
%% Directories
main_dir = cd;
adj_mat_dir = strrep(main_dir, ...
    "scripts", ...
    "data/AdjacencyMatrixCSV/");
msg = strcat("Reading CSV files from ", adj_mat_dir);
disp(msg);

output_dir = strrep(adj_mat_dir, ...
    "AdjacencyMatrixCSV", ...
    "CascadeNumbers");
mkdir(output_dir);

msg = strcat("Saving cascade number files to: ", output_dir);
disp(msg);

%% Calculate cascade numbers
adj_mat_files = dir(fullfile(adj_mat_dir, "*.csv"));

for k = 1:length(adj_mat_files)
    fname = adj_mat_files(k).name;
    full_fname = fullfile(adj_mat_dir, fname);
    fprintf(1, "Now reading %s\n", fname);

    current_file = readtable(full_fname);
    adj_mat = table2array(current_file);

    % Initailize array to store cascade numbers
    sz = size(current_file, 1);
    list_of_cascade_number=zeros(1, sz);

    for i=1:sz

	    affected_nodes = cascade(adj_mat, i);
	    list_of_cascade_number(1,i)=length(affected_nodes);
	
    end

    outname = strrep(fname, "AdjMat", "CascadeNum");
    writematrix(list_of_cascade_number, ...
        strcat(output_dir, outname))
end

disp("Done calculating cascade numbers.")