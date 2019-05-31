%
% Generates .sh files to run on n_cpus in batches
% params - n_cpus, hoc_start, hoc_end, hoc_rootdir, special_dir
% If n_cpus is a vector, then it means that total CPUs should be split up
% over length(n_cpus) batches, and special_dir should be a cell array of strings indicating
% the special executable directory for each machine

function generate_sh_files(n_cpus, hoc_start, hoc_end, hoc_rootdir, special_dir)
  
  % check that the special directories given actually exist
  if iscell(special_dir)
      for i = 1:length(special_dir)
          ss = char(special_dir{i});
          if ( exist([hoc_rootdir '/' ss '/special']) ~= 2 )
              disp(['unable to locate ' ss '/special']);
              return;
          end;
      end
  else
      if ( exist([hoc_rootdir '/' special_dir '/special']) ~= 2 )
          disp(['unable to locate ' special_dir '/special']);
          return
      end    
  end
  
  n_batches = length(n_cpus);
  total_cpus = sum(n_cpus);
  nhoc = hoc_end-hoc_start+1;
  rem = mod(nhoc,total_cpus);
  nhoc_cpu = [ones(1, rem), zeros(1, total_cpus-rem)]; %the number of hocfiles per cpu - distribute with the remainder then add the even
  nhoc_cpu = nhoc_cpu + floor(nhoc/total_cpus)*ones(1,total_cpus);
 
  n_prev_cpu = 0; hi = hoc_start;
  for j=1:n_batches
      cpu_range = (n_prev_cpu + 1):(n_cpus(j) + n_prev_cpu); %the cpu numbers in this batch
      
      % first, write the batch file
      fid = fopen([hoc_rootdir '/runbatch' num2str(j) '.sh'], 'w');
      for n=cpu_range
          fprintf(fid,'./run_%d.sh > /dev/null &\n', n);
      end
      fclose(fid);
      n_prev_cpu = cpu_range(end);
  
      % need to check if the special_dir is a cell array or just a string
      if iscell(special_dir) special_str = char(special_dir{j});
      else special_str = special_dir; end
      
      % now write the individual per cpu files that list the hoc files
      for n=cpu_range
          fid = fopen([hoc_rootdir '/run_' num2str(n) '.sh'], 'w');
          for h=1:nhoc_cpu(n)
              fprintf(fid, '%s/special %d.hoc\n', special_str, hi);
              hi = hi + 1;
          end
          fclose(fid);
      end
  end
  