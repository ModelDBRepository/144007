%
% This function hits up a given template file and fills out the parameters ...
% 
% params should contain params.idstr and params.value, with each instance of
%   <%params.idstr%> replaced by params.value
%
% generates hocfile_path_n.hoc
%
function hoc_idx = generate_hocfiles(template_path, hocfile_path, outfile_path, params, hoc_idx, ofi, n)

  if ( exist(hocfile_path) ~= 7 )
    disp(['hocfiles directory ' hocfile_path ' does not exist']);
    return;
  end;
  
  if ( exist(template_path) ~= 2 )
    disp(['unable to locate template file ' template_path]);
    return;
  end;
  
  %find the last file of that name - and assign a number 1 above it.
  existing = dir([outfile_path '/' params(ofi).value '_*']);
  match_str = [params(ofi).value '_%d.out'];
  for i=1:length(existing)
     temp = textscan(existing(i).name, match_str);
     fn(i) = temp{1};
  end
  if exist('fn', 'var')
      n = max(fn) + n;
  end
  
  
  % output fname chg
  params(ofi).value = [params(ofi).value '_' num2str(n)];
  
  % Display ...
  disp(['Generating ' hocfile_path '/' num2str(hoc_idx) '.hoc']);
  
  % File stuff
  fid = fopen(template_path, 'rt');
  content = fread(fid, 'uint8=>char')';
  fclose(fid);

  % Replace template stuff ...
  for p=1:length(params)
    content = strrep(content, ['<%' params(p).idstr '%>'], params(p).value);
  end

  % And hocfile
  fid = fopen([hocfile_path '/' num2str(hoc_idx) '.hoc'], 'w');
  fprintf(fid,'%s', content);
  fclose(fid);
  
  hoc_idx = hoc_idx  + 1;
