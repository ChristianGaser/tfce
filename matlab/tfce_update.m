function varargout = tfce_update(update)
% check for new updates
%
% FORMAT [sts, msg] = tfce_update(update)
% sts    - status code:
%        NaN - server not accessible
%        Inf - no updates available
%        0   - TFCE installation up-to-date
%        str - new version string (e.g. '1.3.1') available for download
% msg    - string describing outcome, that would otherwise be displayed.
% update - allow installation of update
% 
% This function will connect itself to the SBM server, compare the
% version number of the updates with the one of the TFCE installation 
% currently in the MATLAB path and will display the outcome.
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev$';

if isdeployed
  sts= Inf;
  msg = 'Update function is not working for compiled TFCE. Please check for a new compiled TFCE version.';
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return;
end

if ~nargin
    update = false;
else
    update = true;
end

% get current release number as a version string, e.g. '1.3.1'
A = ver;
rstr = '';
for i=1:length(A)
  if strcmp(A(i).Name,'TFCE Toolbox')
    rstr = strrep(A(i).Release,'(TFCE','');
    rstr = strrep(rstr,')','');
  end
end
r = parse_version(rstr);   % numeric vector, e.g. [1 3 1]

url_github = 'https://api.github.com/repos/ChristianGaser/tfce/releases';

% get new release numbers
try
  [jsonStr, sts] = urlread(url_github,'Timeout',2);
catch
  [jsonStr, sts] = urlread(url_github);
end

if ~sts
  msg = sprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\nYou can download your update at %s\n',url_github,url_github);
  if ~nargout, error(msg); else varargout = {NaN, msg}; end
  return
end

data = jsondecode(jsonStr);
% find the largest release version by lexicographic version comparison
rnew = [];   % parsed version vector of the newest release
ind  = 0;
for i = 1:length(data)
  v = parse_version(data(i).tag_name);
  if isempty(v), continue; end
  if isempty(rnew) || cmp_version(v, rnew) > 0
    rnew = v;
    ind  = i;
  end
end

if ind == 0
  sts = Inf;
  msg = sprintf('No valid release versions found on %s.\n',url_github);
  if ~nargout, fprintf([blanks(9) msg]); else varargout = {sts, msg}; end
  return
end
new_tag_name = data(ind).tag_name;
rnew_str = version_string(rnew);

if isempty(r) || cmp_version(rnew, r) > 0
  sts = rnew_str;
  msg = sprintf('         A new version of TFCE is available on:\n');
  msg = [msg sprintf('   %s\n',url_github)];
  msg = [msg sprintf('        (Your version: %s - New version: %s)\n',rstr,rnew_str)];
  if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
  sts = 0;
  msg = sprintf('Your version of TFCE is up to date.');
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return
end

% prefer the actual .zip asset attached to the release; fall back to the
% conventional name if the release carries no assets
url = '';
if isfield(data,'assets')
  assets = data(ind).assets;
  for k = 1:numel(assets)
    nm = assets(k).name;
    if numel(nm) >= 4 && strcmpi(nm(end-3:end),'.zip')
      url = assets(k).browser_download_url;
      break;
    end
  end
end
if isempty(url)
  url = sprintf('https://github.com/ChristianGaser/tfce/releases/download/%s/tfce_%s.zip',new_tag_name,rnew_str);
end

if update
  overwrite = spm_input('Update',1,'m','Do not update|Download zip-file only|Overwrite old TFCE installation',[-1 0 1],3);
  d0 = spm('Dir');
  d  = fileparts(fileparts(which('spm_TFCE')));
  
  if overwrite
    try
      % list mex-files and delete these files to prevent that old
      % compiled files are used
      mexfiles = dir(fullfile(fileparts(mfilename('fullpath')),'*.mex*'));
      for i=1:length(mexfiles)
        name = fullfile(fileparts(mfilename('fullpath')),mexfiles(i).name);
        spm_unlink(name);
      end

      lastwarn('');
      warning off
      delete(get(0,'Children')); spm('clean'); evalc('spm_rmpath'); drawnow
      m = '          Download and install TFCE...\n';
      if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
      s = unzip(url, d);
      m = sprintf('         Success: %d files have been updated.\n',numel(s));
      if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
      addpath(d0);
      rehash;
      rehash toolboxcache;
      if exist('toolbox_path_cache','file'), toolbox_path_cache; end
      spm fmri; spm_TFCE
      warning on

    catch
      fprintf('Update failed: check file permissions. Download zip-file only.\n');
      web(url,'-browser');
      fprintf('Unzip file to %s\n',d);
    end
  else
    web(url,'-browser');
    fprintf('Unzip file to %s\n',d);
  end
end

%==========================================================================
function v = parse_version(s)
% Convert a version string ('1.3.1', 'v1.3', 'TFCE1.3.1') into a numeric
% row vector [1 3 1]. Returns [] when no numeric component is found.
v = [];
if isempty(s) || ~ischar(s), return; end
% keep only the numeric-and-dot part, e.g. strip a leading 'v' or 'TFCE'
tok = regexp(s,'\d+(\.\d+)*','match','once');
if isempty(tok), return; end
parts = regexp(tok,'\.','split');
v = zeros(1,numel(parts));
for i = 1:numel(parts)
  v(i) = str2double(parts{i});
end

%==========================================================================
function c = cmp_version(a,b)
% Compare version vectors from parse_version. Returns -1 if a<b, 0 if
% a==b, +1 if a>b. The shorter vector is zero-padded, so 1.3 < 1.3.1.
n = max(numel(a),numel(b));
a(end+1:n) = 0;
b(end+1:n) = 0;
d = a - b;
k = find(d~=0,1);
if isempty(k), c = 0; else, c = sign(d(k)); end

%==========================================================================
function s = version_string(v)
% Join a version vector [1 3 1] back into the string '1.3.1'.
s = strjoin(arrayfun(@(x) sprintf('%g',x),v,'UniformOutput',false),'.');