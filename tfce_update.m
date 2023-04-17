function varargout = tfce_update(update)
% check for new updates
%
% FORMAT [sts, msg] = tfce_update(update)
% sts    - status code:
%        NaN - server not accessible
%        Inf - no updates available
%        0   - TFCE installation up-to-date
%        n   - new revision <n> is available for download
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

rev = '$Rev: 189 $';

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

r = 0;

% get current release number
A = ver;
for i=1:length(A)
  if strcmp(A(i).Name,'TFCE Toolbox')
    r = str2double(A(i).Version);
  end
end

url = 'http://www.neuro.uni-jena.de/tfce/';

% get new release numbers
try
  [s,sts] = urlread(url,'Timeout',2);
catch
  [s,sts] = urlread(url);
end

if ~sts
  sts = NaN;
  msg = sprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\nYou can download your update at %s\n',url,url); 
  if ~nargout, error(msg); else varargout = {sts, msg}; end
  return
end

n = regexp(s,'tfce_r(\d.*?)\.zip','tokens');
if isempty(n)
  sts= Inf;
  msg = 'There are no updates available yet.';
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return;
else
  % get largest release number
  rnew = [];
  for i=1:length(n)
    rnew = [rnew str2double(n{i})];
  end
  rnew = max(rnew);
end

if rnew > r
  sts = n;
  msg = sprintf('         A new version of TFCE is available on:\n');
  msg = [msg sprintf('   %s\n',url)];
  msg = [msg sprintf('        (Your version: %d - New version: %d)\n',r,rnew)];
  if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
  sts = 0;
  msg = sprintf('Your version of TFCE is up to date.');
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return
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
      s = unzip([url sprintf('tfce_r%d.zip',rnew)], d);
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
      web([url sprintf('tfce_r%d.zip',rnew)],'-browser');
      fprintf('Unzip file to %s\n',d);
    end
  else
    web([url sprintf('tfce_r%d.zip',rnew)],'-browser');
    fprintf('Unzip file to %s\n',d);
  end
end