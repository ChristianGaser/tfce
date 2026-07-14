function spm_TFCE
% TFCE Toolbox wrapper to call TFCE functions
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

SPMid = spm('FnBanner',mfilename);
[Finter,Fgraph] = spm('FnUIsetup','TFCE1.3');
url = fullfile(fileparts(mfilename('fullpath')),'html','tfce.html');
web(url,'-noaddressbox','-new')
tfcedir = fileparts(mfilename('fullpath')); 


fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
  'Label',  'TFCE',...
  'Separator',  'on',...
  'Tag',    'TFCE',...
  'HandleVisibility','on');
h1  = uimenu(h0,...
  'Label',  'Estimate',...
  'Separator',  'off',...
  'Tag',    'Estimate TFCE statistic',...
  'CallBack','spm_jobman(''interactive'','''',''spm.tools.tfce_estimate'');',...
  'HandleVisibility','on');
h2  = uimenu(h0,...
  'Label',  'Results',...
  'Separator',  'off',...
  'Tag',    'Non-parametric inference',...
  'CallBack','[hReg xSPM SPM] = cat_spm_results_ui(''Setup'');',...
  'HandleVisibility','on');
h3  = uimenu(h0,...
  'Label',  'Check for updates',...
  'Separator',  'off',...
  'Tag',    'Check for updates',...
  'CallBack','tfce_update(1);',...
  'HandleVisibility','on');

% Proactively clear macOS Gatekeeper quarantine before any binary is called.
% We check a single representative binary for the quarantine attribute (~1ms).
% Only if found, we recursively clear quarantine from the entire TFCE directory.
if ismac && ~isdeployed
  [~, archOutput] = system('uname -v');
  if ~isempty(strfind(archOutput, 'ARM64')) %#ok<STREMP>
    mexFiles = dir(fullfile(tfcedir, '*.mexmaca64'));
  else
    mexFiles = dir(fullfile(tfcedir, '*.mexmaci64'));
  end
  if ~isempty(mexFiles)
    testMex = fullfile(tfcedir, mexFiles(1).name);
    cmd = ['xattr -p com.apple.quarantine "' testMex '"'];
    [qST, result] = system(cmd);
    if (qST == 0 && ~isempty(strtrim(result)))
      [fixStatus1, ~] = system(sprintf('xattr -dr com.apple.quarantine "%s"', tfcedir));
      if fixStatus1 ~= 0 || fixStatus2 ~= 0
        fprintf(2, '\n========================================================================\n');
        fprintf(2, 'TFCE: Could not remove macOS quarantine automatically.\n');
        fprintf(2, 'Please run this command in your Terminal to fix this:\n\n');
        fprintf(2, '     sudo xattr -dr com.apple.quarantine "%s"\n\n', tfcedir);
        fprintf(2, 'Then restart MATLAB.\n');
        fprintf(2, '========================================================================\n\n');
      end
    end
  end
end
