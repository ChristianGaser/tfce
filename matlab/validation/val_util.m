function varargout = val_util(action, varargin)
% Helpers for the validation scripts
% FORMAT p    = val_util('extract', name1, name2, ...)
%        hdr  = val_util('header', title)
%        val_util('result', name, passed, detail)
%        ok   = val_util('summary')
%
% 'extract' pulls the named subfunctions out of tfce_estimate_stat.m into a
% temporary folder and puts it on the path. The validation therefore always
% exercises the code that ships, never a copy of it that could go stale.
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% ______________________________________________________________________
% $Id$

persistent results

switch lower(action)

  case 'extract'
    src = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'tfce_estimate_stat.m');
    if ~exist(src,'file')
      error('Cannot find %s', src);
    end
    txt = fileread(src);

    % split at the subfunction boundaries
    lines = regexp(txt, '\r?\n', 'split');
    starts = find(~cellfun(@isempty, regexp(lines, '^function\s', 'once')));

    outdir = fullfile(tempdir, 'tfce_validation_src');
    if ~exist(outdir,'dir'), mkdir(outdir); end

    for k = 1:numel(varargin)
      name = varargin{k};
      idx  = [];
      for i = 1:numel(starts)
        if ~isempty(regexp(lines{starts(i)}, ['[=\s]' name '\s*\('], 'once'))
          idx = i; break
        end
      end
      if isempty(idx)
        error('Subfunction %s not found in tfce_estimate_stat.m', name);
      end
      if idx < numel(starts)
        last = starts(idx+1) - 1;
      else
        last = numel(lines);
      end
      body = lines(starts(idx):last);
      fid  = fopen(fullfile(outdir, [name '.m']), 'w');
      fprintf(fid, '%s\n', body{:});
      fclose(fid);
    end

    addpath(outdir);
    varargout{1} = outdir;

  case 'header'
    results = [];
    fprintf('\n%s\n', repmat('=',1,74));
    fprintf('  %s\n', varargin{1});
    fprintf('%s\n', repmat('=',1,74));

  case 'result'
    name   = varargin{1};
    passed = varargin{2};
    detail = '';
    if nargin > 3, detail = varargin{3}; end
    if passed, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%s] %-44s %s\n', tag, name, detail);
    results(end+1).passed = passed; %#ok<AGROW>

  case 'summary'
    n  = numel(results);
    np = sum([results.passed]);
    fprintf('\n  %d of %d checks passed\n', np, n);
    if np < n
      fprintf('  *** %d FAILED ***\n', n - np);
    end
    varargout{1} = (np == n);

end
