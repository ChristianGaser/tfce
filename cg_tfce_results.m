function varargout = cg_tfce_results(varargin)
% User interface for TFCE results: Display and analysis of regional effects
% FORMAT [hReg,xSPM,SPM] = cg_tfce_results('Setup',[xSPM])
%
% hReg   - handle of MIP XYZ registry object
%          (see spm_XYZreg.m for details)
% xSPM   - structure containing specific SPM, distribution & filtering details
%          (see spm_getSPM.m for contents)
% SPM    - SPM structure containing generic parameters
%          (see spm_spm.m for contents)
%_______________________________________________________________________
% Christian Gaser
% $Id$

% based on spm_results_ui.m

%-Condition arguments
%--------------------------------------------------------------------------
if nargin == 0, Action='Setup'; else Action=varargin{1}; end
 
 
%==========================================================================
switch lower(Action), case 'setup'                         %-Set up results
%==========================================================================
        
    %-Initialise
    %----------------------------------------------------------------------
    spm('FnBanner',mfilename);
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','TFCE: Results');
    spm_clf('Satellite');
  
    %-Get thresholded xSPM data and parameters of design
    %======================================================================
    if nargin > 1
        [SPM,xSPM] = cg_get_tfce_results(varargin{2});
    else
        [SPM,xSPM] = cg_get_tfce_results;
    end
 
    if isempty(xSPM) 
        varargout = {[],[],[]};
        return;
    end
 
    %-Check whether mesh are detected if we use spm12
    %--------------------------------------------------------------------------
    if strcmp(spm('ver'),'SPM12')
        if spm_mesh_detect(SPM.xY.VY)
            mesh_detected = 1;
        else
            mesh_detected = 0;
        end
    else
          mesh_detected = 0;
    end

    %-Ensure pwd = swd so that relative filenames are valid
    %----------------------------------------------------------------------
    cd(SPM.swd)
    
    %-Get space information
    %======================================================================
    M         = SPM.xVol.M;
    DIM       = SPM.xVol.DIM;

    %-Space units
    %----------------------------------------------------------------------
    try
        try
            units = SPM.xVol.units;
        catch
            units = xSPM.units;
        end
    catch
        try
            Modality = spm('CheckModality');
        catch
            Modality = {'PET','FMRI','EEG'};
            selected = spm_input('Modality: ','+1','m',Modality);
            Modality = Modality{selected};
            spm('ChMod',Modality);
        end
        if strcmp(Modality,'EEG')
            datatype = {...
                'Volumetric (2D/3D)',...
                'Scalp-Time',...
                'Scalp-Frequency',...
                'Time-Frequency',...
                'Frequency-Frequency'};
            selected = spm_input('Data Type: ','+1','m',datatype);
            datatype = datatype{selected};
        else
            datatype = 'Volumetric (2D/3D)';
        end
        
        switch datatype
            case 'Volumetric (2D/3D)'
                units    = {'mm' 'mm' 'mm'};
            case 'Scalp-Time'
                units    = {'mm' 'mm' 'ms'};
            case 'Scalp-Frequency'
                units    = {'mm' 'mm' 'Hz'};
            case 'Time-Frequency'
                units    = {'Hz' 'ms' ''};
            case 'Frequency-Frequency'
                units    = {'Hz' 'Hz' ''};
            otherwise
                error('Unknown data type.');
        end
    end
    if mesh_detected
        DIM(3) = Inf; % force 3D coordinates
    elseif DIM(3) == 1
        units{3} = '';
        if DIM(2) == 1
            units{2} = '';
        end
    end
    xSPM.units      = units;
    SPM.xVol.units  = units;
    
     
    %-Setup Results User Interface; Display MIP, design matrix & parameters
    %======================================================================
 
    %-Setup results GUI
    %----------------------------------------------------------------------
    spm_clf(Finter);
    spm('FigName',['SPM{',xSPM.STAT,'}: Results'],Finter,CmdLine);
    hReg      = cg_tfce_results('SetupGUI',M,DIM,xSPM,Finter);
 
    %-Setup design interrogation menu
    %----------------------------------------------------------------------
    hDesRepUI = spm_DesRep('DesRepUI',SPM);
 
    %-Atlas menu
    %----------------------------------------------------------------------
    if isequal(units,{'mm' 'mm' 'mm'}) & strcmp(spm('ver'),'SPM12')
        hAtlasUI = cg_tfce_results('SetupAtlasMenu',Finter);
    end

    %-Setup Maximum intensity projection (MIP) & register
    %----------------------------------------------------------------------
    FS     = spm('FontSizes');
    hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
    if mesh_detected
        hMax = spm_mesh_render('Disp',SPM.xVol.G,'Parent',hMIPax);
        tmp = zeros(1,prod(xSPM.DIM));
        tmp(xSPM.XYZ(1,:)) = xSPM.Z;
        hMax = spm_mesh_render('Overlay',hMax,tmp);
        hMax = spm_mesh_render('Register',hMax,hReg);
    elseif isequal(units(2:3),{'' ''})
        set(hMIPax, 'Position',[0.05 0.65 0.55 0.25]);
        [allS,allXYZmm] = spm_read_vols(xSPM.Vspm);
        plot(hMIPax,allXYZmm(1,:),allS,'Color',[0.6 0.6 0.6]);
        set(hMIPax,'NextPlot','add');
        MIP = NaN(1,xSPM.DIM(1));
        MIP(xSPM.XYZ(1,:)) = xSPM.Z;
        XYZmm = xSPM.M(1,:)*[1:xSPM.DIM(1);zeros(2,xSPM.DIM(1));ones(1,xSPM.DIM(1))];
        plot(hMIPax,XYZmm,MIP,'b-+','LineWidth',2);
        plot(hMIPax,[XYZmm(1) XYZmm(end)],[xSPM.u xSPM.u],'r');
        clim = get(hMIPax,'YLim');
        axis(hMIPax,[sort([XYZmm(1) XYZmm(end)]) 0 clim(2)]);
        %set(hMIPax,'XTick',[],'YTick',[]);
    else
        hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax,units);
        spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
    end
    
    if xSPM.STAT == 'P'
        str = xSPM.STATstr;
    else
        str = ['SPM\{',xSPM.STATstr,'\}'];
    end
    text(240,260,str,...
        'Interpreter','TeX',...
        'FontSize',FS(14),'Fontweight','Bold',...
        'Parent',hMIPax)
 
 
    %-Print comparison title
    %----------------------------------------------------------------------
    hTitAx = axes('Parent',Fgraph,...
        'Position',[0.02 0.96 0.96 0.04],...
        'Visible','off');
 
    text(0.5,0.5,xSPM.title,'Parent',hTitAx,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'FontWeight','Bold','FontSize',FS(14))
 
 
    %-Print SPMresults: Results directory & thresholding info
    %----------------------------------------------------------------------
    hResAx = axes('Parent',Fgraph,...
        'Position',[0.05 0.55 0.45 0.05],...
        'DefaultTextVerticalAlignment','baseline',...
        'DefaultTextFontSize',FS(9),...
        'DefaultTextColor',[1,1,1]*.7,...
        'Units','points',...
        'Visible','off');
    AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
    h     = text(0,24,'SPMresults:','Parent',hResAx,...
        'FontWeight','Bold','FontSize',FS(14));
    text(get(h,'Extent')*[0;0;1;0],24,spm_file(SPM.swd,'short30'),'Parent',hResAx)
    try
        thresDesc = xSPM.thresDesc;
        text(0,12,sprintf('Threshold %s',thresDesc),'Parent',hResAx)
    catch
        text(0,12,sprintf('Height threshold %c = %0.6f',xSPM.STAT,xSPM.u),'Parent',hResAx)
    end
    if mesh_detected, str = 'vertices'; else str = 'voxels'; end
    if xSPM.STAT == 'T', text(0,00,sprintf('Extent threshold k = %0.0f %s',xSPM.k,str), 'Parent',hResAx); end
 
 
    %-Plot design matrix
    %----------------------------------------------------------------------
    hDesMtx   = axes('Parent',Fgraph,'Position',[0.65 0.55 0.25 0.25]);
    hDesMtxIm = image((SPM.xX.nKX + 1)*32,'Parent',hDesMtx);
    xlabel(hDesMtx,'Design matrix','FontSize',FS(10))
    set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')',...
        'UserData',struct(...
        'X',        SPM.xX.xKXs.X,...
        'fnames',   {reshape({SPM.xY.VY.fname},size(SPM.xY.VY))},...
        'Xnames',   {SPM.xX.name}))
 
    %-Plot contrasts
    %----------------------------------------------------------------------
    nPar   = size(SPM.xX.X,2);
    xx     = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
    nCon   = length(xSPM.Ic);
    xCon   = SPM.xCon;
    if nCon
        dy     = 0.15/max(nCon,2);
        hConAx = axes('Parent',Fgraph, 'Position',[0.65 (0.80 + dy*.1) 0.25 dy*(nCon-.1)],...
            'Tag','ConGrphAx','Visible','off');
        if xSPM.invResult
          str    = 'inverse contrast';
        else
          str    = 'contrast';
        end
        if nCon > 1, str = [str 's']; end
        title(hConAx,str)
        htxt   = get(hConAx,'title');
        set(htxt,'FontSize',FS(10),'FontWeight','normal','Visible','on','HandleVisibility','on')
    end
 
    for ii = nCon:-1:1
        hCon = axes('Parent',Fgraph, 'Position',[0.65 (0.80 + dy*(nCon - ii +.1)) 0.25 dy*.9]);
        if xCon(xSPM.Ic(ii)).STAT == 'T' && size(xCon(xSPM.Ic(ii)).c,2) == 1
 
            %-Single vector contrast for SPM{t} - bar
            %--------------------------------------------------------------
            yy = [zeros(1,nPar);repmat(xCon(xSPM.Ic(ii)).c',2,1);zeros(1,nPar)];
            if xSPM.invResult, yy = -yy; end
            h  = patch(xx,yy,[1,1,1]*.5,'Parent',hCon);
            set(hCon,'Tag','ConGrphAx',...
                'Box','off','TickDir','out',...
                'XTick',spm_DesRep('ScanTick',nPar,10) - 0.5,'XTickLabel','',...
                'XLim', [0,nPar],...
                'YTick',[-1,0,+1],'YTickLabel','',...
                'YLim',[min(xCon(xSPM.Ic(ii)).c),max(xCon(xSPM.Ic(ii)).c)] +...
                [-1 +1] * max(abs(xCon(xSPM.Ic(ii)).c))/10  )
 
        else
 
            %-F-contrast - image
            %--------------------------------------------------------------
            h = image((xCon(xSPM.Ic(ii)).c'/max(abs(xCon(xSPM.Ic(ii)).c(:)))+1)*32,...
                'Parent',hCon);
            set(hCon,'Tag','ConGrphAx',...
                'Box','on','TickDir','out',...
                'XTick',spm_DesRep('ScanTick',nPar,10),'XTickLabel','',...
                'XLim', [0,nPar]+0.5,...
                'YTick',[0:size(SPM.xCon(xSPM.Ic(ii)).c,2)]+0.5,...
                'YTickLabel','',...
                'YLim', [0,size(xCon(xSPM.Ic(ii)).c,2)]+0.5 )
 
        end
        ylabel(hCon,num2str(xSPM.Ic(ii)),'FontSize',FS(10),'FontWeight','normal')
        set(h,'ButtonDownFcn','spm_DesRep(''SurfCon_CB'')',...
            'UserData', struct( 'i',    xSPM.Ic(ii),...
                                'h',    htxt,...
                                'xCon', xCon(xSPM.Ic(ii))))
    end
 
 
    %-Store handles of results section Graphics window objects
    %----------------------------------------------------------------------
    H  = get(Fgraph,'Children');
    H  = findobj(H,'flat','HandleVisibility','on');
    H  = findobj(H);
    Hv = get(H,'Visible');
    set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))
 
 
    %-Finished results setup
    %----------------------------------------------------------------------
    varargout = {hReg,xSPM,SPM};
    spm('Pointer','Arrow')

    %======================================================================
    case 'setupgui'                            %-Set up results section GUI
    %======================================================================
        % hReg = cg_tfce_results('SetupGUI',M,DIM,xSPM,Finter)
        if nargin < 5, Finter='Interactive'; else, Finter = varargin{5}; end
        if nargin < 4, error('Insufficient arguments'), end
        M      = varargin{2};
        DIM    = varargin{3};
        Finter = spm_figure('GetWin',Finter);
        WS     = spm('WinScale');
        FS     = spm('FontSizes');
 
        %-Create frame for Results GUI objects
        %------------------------------------------------------------------
        hPan = uipanel('Parent',Finter,'Title','','Units','Pixels',...
                'Position',[001 001 400 190].*WS,...
                'BorderType','Line', 'HighlightColor',[0 0 0],...
                'BackgroundColor',spm('Colour'));
        hReg = uipanel('Parent',hPan,'Title','','Units','Pixels',...
                'BorderType','Etchedin', ...
                'Position',[005 005 390 180].*WS,...
                'BackgroundColor',[179 179 179]/255);

        %-Initialise registry in hReg frame object
        %------------------------------------------------------------------
        [hReg,xyz] = spm_XYZreg('InitReg',hReg,M,DIM,[0;0;0]);
 
        %-Setup editable XYZ widgets & cross register with registry
        %------------------------------------------------------------------
        hFxyz      = cg_tfce_results('DrawXYZgui',M,DIM,varargin{4},xyz,hReg);
        spm_XYZreg('XReg',hReg,hFxyz,'cg_tfce_results');
 
        %-Set up buttons for results functions
        %------------------------------------------------------------------
        cg_tfce_results('DrawButts',hReg,DIM,Finter,WS,FS);
 
        if spm_check_version('matlab','7.11') ~= 0
            drawnow; % required to force "ratio locking"
            set(findobj(hPan),'Units','Normalized','FontUnits','Normalized');
        end

        varargout  = {hReg};
 
 
 
    %======================================================================
    case 'drawbutts'   %-Draw results section buttons in Interactive window
    %======================================================================
        % cg_tfce_results('DrawButts',hReg,DIM,Finter,WS,FS)
        %
        if nargin<3, error('Insufficient arguments'), end
        hReg = varargin{2};
        DIM  = varargin{3};
        if nargin<4,  Finter = spm_figure('FindWin','Interactive');
        else, Finter = varargin{4}; end
        if nargin < 5, WS = spm('WinScale');    else,   WS = varargin{5}; end
        if nargin < 6, FS = spm('FontSizes');   else,   FS = varargin{6}; end
 
        %-p-values
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','p-values','Units','Pixels',...
            'Position',[005 085 110 092].*WS,...
            'BorderType','Beveledout', ...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        uicontrol('Parent',hPan,'Style','PushButton','String','whole brain',...
            'Units','Pixels',...
            'FontSize',FS(10),...
            'ToolTipString',...
            'Tabulate summary of local maxima, p-values & statistics',...
            'Callback','TabDat = cg_tfce_list(''List'',xSPM,hReg);',...
            'Interruptible','on','Enable','on',...
            'Position',[005 055 100 020].*WS);
        uicontrol('Parent',hPan,'Style','PushButton','String','current cluster',...
            'Units','Pixels',...
            'FontSize',FS(10),...
            'ToolTipString',...
            'Tabulate p-values & statistics for local maxima of nearest cluster',...
            'Callback','TabDat = cg_tfce_list(''ListCluster'',xSPM,hReg);',...
            'Interruptible','on','Enable','on',...
            'Position',[005 030 100 020].*WS);
 
        %-SPM area - used for Volume of Interest analyses
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','Multivariate','Units','Pixels',...
            'Position',[120 085 150 092].*WS,...
            'BorderType','Beveledout', ...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        uicontrol('Parent',hPan,'Style','PushButton','String','eigenvariate',...
            'Position',[005 055 069 020].*WS,...
            'ToolTipString',...
            'Responses (principal eigenvariate) in volume of interest',...
            'Callback','[Y,xY] = spm_regions(xSPM,SPM,hReg)',...
            'Interruptible','on','Enable','on',...
            'FontSize',FS(10));
          
        %-Visualisation
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','Display','Units','Pixels',...
            'Position',[275 085 110 092].*WS,...
            'BorderType','Beveledout',...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        uicontrol('Parent',hPan,'Style','PushButton','String','plot',...
            'FontSize',FS(10),...
            'ToolTipString','plot data & contrasts at current voxel',...
            'Callback','[Y,y,beta,Bcov] = spm_graph_ui(xSPM,SPM,hReg);',...
            'Interruptible','on','Enable','on',...
            'Position',[005 055 100 020].*WS,...
            'Tag','plotButton');
 
        str  = { 'overlays...','slices','sections','montage','render','previous sections','previous render'};
        tstr = { 'overlay filtered SPM on another image: ',...
            '3 slices / ','slice overlay /','ortho sections / ','render /','previous ortho sections /','previous surface rendering'};
 
        tmp  = { 'spm_transverse(''set'',xSPM,hReg)',...
            'spm_sections(xSPM,hReg)',...
            {@myslover},...
            ['spm_render(   struct( ''XYZ'',    xSPM.XYZ,',...
            '''t'',     xSPM.Z'',',...
            '''mat'',   xSPM.M,',...
            '''dim'',   xSPM.DIM))'],...
            ['global prevsect;','spm_sections(xSPM,hReg,prevsect)'],...
            ['global prevrend;','if ~isstruct(prevrend)',...
            'prevrend = struct(''rendfile'','''',''brt'',[],''col'',[]); end;',...            
            'spm_render(    struct( ''XYZ'',    xSPM.XYZ,',...
            '''t'',     xSPM.Z'',',...
            '''mat'',   xSPM.M,',...
            '''dim'',   xSPM.DIM),prevrend.brt,prevrend.rendfile)']};
 
        uicontrol('Parent',hPan,'Style','popupmenu','String',str,...
            'FontSize',FS(10),...
            'ToolTipString',cat(2,tstr{:}),...
            'Callback','spm(''PopUpCB'',gcbo)',...
            'UserData',tmp,...
            'Interruptible','on','Enable','on',...
            'Position',[005 030 100 020].*WS);
 
        str = {'save...',...
               'thresholded SPM',...
               'all clusters (binary)',...
               'all clusters (n-ary)',...
               'current cluster'};
        tmp = {{@mysavespm, 'thresh' },...
               {@mysavespm, 'binary' },...
               {@mysavespm, 'n-ary'  },...
               {@mysavespm, 'current'}};
        
        uicontrol('Parent',hPan,'Style','popupmenu','String',str,...
            'FontSize',FS(10),...
            'ToolTipString','Save as image',...
            'Callback','spm(''PopUpCB'',gcbo)',...
            'UserData',tmp,...
            'Interruptible','on','Enable','on',...
            'Position',[005 005 100 020].*WS);
 
        %-ResultsUI controls
        %------------------------------------------------------------------
        uicontrol('Parent',hReg,'Style','PushButton','String','clear',...
            'ToolTipString','Clear results subpane',...
            'FontSize',FS(9),'ForegroundColor','b',...
            'Callback',['cg_tfce_results(''Clear''); ',...
              'spm_input(''!DeleteInputObj''),',...
              'spm_clf(''Satellite'')'],...
            'Interruptible','on','Enable','on',...
            'DeleteFcn','spm_clf(''Graphics'')',...
            'Position',[280 050 048 020].*WS);
 
        uicontrol('Parent',hReg,'Style','PushButton','String','exit',...
            'ToolTipString','Exit the results section',...
            'FontSize',FS(9),'ForegroundColor','r',...
            'Callback','cg_tfce_results(''close'')',...
            'Interruptible','on','Enable','on',...
            'Position',[332 050 048 020].*WS);

    %======================================================================
    case 'setupatlasmenu'                                %-Setup Atlas Menu
    %======================================================================
        % cg_tfce_results('SetupAtlasMenu',Finter)
    
        Finter = varargin{2};
        
        hC   = uimenu(Finter,'Label','Atlas', 'Tag','AtlasUI');
        
        hC1  = uimenu(hC,'Label','Label using');
        
        list = spm_atlas('List','installed');
        for i=1:numel(list)
            uimenu(hC1,'Label',list(i).name,...
                'Callback',sprintf('cg_tfce_list(''label'',''%s'');',list(i).name));
        end
        if isempty(list), set(hC1,'Enable','off'); end
                
        varargout = {hC};
    
 
    %======================================================================
    case 'drawxyzgui'                                   %-Draw XYZ GUI area
    %======================================================================
        % hFxyz = cg_tfce_results('DrawXYZgui',M,DIM,xSPM,xyz,Finter)
        if nargin<6,  hReg=spm_XYZreg('FindReg','Interactive');
        else hReg=varargin{6}; end
        if nargin < 5, xyz=[0;0;0]; else xyz=varargin{5}; end
        if nargin < 4, error('Insufficient arguments'), end
        DIM     = varargin{3};
        M       = varargin{2};
        xyz     = spm_XYZreg('RoundCoords',xyz,M,DIM);
 
        %-Locate windows etc...
        %------------------------------------------------------------------
        WS      = spm('WinScale');
        FS      = spm('FontSizes');
        PF      = spm_platform('fonts');
 
        %-Create XYZ control objects
        %------------------------------------------------------------------
       hFxyz = uipanel('Parent',hReg,'Title','co-ordinates','Units','Pixels',...
            'Position',[005 005 265 040].*WS,...
            'BorderType','Beveledout',...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
 
        uicontrol('Parent',hReg,'Style','Text','String','x =',...
            'Position',[015 010 024 018].*WS,...
            'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
            'HorizontalAlignment','Center');
        hX = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(1)),...
            'ToolTipString','enter x-coordinate',...
            'Position',[039 010 056 020].*WS,...
            'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
            'HorizontalAlignment','Right',...
            'Tag','hX',...
            'Callback','cg_tfce_results(''EdWidCB'')');
 
        uicontrol('Parent',hReg,'Style','Text','String','y =',...
            'Position',[100 010 024 018].*WS,...
            'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
            'HorizontalAlignment','Center');
        hY = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(2)),...
            'ToolTipString','enter y-coordinate',...
            'Position',[124 010 056 020].*WS,...
            'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
            'HorizontalAlignment','Right',...
            'Tag','hY',...
            'Callback','cg_tfce_results(''EdWidCB'')');
 
        if DIM(3) ~= 1
        uicontrol('Parent',hReg,'Style','Text','String','z =',...
            'Position',[185 010 024 018].*WS,...
            'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
            'HorizontalAlignment','Center');
        hZ = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(3)),...
            'ToolTipString','enter z-coordinate',...
            'Position',[209 010 056 020].*WS,...
            'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
            'HorizontalAlignment','Right',...
            'Tag','hZ',...
            'Callback','cg_tfce_results(''EdWidCB'')');
        else
        hZ = [];
        end
        
        %-Statistic value reporting pane
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','statistic','Units','Pixels',...
            'Position',[275 005 110 040].*WS,...
            'BorderType','Beveledout', ...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        hSPM = uicontrol('Parent',hPan,'Style','Text','String','',...
            'Position',[005 001 100 020].*WS,...
            'FontSize',FS(10),...
            'HorizontalAlignment','Center');

 
        %-Store data
        %------------------------------------------------------------------
        set(hFxyz,'Tag','hFxyz','UserData',struct(...
            'hReg', [],...
            'M',    M,...
            'DIM',  DIM,...
            'XYZ',  varargin{4}.XYZmm,...
            'Z',    varargin{4}.Z,...
            'hX',   hX,...
            'hY',   hY,...
            'hZ',   hZ,...
            'hSPM', hSPM,...
            'xyz',  xyz ));
 
        set([hX,hY,hZ],'UserData',hFxyz)
        varargout = {hFxyz};
 
 
    %======================================================================
    case 'edwidcb'                          %-Callback for editable widgets
    %======================================================================
        % cg_tfce_results('EdWidCB')
 
        hC    = gcbo;
        d     = find(strcmp(get(hC,'Tag'),{'hX','hY','hZ'}));
        hFxyz = get(hC,'UserData');
        UD    = get(hFxyz,'UserData');
        xyz   = UD.xyz;
        nxyz  = xyz;
 
        o = evalin('base',['[',get(hC,'String'),']'],'sprintf(''error'')');
        if ischar(o) || length(o)>1
            warning(sprintf('%s: Error evaluating ordinate:\n\t%s',...
                mfilename,lasterr))
        else
            nxyz(d) = o;
            nxyz = spm_XYZreg('RoundCoords',nxyz,UD.M,UD.DIM);
        end
 
        if abs(xyz(d)-nxyz(d))>0
            UD.xyz = nxyz; set(hFxyz,'UserData',UD)
            if ~isempty(UD.hReg), spm_XYZreg('SetCoords',nxyz,UD.hReg,hFxyz); end
            set(hC,'String',sprintf('%.3f',nxyz(d)))
            cg_tfce_results('UpdateSPMval',UD)
        end
 
    %======================================================================
    case 'updatespmval'                           %-Update SPM value in GUI
    %======================================================================
        % cg_tfce_results('UpdateSPMval',hFxyz)
        % cg_tfce_results('UpdateSPMval',UD)
        if nargin<2, error('insufficient arguments'), end
        if isstruct(varargin{2}), UD=varargin{2}; else, UD = get(varargin{2},'UserData'); end
        i  = spm_XYZreg('FindXYZ',UD.xyz,UD.XYZ);
        if isempty(i), str = ''; else, str = sprintf('%6.2f',UD.Z(i)); end
        set(UD.hSPM,'String',str);
 
 
    %======================================================================
    case 'getcoords'             % Get current co-ordinates from XYZ widget
    %======================================================================
        % xyz = cg_tfce_results('GetCoords',hFxyz)
        if nargin<2, hFxyz='Interactive'; else, hFxyz=varargin{2}; end
        hFxyz     = cg_tfce_results('FindXYZframe',hFxyz);
        varargout = {getfield(get(hFxyz,'UserData'),'xyz')};
 
 
 
    %======================================================================
    case 'setcoords'                       % Set co-ordinates to XYZ widget
    %======================================================================
        % [xyz,d] = cg_tfce_results('SetCoords',xyz,hFxyz,hC)
        if nargin<4, hC=0; else, hC=varargin{4}; end
        if nargin<3, hFxyz=cg_tfce_results('FindXYZframe'); else, hFxyz=varargin{3}; end
        if nargin<2, error('Set co-ords to what!'), else, xyz=varargin{2}; end
 
        %-If this is an internal call, then don't do anything
        if hFxyz==hC, return, end
 
        UD = get(hFxyz,'UserData');
 
        %-Check validity of coords only when called without a caller handle
        %------------------------------------------------------------------
        if hC <= 0
            [xyz,d] = spm_XYZreg('RoundCoords',xyz,UD.M,UD.DIM);
            if d>0 && nargout<2, warning(sprintf(...
                '%s: Co-ords rounded to nearest voxel centre: Discrepancy %.2f',...
                mfilename,d))
            end
        else
            d = [];
        end
 
        %-Update xyz information & widget strings
        %------------------------------------------------------------------
        UD.xyz = xyz; set(hFxyz,'UserData',UD)
        set(UD.hX,'String',sprintf('%.2f',xyz(1)))
        set(UD.hY,'String',sprintf('%.2f',xyz(2)))
        set(UD.hZ,'String',sprintf('%.2f',xyz(3)))
        cg_tfce_results('UpdateSPMval',UD)
 
        %-Tell the registry, if we've not been called by the registry...
        %------------------------------------------------------------------
        if (~isempty(UD.hReg) && UD.hReg~=hC)
            spm_XYZreg('SetCoords',xyz,UD.hReg,hFxyz);
        end
 
        %-Return arguments
        %------------------------------------------------------------------
        varargout = {xyz,d};
 
 
 
    %======================================================================
    case 'findxyzframe'                                  % Find hFxyz frame
    %======================================================================
        % hFxyz = cg_tfce_results('FindXYZframe',h)
        % Sorts out hFxyz handles
        if nargin<2, h='Interactive'; else, h=varargin{2}; end
        if ischar(h), h=spm_figure('FindWin',h); end
        if ~ishandle(h), error('invalid handle'), end
        if ~strcmp(get(h,'Tag'),'hFxyz'), h=findobj(h,'Tag','hFxyz'); end
        if isempty(h), error('XYZ frame not found'), end
        if length(h)>1, error('Multiple XYZ frames found'), end
        varargout = {h};
 
 
 
    %======================================================================
    case 'plotui'                               %-GUI for plot manipulation
    %======================================================================
        % cg_tfce_results('PlotUi',hAx)
        if nargin<2, hAx=gca; else, hAx=varargin{2}; end
 
        WS = spm('WinScale');
        FS = spm('FontSizes');
        Finter=spm_figure('FindWin','Interactive');
        figure(Finter)
 
        %-Check there aren't already controls!
        %------------------------------------------------------------------
        hGraphUI = findobj(Finter,'Tag','hGraphUI');
        if ~isempty(hGraphUI)           %-Controls exist
            hBs = get(hGraphUI,'UserData');
            if hAx==get(hBs(1),'UserData')  %-Controls linked to these axes
                return
            else                %-Old controls remain
                delete(findobj(Finter,'Tag','hGraphUIbg'))
            end
        end
 
         %-Frames & text
        %------------------------------------------------------------------
        hGraphUIbg = uipanel('Parent',Finter,'Title','','Tag','hGraphUIbg',...
            'BackgroundColor',spm('Colour'),...
            'BorderType','Line', 'HighlightColor',[0 0 0],...
            'Units','Pixels','Position',[001 195 400 055].*WS);
        hGraphUI   = uipanel('Parent',hGraphUIbg,'Title','','Tag','hGraphUI',...
            'BorderType','Etchedin', ...
            'BackgroundColor',[179 179 179]/255,...
            'Units','Pixels','Position',[005 005 390 046].*WS);
        hGraphUIButtsF = uipanel('Parent',hGraphUI,'Title','plot controls',...
            'Units','Pixels','Position',[005 005 380 039].*WS,...
            'BorderType','Beveledout', ...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
 
        %-Controls
        %------------------------------------------------------------------
        h1 = uicontrol('Parent',hGraphUIButtsF,'Style','CheckBox','String','hold',...
            'ToolTipString','toggle hold to overlay plots',...
            'FontSize',FS(10),...
            'Value',double(strcmp(get(hAx,'NextPlot'),'add')),...
            'Callback',[...
            'if get(gcbo,''Value''), ',...
            'set(get(gcbo,''UserData''),''NextPlot'',''add''), ',...
            'else, ',...
            'set(get(gcbo,''UserData''),''NextPlot'',''replace''), ',...
            'end'],...
            'Interruptible','on','Enable','on',...
            'Tag','holdButton',...
            'Position',[005 005 070 020].*WS);
        h2 = uicontrol('Parent',hGraphUIButtsF,'Style','CheckBox','String','grid',...
            'ToolTipString','toggle axes grid',...
            'FontSize',FS(10),...
            'Value',double(strcmp(get(hAx,'XGrid'),'on')),...
            'Callback',[...
            'if get(gcbo,''Value''), ',...
            'set(get(gcbo,''UserData''),''XGrid'',''on'','...
            '''YGrid'',''on'',''ZGrid'',''on''), ',...
            'else, ',...
            'set(get(gcbo,''UserData''),''XGrid'',''off'','...
            '''YGrid'',''off'',''ZGrid'',''off''), ',...
            'end'],...
            'Interruptible','on','Enable','on',...
            'Position',[080 005 070 020].*WS);
        h3 = uicontrol('Parent',hGraphUIButtsF,'Style','CheckBox','String','Box',...
            'ToolTipString','toggle axes box',...
            'FontSize',FS(10),...
            'Value',double(strcmp(get(hAx,'Box'),'on')),...
            'Callback',[...
            'if get(gcbo,''Value''), ',...
            'set(get(gcbo,''UserData''),''Box'',''on''), ',...
            'else, ',...
            'set(get(gcbo,''UserData''),''Box'',''off''), ',...
            'end'],...
            'Interruptible','on','Enable','on',...
            'Position',[155 005 070 020].*WS);
        h4 = uicontrol('Parent',hGraphUIButtsF,'Style','popupmenu',...
            'ToolTipString','edit axis text annotations',...
            'FontSize',FS(10),...
            'String','text|Title|Xlabel|Ylabel',...
            'Callback','spm_results_ui(''PlotUiCB'')',...
            'Interruptible','on','Enable','on',...
            'Position',[230 005 070 020].*WS);
        h5 = uicontrol('Parent',hGraphUIButtsF,'Style','popupmenu',...
            'ToolTipString','change various axes attributes',...
            'FontSize',FS(10),...
            'String','attrib|LineWidth|XLim|YLim|handle',...
            'Callback','spm_results_ui(''PlotUiCB'')',...
            'Interruptible','off','Enable','on',...
            'Position',[305 005 070 020].*WS);
 
        %-Handle storage for linking, and DeleteFcns for linked deletion
        %------------------------------------------------------------------
        set(hGraphUI,'UserData',[h1,h2,h3,h4,h5])
        set([h1,h2,h3,h4,h5],'UserData',hAx)
 
        set(hGraphUIbg,'UserData',...
            [hGraphUI,hGraphUIButtsF,h1,h2,h3,h4,h5],...
            'DeleteFcn','spm_results_ui(''Delete'',get(gcbo,''UserData''))')
        set(hAx,'UserData',hGraphUIbg,...
            'DeleteFcn','spm_results_ui(''Delete'',get(gcbo,''UserData''))')

 
 
 
    %======================================================================
    case 'plotuicb'
    %======================================================================
        % cg_tfce_results('PlotUiCB')
        hPM = gcbo;
        v   = get(hPM,'Value');
        if v==1, return, end
        str = cellstr(get(hPM,'String'));
        str = str{v};
 
        hAx = get(hPM,'UserData');
        switch str
            case 'Title'
                h = get(hAx,'Title');
                set(h,'String',spm_input('Enter title:',-1,'s+',get(h,'String')))
            case 'Xlabel'
                h = get(hAx,'Xlabel');
                set(h,'String',spm_input('Enter X axis label:',-1,'s+',get(h,'String')))
            case 'Ylabel'
                h = get(hAx,'Ylabel');
                set(h,'String',spm_input('Enter Y axis label:',-1,'s+',get(h,'String')))
            case 'LineWidth'
                lw = spm_input('Enter LineWidth',-1,'e',get(hAx,'LineWidth'),1);
                set(hAx,'LineWidth',lw)
            case 'XLim'
                XLim = spm_input('Enter XLim',-1,'e',get(hAx,'XLim'),[1,2]);
                set(hAx,'XLim',XLim)
            case 'YLim'
                YLim = spm_input('Enter YLim',-1,'e',get(hAx,'YLim'),[1,2]);
                set(hAx,'YLim',YLim)
            case 'handle'
                varargout={hAx};
            otherwise
                warning(['Unknown action: ',str])
        end
 
        set(hPM,'Value',1)
 
 
    %======================================================================
    case {'clear'}                                  %-Clear results subpane
    %======================================================================
        % Fgraph = cg_tfce_results('Clear',F,mode)
        % mode 1 [default] usual, mode 0 - clear & hide Res stuff, 2 - RNP
        
        if nargin<3, mode=1; else, mode=varargin{3}; end
        if nargin<2, F='Graphics'; else, F=varargin{2}; end
        F = spm_figure('FindWin',F);
 
        %-Clear input objects from 'Interactive' window
        %------------------------------------------------------------------
        %spm_input('!DeleteInputObj')
 
 
        %-Get handles of objects in Graphics window & note permanent results objects
        %------------------------------------------------------------------
        H = get(F,'Children');                          %-Get contents of window
        H = findobj(H,'flat','HandleVisibility','on');  %-Drop GUI components
        h = findobj(H,'flat','Tag','PermRes');          %-Look for 'PermRes' object
 
        if ~isempty(h)
            %-Found 'PermRes' object
            % This has handles of permanent results objects in it's UserData
            tmp  = get(h,'UserData');
            HR   = tmp.H;
            HRv  = tmp.Hv;
        else
            %-No trace of permanent results objects
            HR   = [];
            HRv  = {};
        end
        H = setdiff(H,HR);              %-Drop permanent results obj
 
 
        %-Delete stuff as appropriate
        %------------------------------------------------------------------
        if mode==2  %-Don't delete axes with NextPlot 'add'
            H = setdiff(H,findobj(H,'flat','Type','axes','NextPlot','add'));
        end
 
        delete(H)
 
        if mode==0  %-Hide the permanent results section stuff
            set(HR,'Visible','off')
        else
            set(HR,{'Visible'},HRv)
        end
 
 
    %======================================================================
    case 'close'                                            %-Close Results
    %======================================================================
        spm_clf('Interactive');
        spm_clf('Graphics');
        close(spm_figure('FindWin','Satellite'));
        evalin('base','clear');
    

    %======================================================================
    case 'launchmp'                            %-Launch multiplanar toolbox
    %======================================================================
        % hMP = cg_tfce_results('LaunchMP',M,DIM,hReg,hBmp)
        if nargin<5, hBmp = gcbo; else, hBmp = varargin{5}; end
        hReg = varargin{4};
        DIM  = varargin{3};
        M    = varargin{2};
 
        %-Check for existing MultiPlanar toolbox
        hMP  = get(hBmp,'UserData');
        if ishandle(hMP)
            figure(ancestor(hMP,'figure'));
            varargout = {hMP};
            return
        end
 
        %-Initialise and cross-register MultiPlanar toolbox
        hMP = spm_XYZreg_Ex2('Create',M,DIM);
        spm_XYZreg('Xreg',hReg,hMP,'spm_XYZreg_Ex2');
 
        %-Setup automatic deletion of MultiPlanar on deletion of results controls
        set(hBmp,'Enable','on','UserData',hMP)
        set(hBmp,'DeleteFcn','cg_tfce_results(''delete'',get(gcbo,''UserData''))')
 
        varargout = {hMP};
 
  
    %======================================================================
    case 'delete'                           %-Delete HandleGraphics objects
    %======================================================================
        % cg_tfce_results('Delete',h)
        h = varargin{2};
        delete(h(ishandle(h)));
 
 
    %======================================================================
    otherwise
    %======================================================================
        error('Unknown action string')
 
    %======================================================================
end

%==========================================================================
function mychgcon(obj,evt,xSPM)
%==========================================================================
xSPM2.swd   = xSPM.swd;
try, xSPM2.units = xSPM.units; end
xSPM2.Ic    = getfield(get(obj,'UserData'),'Ic');
if isempty(xSPM2.Ic) || all(xSPM2.Ic == 0), xSPM2 = rmfield(xSPM2,'Ic'); end
xSPM2.Im    = xSPM.Im;
xSPM2.pm    = xSPM.pm;
xSPM2.Ex    = xSPM.Ex;
xSPM2.title = '';
if ~isempty(xSPM.thresDesc)
    if strcmp(xSPM.STAT,'P')
        % These are soon overwritten by spm_getSPM
        xSPM2.thresDesc = xSPM.thresDesc;
        xSPM2.u = xSPM.u;
        xSPM2.k = xSPM.k;
        % xSPM.STATstr contains Gamma
    else
        td = regexp(xSPM.thresDesc,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
        if isempty(td)
            td = regexp(xSPM.thresDesc,'\w=(?<u>[\.\d]+)','names');
            td.thresDesc = 'none';
        end
        if strcmp(td.thresDesc,'unc.'), td.thresDesc = 'none'; end
        xSPM2.thresDesc = td.thresDesc;
        xSPM2.u     = str2double(td.u);
        xSPM2.k     = xSPM.k;
    end
end
hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
xyz  = spm_XYZreg('GetCoords',hReg);
[hReg,xSPM,SPM] = cg_tfce_results('setup',xSPM2);
TabDat = cg_tfce_list('List',xSPM,hReg);
spm_XYZreg('SetCoords',xyz,hReg);
assignin('base','hReg',hReg);
assignin('base','xSPM',xSPM);
assignin('base','SPM',SPM);
assignin('base','TabDat',TabDat);
figure(spm_figure('GetWin','Interactive'));

%==========================================================================
function mycheckres(obj,evt,xSPM)
%==========================================================================
spm_check_results([],xSPM);

%==========================================================================
function mysavespm(action)
%==========================================================================
xSPM = evalin('base','xSPM;');
XYZ  = xSPM.XYZ;

switch lower(action)
    case 'thresh'
        Z = xSPM.Z;
        
    case 'binary'
        Z = ones(size(xSPM.Z));
        
    case 'n-ary'
        if ~isfield(xSPM,'G')
            Z       = spm_clusters(XYZ);
            num     = max(Z);
            [n, ni] = sort(histc(Z,1:num), 2, 'descend');
            n       = size(ni);
            n(ni)   = 1:num;
            Z       = n(Z);
        else
            C       = NaN(1,size(xSPM.G.vertices,1));
            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
            C       = spm_mesh_clusters(xSPM.G,C);
            Z       = C(xSPM.XYZ(1,:));
        end
        
    case 'current'
        [xyzmm,i] = spm_XYZreg('NearestXYZ',...
            spm_results_ui('GetCoords'),xSPM.XYZmm);
        spm_results_ui('SetCoords',xSPM.XYZmm(:,i));
        
        if ~isfield(xSPM,'G')
            A   = spm_clusters(XYZ);
            j   = find(A == A(i));
            Z   = ones(1,numel(j));
            XYZ = xSPM.XYZ(:,j);
        else
            C   = NaN(1,size(xSPM.G.vertices,1));
            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
            C   = spm_mesh_clusters(xSPM.G,C);
            C   = C==C(xSPM.XYZ(1,i));
            Z   = C(xSPM.XYZ(1,:));
        end
        
    otherwise
        error('Unknown action.');
end

if isfield(xSPM,'G')
    F     = spm_input('Output filename',1,'s');
    if isempty(spm_file(F,'ext'))
        F = spm_file(F,'ext','.gii');
    end
    F     = spm_file(F,'CPath');
    M     = gifti(xSPM.G);
    C     = zeros(1,size(xSPM.G.vertices,1));
    C(xSPM.XYZ(1,:)) = Z; % or use NODE_INDEX
    M.cdata = C;
    save(M,F);
    cmd   = 'spm_mesh_render(''Disp'',''%s'')';
else
    V   = spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k));
    cmd = 'spm_image(''display'',''%s'')';
    F   = V.fname;
end

fprintf('Written %s\n',spm_file(F,'link',cmd));                         %-#

%==========================================================================
function myslover
%==========================================================================
spm_input('!DeleteInputObj');
xSPM = evalin('base','xSPM;');

so = slover;
[img,sts] = spm_select(1,'image','Select image for rendering on');
if ~sts, return; end
so.img.vol = spm_vol(img);
%obj.img.type = 'truecolour';
%obj.img.cmap = gray;
%[mx,mn] = slover('volmaxmin', obj.img.vol);
%obj.img.range = [mn mx];
so.img.prop = 1;

so = add_spm(so,xSPM);

so.transform = deblank(spm_input('Image orientation', '+1', ...
    'Axial|Coronal|Sagittal', char('axial','coronal','sagittal'), 1));
so = fill_defaults(so);
slices = so.slices;
so.slices = spm_input('Slices to display (mm)', '+1', 'e', ...
    sprintf('%0.0f:%0.0f:%0.0f',slices(1),mean(diff(slices)),slices(end)));

so.figure = spm_figure('GetWin', 'SliceOverlay');
so = paint(so);
assignin('base','so',so);

