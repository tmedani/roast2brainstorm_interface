function Rost2Brainstorm(varargin)
% Rost2Brainstorm : Create a Brainstorm protocol and load the Roast outputs.
% Usages :
% Rost2Brainstorm()
% Rost2Brainstorm(RoastOutputDir, fmriFileName)
% Rost2Brainstorm(RoastOutputDir, fmriFileName, BrainstormDbDir)
% Rost2Brainstorm(RoastOutputDir, fmriFileName, BrainstormDbDir,options)
% Rost2Brainstorm(options)
% Inputs arguments
% RoastOutputDir : Complete path to the Roast outputs, 
%                          default :  roast/example
% fmriFileName : Name of the original fMRI file, 
%                          default :  'MNI152_T1_1mm.nii'
% BrainstormDbDir : Complete path to the Brainstom data base, 
%                           default : default Brainstorm folder 
%                           if the specified folder does not esxit, i will be created           
% options : Matlab structures with the following fields 
% BrainstormDbDir: Similar as above
% RoastOutputDir: Similar as above
% RoastHomeDir: path where roast is located
% BrainstormHomeDir: path where brainstorm is located
% ProtocolName: brainstorm protocol's name, default : 'roast2brainstorm'
% SubjectName: brainstoem subject name, default : 'Subject01'
% inputSubjectMri: name if the input MRI to roast, default: 'MNI152_T1_1mm.nii'
% viewMRI: view the MRI with Brainstorm MRI viewer, default :0
% editMriFiducial: Edit the fiducial point from the MRI viewer, defaul : 0
% computeMNI: compute the MNI transformation, require SPM, default :1
% addFullHeadLf: load the full LF matrix to Brainstorm, default :0
% computeCatCortex: extract the corex with spm-cat, reauires cat, default : 1
% convert2scsCoordinate: not used, default 1
% replaceAirByScalp: replace the air cavities by scalp tissues, default 0
% deletePreviousProtocol: delete the brainstorm protocole, default 1
% extractFemCortex: extract the FEM cortex from the FEM mesh, default 1
% removeElectrodMesh: remove the mesh of the electrodes from the mesh, default :1
% computeSpmSurfaces: compute the spm canonical surfaces, require spm, default 1

% Requirement : 
% Have already run the cmd : roast([],'leadField','simulationTag','MNI152leadField')  
%                                           Need to have brainstorm within the Matlab path
%                                           Some optiion reauires SPM and
%                                           CAT tooolbox. 

% Authors
% Takfarinas MEDANI,  Maximilian Nentwich, Stefan Haufe, Yu Andy Huang
% Discussion : Link to roast github discussion https://github.com/andypotatohy/roast/issues/25

% Read the Roast output data and convet to brainstorm protocol
% Create the roast2brainstorm_protocol with a subject. The anatomy folder contains the MRI,
% the FEM mesh of the head model and the outersurface of the cortex extracted from the FEM mesh
% The functional folder conains the location of the electrods and the LF (forward model)
% for the cortex and for the whole head (if specified in the opts).
% check the available options within the function roast2brainstorm_defaults 


%% ===== START =====
%% Check the if brainstorm is available

defOPTIONS = roast2brainstorm_defaults();
OPTIONS = [];
% Open the bst waitBar
bst_progress('start','roast2brainstorm','Check the inputs')
% Check inputs
if (nargin < 1) || isempty(varargin{1})
    OPTIONS.RoastOutputDir = defOPTIONS.RoastOutputDir;
    OPTIONS.BrainstormDbDir = defOPTIONS.BrainstormDbDir;
end
if (nargin == 2) && ischar(varargin{1}) && ischar(varargin{2})
    OPTIONS.RoastOutputDir = varargin{1};
    OPTIONS.inputSubjectMri = varargin{2};
    OPTIONS.BrainstormDbDir = defOPTIONS.BrainstormDbDir;
end

if (nargin == 2) && isstruct(varargin{1})
    OPTIONS = struct_copy_fields(opts, defOPTIONS, 0);
end

if (nargin == 3) && ischar(varargin{1}) && ischar(varargin{3})
    OPTIONS.RoastOutputDir = varargin{1};
    OPTIONS.BrainstormDbDir = varargin{3};
end

if (nargin == 4) && ischar(varargin{1}) && ischar(varargin{3}) && isstruct(varargin{3})
    OPTIONS.RoastOutputDir = varargin{1};
    OPTIONS.BrainstormDbDir = varargin{2};
    OPTIONS = struct_copy_fields(opts, defOPTIONS, 0);
end

if ~isfolder(OPTIONS.BrainstormDbDir)
    if ~mkdir(OPTIONS.BrainstormDbDir)
        bst_progress('stop');
        error(['Could not create folder: ' OPTIONS.BrainstormDbDir]);
    end
    OPTIONS.BrainstormDbDir =   fullfile(OPTIONS.BrainstormDbDir, 'brainstorm_db') ;
    mkdir(BrainstormDbDir);
end
% Update
OPTIONS = struct_copy_fields(OPTIONS, defOPTIONS, 0);

% Initialize random generator
if exist('rng', 'file')
    rng('default');
else
    sRnd = RandStream.getDefaultStream;
    sRnd.reset();
end

% The used folders
disp([ 'RoastOutputDir : ' OPTIONS.RoastOutputDir])
disp([ 'BrainstormDbDir : ' OPTIONS.BrainstormDbDir])
optionTable = struct2table(OPTIONS,'AsArray',true)

% call the roast LF default computation <== TODO for other subjects
if OPTIONS.runRoastLeadField
    currentDir = pwd;
    cd(RoastHomeDir) 
    roast([],'leadField','simulationTag','MNI152leadField');
    cd(currentDir)
end

%% start brainstorm
% Start profiling
% profile on
% Start Brainstorm with GUI
bst_startup(OPTIONS.BrainstormHomeDir, 1, OPTIONS.BrainstormDbDir);
% The protocol name has to be a valid folder name (no spaces, no weird characters...)
bst_progress('text', 'Rosat2Brainstorm: Create Brainstorm protocol...');
% Delete existing protocol
if OPTIONS.deletePreviousProtocol == 1
    gui_brainstorm('DeleteProtocol', OPTIONS.ProtocolName);
end
% Create new protocol
gui_brainstorm('CreateProtocol', OPTIONS.ProtocolName, 0, 0, OPTIONS.BrainstormDbDir);
% Start a new report
% bst_report('Start');
%% ===== CREATE SUBJECT =====
% Create subject
[sSubject, iSubject] = db_add_subject(OPTIONS.SubjectName);
% Input files
sFiles = [];
SubjectNames = {OPTIONS.SubjectName};
RawFiles = {fullfile(OPTIONS.RoastOutputDir, OPTIONS.inputSubjectMri)};
%%% Create a study
StudyName = 'roast';
iStudy = db_add_condition(OPTIONS.SubjectName, StudyName);

% Process: Import MRI
bst_progress('text', 'Rosat2Brainstorm: Import MRI...');
sFiles = bst_process('CallProcess', 'process_import_mri', sFiles, [], ...
    'subjectname', SubjectNames{1}, ...
    'mrifile',     {RawFiles{1}, 'Nifti1'}, ...
    'nas',         [0, 0, 0], ...
    'lpa',         [0, 0, 0], ...
    'rpa',         [0, 0, 0], ...
    'ac',          [0, 0, 0], ...
    'pc',          [0, 0, 0], ...
    'ih',          [0, 0, 0]);

% Process: Compute MNI transformation
if OPTIONS.computeMNI  == 1
    sFiles = bst_process('CallProcess', 'process_mni_affine', sFiles, [], ...
        'subjectname', SubjectNames{1});
end
% % Edit the MRI and set the fiducial points
% % Get subject
[sSubject, iSubject] = bst_get('Subject', OPTIONS.SubjectName);
filenameRelative = sSubject.Anatomy(sSubject.iAnatomy).FileName;
if OPTIONS.viewMRI == 1
    view_mri(filenameRelative);
end
if OPTIONS.editMriFiducial ==1
    view_mri(filenameRelative, 'EditFiducials');
end

% if OPTIONS.computeCatCortex == 1
%     bst_call(@process_segment_cat12, 'ComputeInteractive', iSubject,sSubject.iAnatomy);
% %     nVertices = 15000;
% end
% if OPTIONS.computeSpmSurfaces == 1
%     bst_call(@process_generate_canonical, 'ComputeInteractive', iSubject, sSubject.iAnatomy)
% end
%% Load the FEM mesh
bst_progress('text', 'Rosat2Brainstorm: Loading Roast FEM Mesh...');
roast_headmodel = load(fullfile(OPTIONS.RoastOutputDir, ...
    strrep(OPTIONS.inputSubjectMri, '.nii', '_MNI152leadField.mat')));
sMri = load(file_fullpath(filenameRelative));
% Write to the Brainstorm format
FemMat = db_template('femmat');
% New surface structure
if OPTIONS.removeElectrodMesh == 1
    [no,el] = removeisolatednode(roast_headmodel.node(:,1:3),roast_headmodel.elem(roast_headmodel.elem(:,5)<=6,:));
    newelem = meshreorient(no,el(:,1:4));
    newelem = [newelem el(:,5)];
%    newelem(newelem(:,end) >= 7,end) = [];
    FemMat.TissueLabels = [{'1-wm'} {'2-gm'} {'3-csf'} {'4-skull'} {'5-scalp'} {'6-air'}];
    if OPTIONS.replaceAirByScalp == 1
        newelem(newelem(:,end) == 6,end) = 5;
        FemMat.TissueLabels = [{'1-wm'} {'2-gm'} {'3-csf'} {'4-skull'} {'5-scalp'}];
    end
else
    no = roast_headmodel.node(:,1:3);
    newelem = meshreorient(roast_headmodel.node(:,1:3),roast_headmodel.elem(:,1:4));
    newelem = [newelem roast_headmodel.elem(:,5)];
    FemMat.TissueLabels = [{'1-wm'} {'2-gm'} {'3-csf'} {'4-skull'} {'5-scalp'} {'6-air'} {'7-electrodes'}] ;
    if OPTIONS.replaceAirByScalp == 1
        newelem(newelem(:,end) == 6,end) = 5;
        FemMat.TissueLabels = [{'1-wm'} {'2-gm'} {'3-csf'} {'4-skull'} {'5-scalp'} {'7-electrodes'}];
    end
end
FemMat.Tissue = newelem(:,end);
% convert the coordinates
if OPTIONS.computeMNI == 1 % will convert to the scs coordinates
    newNode = cs_convert(sMri, 'voxel', 'scs', no(:,1:3));
    FemMat.Comment = ['FEM ' num2str(length(newNode)) '  (roast scs cs,  ' ...
        num2str(length( FemMat.Tissue))  ' layers) '];
else
    newNode =  no(:,1:3);
    FemMat.Comment = ['FEM ' num2str(length(newNode)) '  (roast original cs,  ' ...
        num2str(length( FemMat.Tissue))  ' layers) '];
end
FemMat.Vertices = newNode;
FemMat.Elements = newelem(:,1:4);
FemMat.Tensors = [];
% History: File name
FemMat = bst_history('add', FemMat,  ['FEM mesh imported from roast output']);
% Produce a default surface filename &   Make this filename unique
[filepath,tmp, tmp] = bst_fileparts(file_fullpath(filenameRelative));
FemMatFilename = file_unique(bst_fullfile(filepath, ...
    sprintf(['tess_fem_roast_%dV.mat'], length(FemMat.Vertices))));
% Save new surface in Brainstorm format
bst_progress('text', 'Rosat2Brainstorm: Save Roast FEM Mesh to Brainstorm...');
bst_save(FemMatFilename, FemMat, 'v7');
db_add_surface(iSubject, FemMatFilename, FemMat.Comment);

%% Add the cortex from roast
if OPTIONS.extractFemCortex == 1
    bst_progress('text', 'Rosat2Brainstorm: Extract the Roast FEM cortex...');
    wmGmIndex = 2;
    cortexMeshElement = roast_headmodel.elem(find(roast_headmodel.elem(:,5)<=wmGmIndex),:);
    extractCortexSurface = volface(cortexMeshElement(:,1:4));
    surf.vertices = roast_headmodel.node(:,1:3);
    surf.faces = extractCortexSurface;
    [surf,UsedV] = delete_unused_vertices(surf);
    surf.faces = meshreorient(surf.vertices, surf.faces);
    if OPTIONS.computeMNI == 1 % will convert to the scs coordinates
        % Tramsform the coordinates
        surf.vertices = bsxfun(@plus, surf.vertices * sMri.SCS.R', sMri.SCS.T');
        surf.faces = meshreorient(surf.vertices, surf.faces);
        % Convert the units
        surf.vertices = surf.vertices/1e3;
    end
    bst_progress('text', 'Rosat2Brainstorm: Save FEM cortex to Brainstorm...');
    % Save to braisntoem database
    CortexMat = db_template('surfacemat');
    CortexMat.Comment = 'cortex (outer fem)';
    CortexMat.Vertices = surf.vertices;
    CortexMat.Faces = surf.faces;
    % History: File name
    CortexMat = bst_history('add', CortexMat, 'create', 'Cortex extracted from Roast FEM model');
    % Produce a default surface filename &   Make this filename unique
    CortexFile = file_unique(bst_fullfile(filepath, ...
        sprintf('tess_cortex_roast_%dV.mat', length(CortexMat.Vertices))));
    bst_save(CortexFile, CortexMat, 'v7');
    db_add_surface(iSubject, CortexFile, 'cortex from roast');
    % Update subject node
    panel_protocols('UpdateNode', 'Subject', iSubject);
end

% The cortical surface is already present in the roast output and can be converted to brainstorm format
% The index 2 labels grey matter; extracting all faces with label two gives the cortical surface
% This also includes the surface between grey matter and white matter

%% Add the channel
% TODO : get these coordinates from the Roast
bst_progress('text', 'Rosat2Brainstorm: Load the channels file...');
electrode_coord = load('lf_electrode_coord.mat');
load('MNI152_T1_1mm_indInRoastCore.mat', 'indInRoastCore');
% Find the real electrode names
fid = fopen(fullfile(OPTIONS.RoastHomeDir,'./elec72.loc')); C = textscan(fid,'%d %f %f %s'); fclose(fid);
elecName = C{4}; for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end
if OPTIONS.computeMNI == 1 % will convert to the scs coordinates
    % Tramsform the coordinates
    electrode_coord = bsxfun(@plus, electrode_coord.electrode_coord * sMri.SCS.R', sMri.SCS.T');
    % Convert the units
    electrode_coord = electrode_coord/1e3;
end
% Load the indices to reorder the electrodes correctly
indexThatReorderLF = indInRoastCore;
indexThatReorderLF(end) = 72; % 72 is the ref elec
electrode_coordReorder = electrode_coord(indexThatReorderLF,:);
for iCh = 1:length(electrode_coordReorder)
    ChannelMat.Channel(iCh).Type = 'EEG';
    ChannelMat.Channel(iCh).Name = elecName{iCh};
    ChannelMat.Channel(iCh).Loc = electrode_coordReorder(iCh,:)';
    ChannelMat.Channel(iCh).Orient = [];
    ChannelMat.Channel(iCh).Comment = [];
    ChannelMat.Channel(iCh).Weight = 1;
    ChannelMat.Channel(iCh).Group = [];
end
ChannelMat.Comment =  ['roast channel (' num2str(length(electrode_coord)) ')'];
ChannelMat.HeadPoints.Loc = [];
ChannelMat.HeadPoints.Label = [];
ChannelMat.HeadPoints.Type = [];
ChannelMat = bst_history('add', ChannelMat, 'import', ['Import from: rosat']);
ChannelMat.IntraElectrodes = [];
ChannelMat.MegRefCoef = [];
ChannelMat.Projector = struct('Comment', [], 'Components', [], 'CompMask', [], 'Status', [], 'StringVal', []);
ChannelMat.SCS = struct('NAS', [], 'LPA', [], 'RPA', [], 'R', [], 'T', []);
ChannelMat.TransfEeg = [];
ChannelMat.TransfEegLabels = {'manual correction'};
ChannelMat.TransfMeg = [];
ChannelMat.TransfMegLabels = [];

bst_progress('text', 'Rosat2Brainstorm: Save to Brainstorm...');
ChannelMatFile = file_unique(bst_fullfile([strrep(filepath, 'anat', 'data'), '/', StudyName], 'channel_roast.mat'));
bst_save(ChannelMatFile, ChannelMat, 'v7');
db_add(iStudy, ChannelMatFile, 'electrodes from roast');

%% Import the Leadfield and save to data base
% FEM cortex
% Load the leadfield
bst_progress('text', 'Rosat2Brainstorm: Load the Roast Leadfield...');
roast_leadfield = load([OPTIONS.RoastOutputDir, '/', ...
    strrep(OPTIONS.inputSubjectMri, '.nii', '_MNI152leadField_roastResult.mat')]);
% Add zeros for the reference
roast_leadfield.A_all = cat(3, ...
    roast_leadfield.A_all, zeros(size(roast_leadfield.A_all,1), size(roast_leadfield.A_all,2)));
% Permute the array from vertex x dimension x channel to channel x vertex x dimension
% (dimension is x,y,z in which coordinate system?)
roast_leadfield.A_all2 = permute(roast_leadfield.A_all,[3 1 2]);

if OPTIONS.extractFemCortex == 1
    % Select only the vertices on the cortex found in the previous steps
    roast_leadfield.A_cortex = roast_leadfield.A_all2(:,UsedV,:);
    if OPTIONS.computeMNI == 1 % will convert to the scs coordinates
        % Transform the leadfield to match the coordinates of the MRI
        for ch = 1:size(roast_leadfield.A_cortex,1)
            roast_leadfield.A_trans(ch,:,:) = squeeze(roast_leadfield.A_cortex(ch,:,:)) * sMri.SCS.R';
            %         roast_leadfield.Ahead_trans(ch,:,:) = squeeze(roast_leadfield.All_head(ch,:,:)) * sMri.SCS.R';
        end
    end
    lfAll_cortex = zeros(size(roast_leadfield.A_trans, 1), size(roast_leadfield.A_trans, 2)*3);
    lfAll_cortex(:, 1:3:end) = roast_leadfield.A_trans(:,:,1);    % x
    lfAll_cortex(:, 2:3:end) = roast_leadfield.A_trans(:,:,2);    % y
    lfAll_cortex(:, 3:3:end) = roast_leadfield.A_trans(:,:,3);    % z
    
    %% Save another file that saves all variables individually
    bst_cortexLF.Comment = 'Headmodel from Roast';
    bst_cortexLF.ECOGMethod = '';
    bst_cortexLF.EEGMethod = 'Roast FEM';
    bst_cortexLF.Gain = lfAll_cortex;
    bst_cortexLF.GridAtlas = [];
    bst_cortexLF.GridLoc = surf.vertices;
    bst_cortexLF.GridOptions = [];
    bst_cortexLF.GridOrient = [];
    bst_cortexLF.HeadModelType = 'surface';
    bst_cortexLF.History = {'07-Sep-2020 13:08:04', 'compute', 'Roast FEM'};
    bst_cortexLF.MEGMethod = '';
    bst_cortexLF.Param = [];
    bst_cortexLF.SEEGMethod = '';
    bst_cortexLF.SurfaceFile = CortexFile;
    % Define the filename and path
    headmodel_file = file_unique(bst_fullfile([strrep(filepath, 'anat', 'data'), '/', StudyName], 'headmodel_surf_roast.mat'));
    bst_save(headmodel_file, bst_cortexLF, 'v7');
    db_add_data(iStudy, headmodel_file, bst_cortexLF);
    panel_protocols('UpdateTree')
end
%
%% add the whole head
if OPTIONS.addFullHeadLf == 1
    UsedVhead = roast_headmodel.elem(roast_headmodel.elem(:,5)<=6,:);
    UsedVhead=unique(UsedVhead(:));
    roast_leadfield.All_head = roast_leadfield.A_all2(:,UsedVhead,:);
    if OPTIONS.computeMNI == 1 % will convert to the scs coordinates
        % Transform the leadfield to match the coordinates of the MRI
        for ch = 1:size(roast_leadfield.A_cortex,1)
            roast_leadfield.A_trans(ch,:,:) = squeeze(roast_leadfield.A_cortex(ch,:,:)) * sMri.SCS.R';
        end
    end
    
    lfAll_head = zeros(size(roast_leadfield.Ahead_trans, 1), size(roast_leadfield.Ahead_trans, 2)*3);
    lfAll_head(:, 1:3:end) = roast_leadfield.Ahead_trans(:,:,1);    % x
    lfAll_head(:, 2:3:end) = roast_leadfield.Ahead_trans(:,:,2);    % y
    lfAll_head(:, 3:3:end) = roast_leadfield.Ahead_trans(:,:,3);    % z
    bst_HeasLF =[];
    bst_HeasLF.Comment = 'Headmodel from Roast all head';
    bst_HeasLF.ECOGMethod = '';
    bst_HeasLF.EEGMethod = 'Roast FEM';
    bst_HeasLF.Gain = lfAll_head;
    bst_HeasLF.GridAtlas = [];
    bst_HeasLF.GridLoc = FemMat.Vertices;
    bst_HeasLF.GridOptions = [];
    % bst_HeasLF.GridOrient = [];
    bst_HeasLF.HeadModelType = 'surface';
    bst_HeasLF.History = {'07-Sep-2020 13:08:04', 'compute', 'Roast FEM'};
    bst_HeasLF.MEGMethod = '';
    bst_HeasLF.Param = [];
    bst_HeasLF.SEEGMethod = '';
    bst_HeasLF.SurfaceFile = FemMatFilename;
    % Define the filename and path
    headmodel_file = file_unique(bst_fullfile([strrep(filepath, 'anat', 'data'), '/', StudyName], 'headmodel_surf_roastFullHead.mat'));
    bst_save(headmodel_file, bst_HeasLF, 'v7');
    db_add_data(iStudy, headmodel_file, bst_HeasLF);
    panel_protocols('UpdateTree')
end


% % Save and display report
% ReportFile = bst_report('Save', sFiles);
% bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);
end

function defOPTIONS = roast2brainstorm_defaults()
    % find roast dir
    mfilename = ('brainstorm.m');
    BrainstormHomeDir = fileparts(which(mfilename));
    if isempty(BrainstormHomeDir)
        error('Brainstorm should be available in your computer ( https://neuroimage.usc.edu/bst/download.php)');
    else
        disp('roast2brainstorm >> Great, you have brainstorm in your computer, let''s continue...')
    end
    % find roast dir
    mfilename = ('roast.m');
    RoastHomeDir = fileparts(which(mfilename));    
    RoastOutputDir = fullfile(RoastHomeDir,'example');
    % Should be defined by the user
%     RoastOutputDir = fullfile('/home/max/Desktop/MyMRI/Registered/');
    
    % start braisntom 
    % Start Brainstorm with GUI
    % bst_startup(BrainstormHomeDir, 0, BrainstormDbDir); < == TODO 
    brainstorm     
    
    % find Brainstorm Db dir
    BrainstormDbDir = bst_get('BrainstormDbDir');

    % outputs
    defOPTIONS.RoastHomeDir = RoastHomeDir;
    defOPTIONS.BrainstormDbDir = BrainstormDbDir;
    defOPTIONS.RoastOutputDir = RoastOutputDir;
    defOPTIONS.BrainstormHomeDir = BrainstormHomeDir;
    defOPTIONS.runRoastLeadField = 0; 
    defOPTIONS.ProtocolName = 'roast2brainstorm';
    defOPTIONS.SubjectName = 'Subject01';
    defOPTIONS.inputSubjectMri = 'MNI152_T1_1mm.nii';
    defOPTIONS.viewMRI = 0;
    defOPTIONS.editMriFiducial =0;
    defOPTIONS.computeMNI = 1;
    defOPTIONS.addFullHeadLf = 0;
    defOPTIONS.computeCatCortex = 0;
    defOPTIONS.convert2scsCoordinate = 1;
    defOPTIONS.replaceAirByScalp  = 0;
    defOPTIONS.deletePreviousProtocol = 1 ;
    defOPTIONS.extractFemCortex = 1;
    defOPTIONS.removeElectrodMesh = 1;
    defOPTIONS.computeSpmSurfaces =1;
    defOPTIONS.computeCatCortex =1;
  
    end

function [surf,UsedV]=delete_unused_vertices(surf)
%Author Anand A. Joshi ajoshi@sipi.usc.edu
%UsedV=reshape(surf.faces,size(surf.faces,1)*3,1);
UsedV=surf.faces(:);
UsedV=unique(UsedV);

%usedmarker=zeros(size(surf.vertices,1),1); usedmarker(UsedV)=1;
new_indx=zeros(size(surf.vertices,1),1);
%num_used=sum(usedmarker);
num_used=size(UsedV,1);
new_indx(UsedV)=[1:num_used];

surf.vertices=surf.vertices(UsedV,:);
surf.faces=new_indx(surf.faces);
% surf.tissue=new_indx(surf.tissue);
end
