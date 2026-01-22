% generate_PIPELINE_simulation_report.m
%
% Purpose:
%   Generate a PDF report using MATLAB Report Generator from precomputed
%   simulation assets. This script does not recompute dRSA; it loads
%   metadata and figures saved by simulations/build_simulation_report_assets.m.
%
% Example usage (from repo root in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   run('simulations/build_simulation_report_assets.m');
%   run('simulations/generate_PIPELINE_simulation_report.m');
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/build_simulation_report_assets.m');
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/generate_PIPELINE_simulation_report.m');
%
% Inputs:
%   - simulations/output/subXX/report_assets_subXX.mat
%   - simulations/output/subXX/report_figures/*.png
%
% Outputs:
%   - simulations/report/subXX/PIPELINE_simulation_report_RG_subXX.pdf
%
% Key assumptions:
%   - Report Generator is licensed (license('test','MATLAB_Report_Gen') == 1).
%   - Report assets were generated with build_simulation_report_assets.m.

%% Resolve paths and add dependencies
% Data flow: script location -> simulations dir -> repo root -> addpath.
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    % Fallback for unsaved or moved scripts: resolve by name on MATLAB path.
    scriptPath = which('generate_PIPELINE_simulation_report.m');
end
scriptDir = fileparts(scriptPath);
simDir = scriptDir;
repoRoot = fileparts(simDir);
addpath(simDir);
addpath(repoRoot);

%% Report configuration
% Data flow: load precomputed assets -> report sections.
participantNumber = 98;
subjectLabel = sprintf('sub%02d', participantNumber);
reportDir = fullfile(scriptDir, 'report', subjectLabel);
if ~exist(reportDir, 'dir')
    mkdir(reportDir);
end
reportOutput = fullfile(reportDir, sprintf('PIPELINE_simulation_report_RG_%s.pdf', subjectLabel));
assetsFile = fullfile(simDir, 'output', subjectLabel, ...
    sprintf('report_assets_%s.mat', subjectLabel));
if ~exist(assetsFile, 'file')
    error(['Report assets not found. Run simulations/build_simulation_report_assets.m ', ...
        'before generating the report. Missing file: %s'], assetsFile);
end
assets = load(assetsFile);
conditions = assets.conditions;
params = assets.params;
paramsPCR = assets.paramsPCR;
figureOutputDir = assets.figureOutputDir;
includePCR = true;
if isfield(assets, 'includePCR')
    includePCR = assets.includePCR;
end

%% Initialize the report
% Data flow: report metadata -> Report Generator document skeleton.
import mlreportgen.report.*
import mlreportgen.dom.*

rpt = Report(reportOutput, 'pdf');
rpt.Layout.Landscape = true;

titlePg = TitlePage();
titlePg.Title = 'PIPELINE Simulation Report';
titlePg.Subtitle = sprintf('Participant %02d', participantNumber);
titlePg.Author = 'Report Generator (MATLAB)';
add(rpt, titlePg);
add(rpt, TableOfContents);

%% Add parameter summary section
% Data flow: params struct -> report table.
paramSection = Section('Title', 'Key Parameters');
add(paramSection, Paragraph('Important parameters used in the simulations:'));

paramTableData = {
    'Parameter', 'Value';
    'modelDistMeasure', strjoin(string(params.modelDistMeasure), ', ');
    'neuralDistMeasure', string(params.neuralDistMeasure);
    'dRSAtype (corr)', string(params.dRSAtype);
    'dRSAtype (PCR)', string(paramsPCR.dRSAtype);
    'includePCR', string(includePCR);
};
paramTable = Table(paramTableData);
paramTable.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
paramTable.TableEntriesStyle = {VAlign('middle')};
add(paramSection, paramTable);
add(rpt, paramSection);

%% Add grouped dot-path plots
% Data flow: precomputed figures -> report section.
pathsSection = Section('Title', 'Dot paths (grouped)');
add(pathsSection, Paragraph('Each subplot shows dot trajectories per condition.'));
add_image(pathsSection, fullfile(figureOutputDir, ...
    sprintf('dot_paths_grouped_%s.png', subjectLabel)), ...
    sprintf('Dot paths grouped by condition and dot (%s)', subjectLabel));
add(pathsSection, Paragraph('Mean distance-from-center profiles averaged across conditions.'));
add_image(pathsSection, fullfile(figureOutputDir, ...
    sprintf('dot_paths_position_distribution_%s.png', subjectLabel)), ...
    sprintf('Distance from center (mean across conditions, %s)', subjectLabel));
add(rpt, pathsSection);

%% Add grouped position-time distance matrices
% Data flow: precomputed figures -> report section.
distSection = Section('Title', 'Position-time distance matrices (grouped)');
add(distSection, Paragraph('Each subplot shows the time-by-time distance matrix per dot and condition.'));
add_image(distSection, fullfile(figureOutputDir, ...
    sprintf('distance_matrices_grouped_%s.png', subjectLabel)), ...
    sprintf('Position-time distance matrices grouped by condition and dot (%s)', subjectLabel));
add(rpt, distSection);

%% Add combined dRSA matrices page
% Data flow: precomputed figures -> report section.
drsaSection = Section('Title', 'dRSA matrices');
add(drsaSection, Paragraph('Top row: position data (dot 1). Bottom row: direction data (dot 1).'));
add_image(drsaSection, fullfile(figureOutputDir, ...
    sprintf('drsa_matrices_corr_%s.png', subjectLabel)), ...
    sprintf('dRSA matrices (corr, %s)', subjectLabel));
if includePCR
    add(drsaSection, Paragraph('PCR matrices shown below.'));
    add_image(drsaSection, fullfile(figureOutputDir, ...
        sprintf('drsa_matrices_pcr_%s.png', subjectLabel)), ...
        sprintf('dRSA matrices (PCR, %s)', subjectLabel));
end
add(rpt, drsaSection);

%% Add combined lagged plots page
% Data flow: precomputed figures -> report section.
lagSection = Section('Title', 'Lagged dRSA');
add(lagSection, Paragraph('Two panels: position data and direction data from dot 1.'));
add_image(lagSection, fullfile(figureOutputDir, ...
    sprintf('lagged_drsa_%s.png', subjectLabel)), ...
    sprintf('Lagged dRSA plots combined (%s)', subjectLabel));
add(rpt, lagSection);

%% Finalize the report
% Data flow: report content -> PDF on disk.
close(rpt);
rptview(rpt);

%% Local helper functions
% Data flow: add precomputed images into the report.

function add_image(sectionHandle, imgFile, captionText)
% add_image
%
% Purpose:
%   Add a precomputed image to a report section with a caption.
%
% Inputs:
%   sectionHandle : mlreportgen.report.Section
%   imgFile       : full path to the image file
%   captionText   : caption for the image

    import mlreportgen.dom.*
    if ~exist(imgFile, 'file')
        error('Report image not found: %s', imgFile);
    end
    img = Image(imgFile);
    % Scale image to fit the page width and preserve aspect ratio.
    img.Style = {ScaleToFit(true)};
    add(sectionHandle, img);
    add(sectionHandle, Paragraph(captionText));
end
