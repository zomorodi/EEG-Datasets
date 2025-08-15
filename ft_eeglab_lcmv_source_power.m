% LCMV source analysis for 5 EEGLAB datasets using FieldTrip
% - Loads .set files
% - Prepares template headmodel and sourcemodel in MNI
% - Builds common LCMV filter on broadband (1-50 Hz)
% - Applies filter to band-limited data to estimate source power in 5 bands
% - Saves orthogonal slice figures per band per subject

% ------------- USER CONFIGURABLE -------------
% Update these paths before running
%eeglabDir = '/absolute/path/to/eeglab';
%fieldtripDir = '/absolute/path/to/fieldtrip';

% Absolute paths to exactly 5 EEGLAB .set files
filelist = {
	'/absolute/path/to/subj01.set'
	'/absolute/path/to/subj02.set'
	'/absolute/path/to/subj03.set'
	'/absolute/path/to/subj04.set'
	'/absolute/path/to/subj05.set'
};

% Where to write figures
outputDir = fullfile(pwd, 'source_figs');

% Analysis parameters
segmentLengthSec = 2;   % segment continuous data into 2 s trials
resampleHz = 250;       % resample for speed (set [] to skip)
freqBands = [1 4; 4 8; 8 12; 12 30; 30 50];
bandNames = {'delta','theta','alpha','beta','gamma'};
% --------------------------------------------

% Optional: add EEGLAB/FieldTrip to path if not already
if exist('eeglab', 'file') ~= 2 && exist('eeglabDir', 'var') && ~isempty(eeglabDir)
	addpath(genpath(eeglabDir));
end
if exist('ft_defaults', 'file') ~= 2 && exist('fieldtripDir', 'var') && ~isempty(fieldtripDir)
	addpath(fieldtripDir);
end

% Initialize FieldTrip
if exist('ft_defaults', 'file') ~= 2
	error('FieldTrip not on path. Define fieldtripDir above or add it to MATLAB path.');
end
ft_defaults;

% Make figure output directory
if ~exist(outputDir, 'dir')
	mkdir(outputDir);
end

% Render figures headless by default
set(0, 'DefaultFigureVisible', 'off');

% Resolve FieldTrip template paths
ftPath = fileparts(which('ft_defaults'));
templateDir = fullfile(ftPath, 'template');

% Load template headmodel (BEM) and sourcemodel (8 mm grid) in MNI coords
H = load(fullfile(templateDir, 'headmodel', 'standard_bem.mat'));
if isfield(H, 'vol')
	headmodel = H.vol;
elseif isfield(H, 'headmodel')
	headmodel = H.headmodel;
else
	error('Unexpected contents in standard_bem.mat');
end
S = load(fullfile(templateDir, 'sourcemodel', 'standard_sourcemodel3d8mm.mat'));
if isfield(S, 'sourcemodel')
	sourcemodel = S.sourcemodel;
elseif isfield(S, 'grid')
	sourcemodel = S.grid;
else
	error('Unexpected contents in standard_sourcemodel3d8mm.mat');
end
headmodel = ft_convert_units(headmodel, 'mm');
sourcemodel = ft_convert_units(sourcemodel, 'mm');

% Iterate over datasets
for si = 1:numel(filelist)
	setfile = filelist{si};
	if ~isfile(setfile)
		error('File does not exist: %s', setfile);
	end
	[setPath, setBase, setExt] = fileparts(setfile);
	if isempty(setPath)
		setPath = pwd;
	end

	% Load with EEGLAB
	if exist('pop_loadset', 'file') ~= 2
		error('EEGLAB not on path. Define eeglabDir above or add it to MATLAB path.');
	end
	EEG = pop_loadset('filename', [setBase setExt], 'filepath', setPath);
	EEG = eeg_checkset(EEG);

	% Convert to FieldTrip
	if exist('eeglab2fieldtrip', 'file') ~= 2
		error('eeglab2fieldtrip not found. Ensure FieldTrip''s EEGLAB interface is available.');
	end
	data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

	% Preprocess: average reference, 1-50 Hz, demean/detrend
	cfg = [];
	cfg.reref = 'yes';
	cfg.refchannel = 'all';
	cfg.demean = 'yes';
	cfg.detrend = 'yes';
	cfg.bpfilter = 'yes';
	cfg.bpfreq = [1 50];
	data = ft_preprocessing(cfg, data);

	% Resample (optional)
	if ~isempty(resampleHz) && isnumeric(resampleHz) && resampleHz > 0
		cfg = [];
		cfg.resamplefs = resampleHz;
		cfg.detrend = 'no';
		data = ft_resampledata(cfg, data);
	end

	% Segment continuous data into fixed-length trials
	cfg = [];
	cfg.length = segmentLengthSec;
	cfg.overlap = 0;
	data = ft_redefinetrial(cfg, data);

	% Build electrode structure from EEGLAB if 3D coords are present; else use template 10-20
	elec = [];
	if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'X') && ~isempty([EEG.chanlocs.X])
		labels = string({EEG.chanlocs.labels})';
		positions = [[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]'];
		elec.label = cellstr(labels);
		elec.chanpos = positions;
		elec.elecpos = positions;
		elec.unit = 'mm';
	else
		elec = ft_read_sens(fullfile(templateDir, 'electrode', 'standard_1020.elc'));
	end
	elec = ft_convert_units(elec, 'mm');

	% Align labels between data and electrodes
	common = intersect(data.label, elec.label, 'stable');
	if numel(common) < 16
		error('Too few overlapping channels between data and electrode montage for %s.', [setBase setExt]);
	end
	cfg = [];
	cfg.channel = common;
	data = ft_selectdata(cfg, data);

	sel = match_str(elec.label, common);
	elec.label = elec.label(sel);
	elec.chanpos = elec.chanpos(sel, :);
	elec.elecpos = elec.elecpos(sel, :);
	data.elec = elec;

	% Compute broadband covariance
	cfg = [];
	cfg.covariance = 'yes';
	cfg.covariancewindow = 'all';
	timelockBroad = ft_timelockanalysis(cfg, data);

	% Prepare leadfield (user-specified settings)
	cfg = [];
	cfg.normalize = 'yes';
	cfg.channel = 'EEG';
	cfg.coordsys = 'mni';
	cfg.elec = data.elec;
	cfg.headmodel = headmodel;
	cfg.sourcemodel = sourcemodel;
	lf = ft_prepare_leadfield(cfg, timelockBroad);

	% Ensure the leadfield carries an explicit, canonical channel list
	if ~isfield(lf, 'label') || isempty(lf.label)
		lf.label = timelockBroad.label;
	end
	if isstring(lf.label)
		lf.label = cellstr(lf.label);
	end
	lf.label = lf.label(:);

	% LCMV common spatial filter on broadband
	cfg = [];
	cfg.method = 'lcmv';
	cfg.keeptrials = 'yes';
	cfg.lcmv.keepfilter = 'yes';
	cfg.lcmv.lambda = '5%';
	cfg.lcmv.fixedori = 'yes';
	cfg.lcmv.projectnoise = 'yes';
	cfg.lcmv.weightnorm = 'arraygain'; % Alternative: 'nai'
	cfg.sourcemodel = lf;
	cfg.headmodel = headmodel;
	sourceBroad = ft_sourceanalysis(cfg, timelockBroad);

	% Build a sourcemodel with common filters and correct channel labels
	lfFilter = lf;
	if isfield(sourceBroad, 'avg') && isfield(sourceBroad.avg, 'filter')
		lfFilter.filter = sourceBroad.avg.filter;
		lfFilter.avg = struct();
		lfFilter.avg.filter = sourceBroad.avg.filter; % also store where FieldTrip expects it
	else
		error('Common LCMV filters not found in sourceBroad.avg.filter');
	end
	if isfield(sourceBroad, 'inside'); lfFilter.inside = sourceBroad.inside; end
	if isfield(sourceBroad, 'pos'); lfFilter.pos = sourceBroad.pos; end

	% Canonical channel list used for filters and all subsequent analyses
	chanList = lf.label;
	lfFilter.label = chanList;

	% Power per frequency band using the common filter
	for bi = 1:size(freqBands, 1)
		band = freqBands(bi, :);

		% Band-pass filter data for this band
		cfg = [];
		cfg.bpfilter = 'yes';
		cfg.bpfreq = band;
		cfg.demean = 'yes';
		dataBand = ft_preprocessing(cfg, data);

		% Covariance for this band
		cfg = [];
		cfg.covariance = 'yes';
		cfg.covariancewindow = 'all';
		tlBand = ft_timelockanalysis(cfg, dataBand);

		% Compute canonical channel selection order consistent with FieldTrip
		canonChan = ft_channelselection(chanList, tlBand.label);

		% Align band-limited data channel order to the canonical list
		cfg = [];
		cfg.channel = canonChan;
		tlBand = ft_selectdata(cfg, tlBand);

		% Ensure the sourcemodel and cfg use the same labels/order
		lfFilter.label = canonChan(:);
 
		% Apply common filter to band-limited data
		cfg = [];
		cfg.method = 'lcmv';
		cfg.keeptrials = 'no';
		cfg.lcmv.keepfilter = 'no';
		cfg.lcmv.lambda = '5%';
		cfg.lcmv.fixedori = 'yes';
		cfg.lcmv.projectnoise = 'yes';
		cfg.lcmv.weightnorm = 'arraygain';
		cfg.headmodel = headmodel;
		cfg.sourcemodel = lfFilter; % reuse common filters with labels
		cfg.grid = lfFilter;        % support older FT versions expecting cfg.grid
		if isstring(canonChan); canonChan = cellstr(canonChan); end
		cfg.channel = canonChan(:); % must match sourcemodel.label
		sourceBand = ft_sourceanalysis(cfg, tlBand);

		% Plot orthogonal slices at global maximum power
		cfgp = [];
		cfgp.method = 'ortho';
		cfgp.funparameter = 'pow';
		cfgp.location = 'max';
		cfgp.funcolorlim = 'maxabs';
		cfgp.maskparameter = 'pow';
		cfgp.camlight = 'no';
		ft_sourceplot(cfgp, sourceBand);
		title(sprintf('%s | %s [%d-%d] Hz', setBase, bandNames{bi}, band(1), band(2)));

		outPng = fullfile(outputDir, sprintf('%s_%s_%d-%dHz.png', setBase, bandNames{bi}, band(1), band(2)));
		set(gcf, 'Color', 'w');
		try
			exportgraphics(gcf, outPng, 'Resolution', 200);
		catch
			saveas(gcf, outPng);
		end
		close(gcf);
	end
end

disp(['Done. Figures written to: ' outputDir]);