% fnd_DTI_proc.m / Nicolas Gninenko / Dec 2022
%
% Preprocessing script for DTI data of the BioGen project
% Original data is located in
% X:\NRLK_FND\Experiments\Original Results\18 - BioGen\02_MRI
% for each participant (PXXX codes)
%

tmp_getLocalHost = char(java.net.InetAddress.getLocalHost.getHostName);
if ismac && strcmp(getenv('USER'),'nicogn')
    basepath_data = '/Volumes/Data/Nico/BioGen/02_MRI/';
    outputnii_dir = '/Volumes/Data/Nico/BioGen/DTIproc/';
    dcm2nii_tmppath = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';
    tpm_spm12path = '/Users/nicogn/Documents/MATLAB/spm12/tpm/TPM.nii';
    addpath('/Users/nicogn/Documents/MATLAB/NIfTI_20140122/');
elseif isunix && strcmp(tmp_getLocalHost(1:end-1),'miplabsrv')
    basepath_data = [filesep 'media' filesep 'miplab-nas2' filesep 'Data' filesep 'Nico' ...
        filesep 'BioGen' filesep '02_MRI' filesep];
    outputnii_dir = [filesep 'media' filesep 'miplab-nas2' filesep 'Data' filesep 'Nico' ...
        filesep 'BioGen' filesep 'DTIproc' filesep];
    addpath([filesep 'media' filesep 'miplab-nas2' filesep 'NeuroTin' filesep 'MATLAB' filesep 'spm12' filesep]);
    addpath([filesep 'media' filesep 'miplab-nas2' filesep 'NeuroTin' filesep 'MATLAB' filesep 'NIfTI_20140122' filesep]);
    tpm_spm12path = [filesep 'media' filesep 'miplab-nas2' filesep 'NeuroTin' filesep ...
        'MATLAB' filesep 'spm12' filesep 'tpm' filesep 'TPM.nii'];
end

patients_list = dir([basepath_data 'P*']);
patients_list = {patients_list.name}';
patients_nb   = length(patients_list);
patients_list(:,2:3) = {false};

for j = 1:patients_nb
    tmp_dti_dir1 = [basepath_data patients_list{j,1} filesep '01_T1' filesep 'DTI' filesep 'DQS_SCAN'];
    tmp_dti_dir2 = [basepath_data patients_list{j,1} filesep '02_FUP' filesep 'DTI' filesep 'DQS_SCAN'];
    if exist(tmp_dti_dir1,'dir'), if ~isempty(tmp_dti_dir1), patients_list{j,2} = true; end; end
    if exist(tmp_dti_dir2,'dir'), if ~isempty(tmp_dti_dir2), patients_list{j,3} = true; end; end
end


%% import DCM to NII for all DTI data

spm('defaults','fmri');
spm_jobman('initcfg');

skip_step = true;
if ~skip_step
    for p = 1:patients_nb
        if ~exist([outputnii_dir patients_list{p,1} filesep],'dir')
            mkdir([outputnii_dir patients_list{p,1} filesep]);
            clear matlabbatch cc;
            cc = 1;
            if patients_list{p,2} % if true, the T01 (first time point) DTI exists
                mkdir([outputnii_dir patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep 'DQS_SCAN']);
                mkdir([outputnii_dir patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep 'DQS_SCAN_ADC']);
                mkdir([outputnii_dir patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep 'DQS_SCAN_CoIFA']);
                mkdir([outputnii_dir patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep 'DQS_SCAN_FA']);
                mkdir([outputnii_dir patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep 'DQS_TRACEW']);
                if exist([basepath_data patients_list{p,1} filesep '01_T1' ...
                        filesep 'DTI' filesep 'DQS_SCAN_TRACEW'],'dir')
                    movefile([basepath_data patients_list{p,1} filesep '01_T1' ...
                        filesep 'DTI' filesep 'DQS_SCAN_TRACEW'],...
                        [basepath_data patients_list{p,1} filesep '01_T1' ...
                        filesep 'DTI' filesep 'DQS_TRACEW']);
                    pause(1);
                end
                for f_folders = {'DQS_SCAN','DQS_SCAN_ADC','DQS_SCAN_FA','DQS_TRACEW'} % 'DQS_SCAN_CoIFA'
                    tmp_flist = dir([basepath_data patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep ...
                        f_folders{1} filesep '*.dcm']);
                    matlabbatch{cc}.spm.util.import.dicom.data = cell(length(tmp_flist),1);
                    for f = 1:length(tmp_flist)
                        matlabbatch{cc}.spm.util.import.dicom.data{f,1} = [basepath_data patients_list{p,1} filesep ...
                            '01_T1' filesep 'DTI' filesep f_folders{1} filesep tmp_flist(f).name];
                    end
                    matlabbatch{cc}.spm.util.import.dicom.root = 'flat';
                    matlabbatch{cc}.spm.util.import.dicom.outdir = {[outputnii_dir patients_list{p,1} filesep ...
                        '01_T1' filesep 'DTI' filesep f_folders{1}]};
                    matlabbatch{cc}.spm.util.import.dicom.protfilter = '.*';
                    matlabbatch{cc}.spm.util.import.dicom.convopts.format = 'nii';
                    matlabbatch{cc}.spm.util.import.dicom.convopts.meta = 0;
                    matlabbatch{cc}.spm.util.import.dicom.convopts.icedims = 0;
                    cc = cc+1;
                end
            end
            if patients_list{p,3} % if true, the F02 (follow-up) DTI exists
                mkdir([outputnii_dir patients_list{p,1} filesep '02_FUP' filesep 'DTI' filesep 'DQS_SCAN']);
                mkdir([outputnii_dir patients_list{p,1} filesep '02_FUP' filesep 'DTI' filesep 'DQS_SCAN_ADC']);
                mkdir([outputnii_dir patients_list{p,1} filesep '02_FUP' filesep 'DTI' filesep 'DQS_SCAN_CoIFA']);
                mkdir([outputnii_dir patients_list{p,1} filesep '02_FUP' filesep 'DTI' filesep 'DQS_SCAN_FA']);
                mkdir([outputnii_dir patients_list{p,1} filesep '02_FUP' filesep 'DTI' filesep 'DQS_TRACEW']);
                if exist([basepath_data patients_list{p,1} filesep '02_FUP' ...
                        filesep 'DTI' filesep 'DQS_SCAN_TRACEW'],'dir')
                    movefile([basepath_data patients_list{p,1} filesep '02_FUP' ...
                        filesep 'DTI' filesep 'DQS_SCAN_TRACEW'],...
                        [basepath_data patients_list{p,1} filesep '02_FUP' ...
                        filesep 'DTI' filesep 'DQS_TRACEW']);
                    pause(1);
                end
                for f_folders = {'DQS_SCAN','DQS_SCAN_ADC','DQS_SCAN_FA','DQS_TRACEW'} % 'DQS_SCAN_CoIFA'
                    tmp_flist = dir([basepath_data patients_list{p,1} filesep '02_FUP' filesep 'DTI' filesep ...
                        f_folders{1} filesep '*.dcm']);
                    matlabbatch{cc}.spm.util.import.dicom.data = cell(length(tmp_flist),1);
                    for f = 1:length(tmp_flist)
                        matlabbatch{cc}.spm.util.import.dicom.data{f,1} = [basepath_data patients_list{p,1} filesep ...
                            '02_FUP' filesep 'DTI' filesep f_folders{1} filesep tmp_flist(f).name];
                    end
                    matlabbatch{cc}.spm.util.import.dicom.root = 'flat';
                    matlabbatch{cc}.spm.util.import.dicom.outdir = {[outputnii_dir patients_list{p,1} filesep ...
                        '02_FUP' filesep 'DTI' filesep f_folders{1}]};
                    matlabbatch{cc}.spm.util.import.dicom.protfilter = '.*';
                    matlabbatch{cc}.spm.util.import.dicom.convopts.format = 'nii';
                    matlabbatch{cc}.spm.util.import.dicom.convopts.meta = 0;
                    matlabbatch{cc}.spm.util.import.dicom.convopts.icedims = 0;
                    cc = cc+1;
                end
            end
            if exist('matlabbatch','var')
                spm_jobman('run',matlabbatch);
            end
        end
    end
end

%% separately import CoIFA via dcm2niix bc of dcm metadata problem there
% only runs locally on OS X, because of dcm2niix version for now

skip_step = true;
if ~skip_step
    if ismac
        for p = 80:patients_nb
            if patients_list{p,2}
                system([dcm2nii_tmppath ' -f "%f_%p_%t_%s" -p y -z n -o ' ...
                    '"' outputnii_dir patients_list{p,1} filesep '01_T1' ...
                    filesep 'DTI' filesep 'DQS_SCAN_CoIFA" ' ...
                    '"' basepath_data patients_list{p,1} filesep '01_T1' ...
                    filesep 'DTI' filesep 'DQS_SCAN_CoIFA"']);
            end
            if patients_list{p,3}
                system([dcm2nii_tmppath ' -f "%f_%p_%t_%s" -p y -z n -o ' ...
                    '"' outputnii_dir patients_list{p,1} filesep '02_FUP' ...
                    filesep 'DTI' filesep 'DQS_SCAN_CoIFA" ' ...
                    '"' basepath_data patients_list{p,1} filesep '02_FUP' ...
                    filesep 'DTI' filesep 'DQS_SCAN_CoIFA"']);
            end
        end
    end
end

%% (Re-)import T1 structurals for later-on reg/coreg.

skip_step = true;
if ~skip_step
for p = 1:patients_nb
    clear matlabbatch cc;
    cc = 1;
    if patients_list{p,2} && isempty(dir([outputnii_dir patients_list{p,1} filesep '01_T1' filesep '*.nii']))
        tmp_flist = dir([basepath_data patients_list{p,1} filesep '01_T1' filesep ...
            'structurals' filesep '*.dcm']);
        matlabbatch{cc}.spm.util.import.dicom.data = cell(length(tmp_flist),1);
        for f = 1:length(tmp_flist)
            matlabbatch{cc}.spm.util.import.dicom.data{f,1} = [basepath_data patients_list{p,1} filesep ...
                '01_T1' filesep 'structurals' filesep tmp_flist(f).name];
        end
        matlabbatch{cc}.spm.util.import.dicom.root = 'flat';
        matlabbatch{cc}.spm.util.import.dicom.outdir = {[outputnii_dir patients_list{p,1} filesep '01_T1']};
        matlabbatch{cc}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{cc}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{cc}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{cc}.spm.util.import.dicom.convopts.icedims = 0;
        cc = cc+1;
    end
    if patients_list{p,3} && isempty(dir([outputnii_dir patients_list{p,1} filesep '02_FUP' filesep '*.nii']))
        tmp_flist = dir([basepath_data patients_list{p,1} filesep '02_FUP' filesep ...
            'structurals' filesep '*.dcm']);
        matlabbatch{cc}.spm.util.import.dicom.data = cell(length(tmp_flist),1);
        for f = 1:length(tmp_flist)
            matlabbatch{cc}.spm.util.import.dicom.data{f,1} = [basepath_data patients_list{p,1} filesep ...
                '02_FUP' filesep 'structurals' filesep tmp_flist(f).name];
        end
        matlabbatch{cc}.spm.util.import.dicom.root = 'flat';
        matlabbatch{cc}.spm.util.import.dicom.outdir = {[outputnii_dir patients_list{p,1} filesep '02_FUP']};
        matlabbatch{cc}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{cc}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{cc}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{cc}.spm.util.import.dicom.convopts.icedims = 0;
        cc = cc+1;
    end
    if exist('matlabbatch','var')
        spm_jobman('run',matlabbatch);
    end
end
end

%% Co-register all AFC(, ColFA), FA, and TRACEW maps
% to respective anatomical (imported T1 above) and
% perform subsequent normalization of T1 into MNI space,
% and then finally apply the y_*.nii transform from
% resulting mapping to all the AFC(, ColFA), FA, and
% TRACEW co-registered maps, to have them all in the
% same space (MNI standard) + smooth all output files
% for further group analysis (next steps below)

skip_step = true;
if ~skip_step
    for p = 88:patients_nb % check later on if P295 (p = 87) has an updated DQS_SCAN_ADC
        subfold_id = 2; % 2 for 01_T1 and 3 for 02_FUP
        for subfold = {'01_T1','02_FUP'}
            if patients_list{p,subfold_id} % if true, the T01 (first time point) DTI exists & also for 02_FUP
                clear matlabbatch tmp_file* cc;
                cc = 1;
                tmp_file_T1 = dir([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 's*.nii']);
                matlabbatch{cc}.spm.spatial.coreg.estwrite.ref = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep tmp_file_T1(1).name ',1']};
                tmp_file = dir([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_SCAN_ADC' filesep 's*.nii']);
                matlabbatch{cc}.spm.spatial.coreg.estwrite.source = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_SCAN_ADC' filesep tmp_file(1).name ',1']};
                matlabbatch{cc}.spm.spatial.coreg.estwrite.other = {''};
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 ...
                    0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                cc = cc+1;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.ref = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep tmp_file_T1(1).name ',1']};
                tmp_file = dir([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_SCAN_FA' filesep 's*.nii']);
                matlabbatch{cc}.spm.spatial.coreg.estwrite.source = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_SCAN_FA' filesep tmp_file(1).name ',1']};
                matlabbatch{cc}.spm.spatial.coreg.estwrite.other = {''};
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 ...
                    0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                cc = cc+1;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.ref = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep tmp_file_T1(1).name ',1']};
                tmp_file = dir([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_TRACEW' filesep 's*.nii']);
                matlabbatch{cc}.spm.spatial.coreg.estwrite.source = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_TRACEW' filesep tmp_file(1).name ',1']};
                matlabbatch{cc}.spm.spatial.coreg.estwrite.other = {
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_TRACEW' filesep tmp_file(2).name ',1']
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                    'DQS_TRACEW' filesep tmp_file(3).name ',1']
                    };
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 ...
                    0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                matlabbatch{cc}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                cc = cc+1;
                matlabbatch{cc}.spm.spatial.normalise.estwrite.subj.vol = {...
                    [outputnii_dir patients_list{p,1} filesep subfold{1} filesep tmp_file_T1(1).name ',1']};
                matlabbatch{cc}.spm.spatial.normalise.estwrite.subj.resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', ...
                    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
                matlabbatch{cc}.spm.spatial.normalise.estwrite.subj.resample(2) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', ...
                    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
                matlabbatch{cc}.spm.spatial.normalise.estwrite.subj.resample(3) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', ...
                    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.tpm = {tpm_spm12path};
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
                matlabbatch{cc}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
                matlabbatch{cc}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                    78 76 85];
                matlabbatch{cc}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
                matlabbatch{cc}.spm.spatial.normalise.estwrite.woptions.interp = 4;
                matlabbatch{cc}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
                cc = cc+1;
                matlabbatch{cc}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', ...
                    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
                matlabbatch{cc}.spm.spatial.smooth.fwhm = [4 4 4];
                matlabbatch{cc}.spm.spatial.smooth.dtype = 0;
                matlabbatch{cc}.spm.spatial.smooth.im = 0;
                matlabbatch{cc}.spm.spatial.smooth.prefix = 's';
                %
                spm_jobman('run',matlabbatch);
            end
            subfold_id = subfold_id + 1;
        end
    end
end
clear matlabbatch tmp_file* subfold* cc;

%% Create intersect mask for all P** (patients)
% for normalized ADC, FA, and TRACEW maps

skip_step = true;
if ~skip_step
    subfold_id = 2;
    for subfold = {'01_T1','02_FUP'}
        for ftype = {'SCAN_ADC','SCAN_FA','TRACEW'} % ColFA exlc. for now [not pre-proc]
            tmp_FOI_intersect = dir([outputnii_dir patients_list{2,1} filesep subfold{1} filesep 'DTI' filesep ...
                'DQS_' ftype{1} filesep 'wrs*.nii']); % we chose patients_list{2,1} because P117 has both T1 and FUP data
            tmp_FOI_intersect = load_untouch_nii([outputnii_dir patients_list{2,1} filesep subfold{1} filesep ...
                'DTI' filesep 'DQS_' ftype{1} filesep tmp_FOI_intersect(1).name]);
            %
            % P295 has no ADC (for now...) -> skip it here temporarily
            if isequal(subfold{1},'01_T1') && isequal(ftype{1},'SCAN_ADC'), tmp_patients_p = [1:86 88:patients_nb];
            else, tmp_patients_p = 1:patients_nb;
            end
            for p = tmp_patients_p
                if patients_list{p,subfold_id}
                    tmp_file = dir([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                        'DQS_' ftype{1} filesep 'wrs*.nii']);
                    tmp_file = load_untouch_nii([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                        'DQS_' ftype{1} filesep tmp_file(1).name]);
                    tmp_FOI_intersect.img = reshape(max(min(tmp_FOI_intersect.img(:),tmp_file.img(:)),0),...
                        size(tmp_FOI_intersect.img));
                end
            end
            tmp_FOI_intersect.fileprefix = ['intersect_' num2str(patients_nb) 'P_' ftype{1} '_' subfold{1}];
            save_untouch_nii(tmp_FOI_intersect,[outputnii_dir tmp_FOI_intersect.fileprefix '.nii']);
        end
        subfold_id = subfold_id + 1;
    end
end

%% mtrix separate file



