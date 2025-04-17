% FND_DTI_proc.m / Nicolas Gninenko / Dec 2022
%
% Preprocessing script for DTI (DWI) data of the BioGen project
% Original data is located in
% X:\NRLK_FND\Experiments\Original Results\18 - BioGen\02_MRI
% for each participant (PXXX codes)
%

tmp_getLocalHost = char(java.net.InetAddress.getLocalHost.getHostName);
if ismac && strcmp(getenv('USER'),'nicogn')
    basepath_data = '/Volumes/Data/Nico/BioGen/02_MRI/';
    outputnii_dir = '/Volumes/Data/Nico/BioGen/DTIproc/';
    fs_reconall_dir = [outputnii_dir '_fs_reconall/'];
    dcm2nii_tmppath = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';
    tpm_spm12path = '/Users/nicogn/Documents/MATLAB/spm12/tpm/TPM.nii';
    addpath('/Users/nicogn/Documents/MATLAB/NIfTI_20140122/');
% elseif isunix && strcmp(tmp_getLocalHost(1:end-1),'miplabsrv')
%     basepath_data = [filesep 'media' filesep 'miplab-nas2' filesep 'Data' filesep 'Nico' ...
%         filesep 'BioGen' filesep '02_MRI' filesep];
%     outputnii_dir = [filesep 'media' filesep 'miplab-nas2' filesep 'Data' filesep 'Nico' ...
%         filesep 'BioGen' filesep 'DTIproc' filesep];
%     fs_reconall_dir = [outputnii_dir '_fs_reconall' filesep];
%     addpath([filesep 'media' filesep 'miplab-nas2' filesep 'NeuroTin' filesep 'MATLAB' filesep 'spm12' filesep]);
%     addpath([filesep 'media' filesep 'miplab-nas2' filesep 'NeuroTin' filesep 'MATLAB' filesep 'NIfTI_20140122' filesep]);
%     tpm_spm12path = [filesep 'media' filesep 'miplab-nas2' filesep 'NeuroTin' filesep ...
%         'MATLAB' filesep 'spm12' filesep 'tpm' filesep 'TPM.nii'];
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


%% Import DCM to NII for all DTI data

spm('defaults','fmri');
spm_jobman('initcfg');

skip_step = true;
if ~skip_step
    for p = 1:patients_nb
    %for p = 39:40 % P235:P236
        if ~exist([outputnii_dir patients_list{p,1} filesep],'dir')
            %mkdir([outputnii_dir patients_list{p,1} filesep]);
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

%% Separately import CoIFA via dcm2niix bc of DCM metadata problem there
%  only runs locally on OS X, because of dcm2niix version for now

skip_step = true;
if ~skip_step
    if ismac
        for p = 1:patients_nb %39:40
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
% TRACEW co-registered maps

skip_step = true;
if ~skip_step
    for p = 1:patients_nb %39:40%:patients_nb % 87 refers to P295, a problematic
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
    subfold_id = 1;
    for subfold = {'01_T1','02_FUP'}
        for ftype = {'SCAN_ADC','SCAN_FA','TRACEW'} % ColFA exlc. for now [not pre-proc]
            tmp_FOI_intersect = dir([outputnii_dir patients_list{2,1} filesep subfold{1} filesep 'DTI' filesep ...
                'DQS_' ftype{1} filesep 'wrs*.nii']); % we chose patients_list{2,1} because P117 has both T1 and FUP data
            tmp_FOI_intersect = load_untouch_nii([outputnii_dir patients_list{2,1} filesep subfold{1} filesep ...
                'DTI' filesep 'DQS_' ftype{1} filesep tmp_FOI_intersect(1).name]);
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

%% MRTRIX pipeline (depth 2/3, without tractography yet)
% [   BEFORE RUNNING THIS CODE IN MATLAB, MAKE SURE YOU HAVE EXPORTED THE
% ENV VARIABLE TO HAVE THE MRTRIX3 BIN DIRECTORY IN THE PATH (BASH) :
%
%   >> bash
%   >> :~$ PATH="/usr/local/mrtrix3/bin:${PATH}"
%   >> :~$ export PATH
%
% ...AND FSL PATH IS ALSO SET UP (all of the below in bash as well) :
%
%   >> :~$ FSLDIR=/usr/local/fsl
%   >> :~$ PATH=${FSLDIR}/bin:${PATH}
%   >> :~$ export FSLDIR PATH
%   >> :~$ . ${FSLDIR}/etc/fslconf/fsl.sh
%
% ...before launching MATLAB in the desired screen   ]

skip_step = true;
if ~skip_step
    subfold_id = 1;
    for subfold = {'01_T1','02_FUP'}
        for p = 1:patients_nb
            if patients_list{p,subfold_id} % i.e. if data exists
                tmp_wrkdir = [outputnii_dir patients_list{p,1} filesep subfold{1} ...
                    filesep 'DTI' filesep 'mrtrix_pipe' filesep];
                if exist(tmp_wrkdir,'dir'), system(['rm -r ' tmp_wrkdir]); end % tmp to clear all when problematic
                if ~exist(tmp_wrkdir,'dir'), mkdir(tmp_wrkdir); end
                if ~exist([tmp_wrkdir 'dwi.mif'],'file')
                    system(['mrconvert ' basepath_data patients_list{p,1} ...
                        filesep subfold{1} filesep 'DTI' filesep 'DQS_SCAN' ' -datatype float32 ' ...
                        tmp_wrkdir 'dwi.mif -force']);
                end
                cd(tmp_wrkdir);
                if ~exist([tmp_wrkdir 'preproc_dwi_mask.mif'],'file')
                    !dwi2mask dwi.mif preproc_dwi_mask.mif
                    !maskfilter preproc_dwi_mask.mif dilate preproc_dwi_mask.mif -npass 5 -force
                end
                if ~exist([tmp_wrkdir 'dwi_noise.mif'],'file')
                    !dwidenoise dwi.mif dwi_d.mif -noise dwi_noise.mif -mask preproc_dwi_mask.mif -force
                end
                if ~exist([tmp_wrkdir 'dwi_dg.mif'],'file')
                    !mrdegibbs dwi_d.mif dwi_dg.mif
                end
                if ~exist([tmp_wrkdir 'dwi_EC.mif'],'file')
                    !dwipreproc dwi_dg.mif dwi_EC.mif -rpe_none -pe_dir AP -nocleanup -force -eddy_options=--data_is_shelled
                    % [DR: disregard above if data is not shelled][acq in q space]
                end
                if ~exist([tmp_wrkdir 'b0AP.nii'],'file')
                    !dwiextract dwi_EC.mif -bzero b0AP.mif -force
                    !mrconvert b0AP.mif b0AP.nii -force
                end
                if ~exist([tmp_wrkdir 'b0AP_mask.mif'],'file')
                    !bet2 b0AP.nii b0AP_mask.nii.gz -m -f 0.1
                    !mrconvert b0AP_mask.nii.gz b0AP_mask.mif -force
                end
                if ~exist([tmp_wrkdir 'dwi_ECn.nii'],'file')
                    !dwinormalise individual dwi_EC.mif b0AP_mask.mif dwi_ECn.mif -force
                    !mrconvert dwi_ECn.mif dwi_ECn.nii -force
                end
                % T1 betting step, will be checked manually
                % [to improve betting, use custom cut...]
                if ~exist([tmp_wrkdir 'T1_betted.nii.gz'],'file')
                    T1_nii_to_bet = dir([outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        filesep 's*.nii']); % should return the T1's filename
                    if isempty(T1_nii_to_bet)
                        error(['Unable to find appropriate T1 file for ' patients_list{p,1}]);
                    end
                    T1_nii_to_bet = load_untouch_nii([outputnii_dir patients_list{p,1} filesep ...
                        subfold{1} filesep T1_nii_to_bet(1).name]);
                    T1_nii_to_bet.img(:,1:70,:) = 0; % hard-cropping is defined by used here, let it be attempted
                    %T1_nii_to_bet.img(:,1:84,:) = 0;
                    T1_nii_to_bet.img(1:136,70:106,:) = 0;
                    T1_nii_to_bet.img(1:70,106:136,:) = 0;
                    %tmp_data.img(1:60,:,:) = 0;
                    %tmp_data.img(1:113,100:113,:) = 0;
                    %tmp_data.img(1:113,100:113,:) = 0;
                    T1_nii_to_bet.fileprefix = [T1_nii_to_bet.fileprefix '_cut'];
                    save_untouch_nii(T1_nii_to_bet,[T1_nii_to_bet.fileprefix '.nii']);
                    %!bet T1.nii T1_betted.nii.gz -f 0.1 -B % original but use bet2 instead [faster, optimized]
                    system(['bet2 ' T1_nii_to_bet.fileprefix '.nii' ' ' ...
                        tmp_wrkdir 'T1_betted.nii.gz -f 0.25 -m']); % default params for -f and keep mask [-m]
                end
                if ~exist([tmp_wrkdir 'T1reo_betted.nii.gz'],'file')
                    !flirt -dof 6 -cost normmi -in T1_betted.nii.gz -ref b0AP.nii -omat T_fsl.txt
                    % [DR: maybe use b0AP_mask or so instead]
                    !transformconvert T_fsl.txt T1_betted.nii.gz b0AP.nii flirt_import T_T1toDWI.txt -force && rm T_fsl.txt
                    !mrtransform -linear T_T1toDWI.txt T1_betted.nii.gz T1reo_betted.nii.gz -force
                end
                if ~exist([tmp_wrkdir '5ttvis.nii'],'file')
                    !5ttgen fsl T1reo_betted.nii.gz 5ttseg.mif -premasked -force
                    !mrconvert 5ttseg.mif 5ttseg.nii -force
                    !5tt2vis 5ttseg.mif 5ttvis.nii -force
                    !mrconvert 5ttseg.mif -coord 3 2 -axes 0,1,2 - | mrthreshold - -abs 0.001 wm_seed.mif -force
                end
                % voxel level modeling [vlm]
                tmp_wrkdir2 = [tmp_wrkdir 'vlm' filesep];
                if ~exist(tmp_wrkdir2,'dir'), mkdir(tmp_wrkdir2); end
                if ~exist([tmp_wrkdir2 'dwi_ECn.mif'],'file') && ~exist([tmp_wrkdir2 'b0AP_mask.mif'],'file') && ...
                        ~exist([tmp_wrkdir2 'b0AP_mask.mif'],'file') && ~exist([tmp_wrkdir2 'T1reo_betted.nii.gz'],'file')
                    !cp b0AP_mask.mif dwi_ECn.mif 5ttseg.mif T1reo_betted.nii.gz vlm/
                end
                cd(tmp_wrkdir2);
                if ~exist([tmp_wrkdir2 'T1reo_betted.mif'],'file')
                    !mrconvert T1reo_betted.nii.gz T1reo_betted.mif -force
                end
                if ~exist([tmp_wrkdir2 'dt.mif'],'file')
                    !dwi2tensor -mask b0AP_mask.mif dwi_ECn.mif dt.mif -force
                end
                if ~exist([tmp_wrkdir2 'dt_ad.nii'],'file') && ~exist([tmp_wrkdir2 'dt_fa.nii'],'file') && ...
                        ~exist([tmp_wrkdir2 'dt_rd.nii'],'file') && ~exist([tmp_wrkdir2 'dt_ev.nii'],'file') && ...
                        ~exist([tmp_wrkdir2 'dt_adc.nii'],'file')
                    !tensor2metric dt.mif -adc dt_adc.nii -fa dt_fa.nii -ad dt_ad.nii -rd dt_rd.nii -vector dt_ev.nii -force
                end
                if ~exist([tmp_wrkdir2 'ss_wm.mif'],'file') && ~exist([tmp_wrkdir2 'ss_voxels.mif'],'file')
                    !dwi2response tax dwi_ECn.mif ss_wm.txt -voxels ss_voxels.mif -force
                    !dwi2fod csd -mask b0AP_mask.mif dwi_ECn.mif ss_wm.txt ss_wm.mif -force
                end
                %!mrconvert -coord 3 0 ss_wm.mif - | mrcat zeros.mif zeros.mif - ss_rgb.mif -force
            end
        end
        subfold_id = subfold_id + 1;
    end
end

%% MRTRIX tractography & SIFT (depth 3/3)

skip_step = true;
if ~skip_step
    subfold_id = 1;
    for subfold = {'01_T1','02_FUP'}
        for p = 1:patients_nb
            if patients_list{p,subfold_id}
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                        'mrtrix_pipe' filesep 'vlm' filesep 'tracks_ACT_seeddynamic.tck'],'file')
                    cd([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                        'mrtrix_pipe' filesep 'vlm' filesep]);
                    fprintf(['\nRunning tckgen for ' patients_list{p,1} ' (' subfold{1} ')...\n']);
                    system(['tckgen ss_wm.mif tracks_ACT_seeddynamic.tck ' ...
                        '-act 5ttseg.mif -seed_dynamic ss_wm.mif -backtrack ' ...
                        '-crop_at_gmwmi -select 100M -maxlength 250 -cutoff 0.06 -force']);
                end
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                        'mrtrix_pipe' filesep 'vlm' filesep 'tracks10M_SIFT.tck'],'file')
                    cd([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 'DTI' filesep ...
                        'mrtrix_pipe' filesep 'vlm' filesep]);
                    fprintf(['\nRunning tcksift [1] for ' patients_list{p,1} ' (' subfold{1} ')...\n']);
                    system(['tcksift tracks_ACT_seeddynamic.tck ss_wm.mif ' ...
                        'tracks10M_SIFT.tck -act 5ttseg.mif -term_number 10M ' ...
                        '-out_mu SIFT1_mu.txt -force']);
                end
            end
        end
        subfold_id = subfold_id + 1; % then proceed with 02_FUP data
    end
    %clear subfold_id;
end

%% FreeSurfer parcellation
% [   run all P1** first, without P16
% then separate code for P200-P249, without P28, and another one for
% P250-P299, and P300+   ]

skip_step = true;
if ~skip_step
    if ~strcmp(getenv('SUBJECTS_DIR'),fs_reconall_dir(1:end-1))
        fprintf('\nRe-setting env. variable SUBJECTS_DIR...\n');
        setenv('SUBJECTS_DIR',fs_reconall_dir(1:end-1));
    end
    subfold_id = 1;
    for subfold = {'01_T1','02_FUP'}
        for p = 1:patients_nb %[2:12 14:20] % only P1* without P16
            if patients_list{p,subfold_id} % i.e. if data exists
                if ~exist([fs_reconall_dir subfold{1} filesep patients_list{p,1} filesep],'dir')
                    cd([outputnii_dir patients_list{p,1} filesep subfold{1} filesep]);
                    tmp_T1_input = dir([outputnii_dir patients_list{p,1} filesep subfold{1} filesep 's*01.nii']);
                    if isempty(tmp_T1_input)
                        error(['Unable to find T1 input file for ' patients_list{p,1} ' (' subfold{1} ').']);
                    end
                    fprintf(['\n\nRunning recon-all for ' patients_list{p,1} ' (' subfold{1} ')...\n\n']);
                    system(['recon-all -autorecon-all -i ' tmp_T1_input(1).name ' -s ' patients_list{p,1}]);
                end
            end
        end
        subfold_id = subfold_id + 1; % then proceed with 02_FUP data...
    end
    %clear subfold_id;
end

%% Connectomics
%   Here we will use the parcellation created above (from FreeSurfer's
%   recon-all), namely the Desikan-Killiany one (aparc+aseg.mgz), to
%   create the connectomes (other parcellations can be added later on)

skip_step = true;
if ~skip_step
    subfold_id = 1;
    tmp_vlm_subpath = [filesep 'DTI' filesep 'mrtrix_pipe' filesep 'vlm' filesep];
    for subfold = {'01_T1','02_FUP'} % [split 02_FUP on srv3 and 01_T1 on srv4]
        for p = 1:length(patients_list) % run this for all patients
            if patients_list{p,subfold_id} % i.e. if data exists
                cd([outputnii_dir patients_list{p,1} filesep subfold{1} tmp_vlm_subpath]);
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'parc_fs.mif'],'file')
                    % first, labelconvert the parcellation into mrtrix3's convention
                    system(['labelconvert ' fs_reconall_dir subfold{1} filesep patients_list{p,1} filesep ...
                        'mri' filesep 'aparc+aseg.mgz ' outputnii_dir 'FreeSurferColorLUT.txt ' ...
                        outputnii_dir 'fs_default.txt ' outputnii_dir patients_list{p,1} filesep ...
                        subfold{1} tmp_vlm_subpath 'parc_fs.mif -force']);
                end
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'parc_fs_reo.mif'],'file')
                    % because mrtransform does not like the
                    % aparc+aseg.mgz directly, we reorient rather the
                    % resulting parcellation file (.mif)
                    system(['mrtransform -linear ' outputnii_dir patients_list{p,1} filesep subfold{1} filesep ...
                        'DTI' filesep 'mrtrix_pipe' filesep 'T_T1toDWI.txt ' outputnii_dir patients_list{p,1} ...
                        filesep subfold{1} tmp_vlm_subpath 'parc_fs.mif ' outputnii_dir patients_list{p,1} filesep ...
                        subfold{1} tmp_vlm_subpath 'parc_fs_reo.mif -force']);
                end
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'connectome.csv'],'file')
                    % generating the basic connectome (number of crossing
                    % tracks w.r.t. each ROI as per the SIFTed file (10M tracks))
                    fprintf(['\nGenerating connectome.csv for ' patients_list{p,1} ' (' subfold{1} ')...\n']);
                    system(['tck2connectome ' outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        tmp_vlm_subpath 'tracks10M_SIFT.tck ' outputnii_dir patients_list{p,1} ...
                        filesep subfold{1} tmp_vlm_subpath 'parc_fs_reo.mif ' outputnii_dir ...
                        patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'connectome.csv']);
                end
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'connectome_ml.csv'],'file')
                    % generate mean length connectome
                    fprintf(['\nGenerating connectome_ml.csv for ' patients_list{p,1} ' (' subfold{1} ')...\n']);
                    system(['tck2connectome ' outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        tmp_vlm_subpath 'tracks10M_SIFT.tck ' outputnii_dir patients_list{p,1} ...
                        filesep subfold{1} tmp_vlm_subpath 'parc_fs_reo.mif ' outputnii_dir ...
                        patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'connectome_ml.csv ' ...
                        '-scale_length -stat_edge mean']);
                end
                if ~exist([outputnii_dir patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'connectome_fa.csv'],'file')
                    % generate FA connectome (3 steps)
                    system(['mrconvert ' outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        tmp_vlm_subpath 'dt_fa.nii ' outputnii_dir patients_list{p,1} ...
                        filesep subfold{1} tmp_vlm_subpath 'dt_fa.mif']);
                    system(['tcksample ' outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        tmp_vlm_subpath 'tracks10M_SIFT.tck ' outputnii_dir patients_list{p,1} ...
                        filesep subfold{1} tmp_vlm_subpath 'dt_fa.mif ' outputnii_dir ...
                        patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'tck_meanFA.txt -stat_tck mean']);
                    fprintf(['\nGenerating connectome_fa.csv for ' patients_list{p,1} ' (' subfold{1} ')...\n']);
                    system(['tck2connectome ' outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        tmp_vlm_subpath 'tracks10M_SIFT.tck ' outputnii_dir patients_list{p,1} ...
                        filesep subfold{1} tmp_vlm_subpath 'parc_fs_reo.mif ' outputnii_dir ...
                        patients_list{p,1} filesep subfold{1} tmp_vlm_subpath 'connectome_fa.csv ' ...
                        '-scale_file ' outputnii_dir patients_list{p,1} filesep subfold{1} ...
                        tmp_vlm_subpath 'tck_meanFA.txt -stat_edge mean']);
                end
            end
        end
        subfold_id = subfold_id + 1; % then proceed with 02_FUP data
    end
end



