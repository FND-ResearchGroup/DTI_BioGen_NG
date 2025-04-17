% FND_DTI_connectomics.m / Nicolas Gninenko / Oct 2023
%
% Script to compute graph metrics and other stats & figure plots
% on obtained connectomes from the DWI data set of the BioGen project.
%
% Original data is located in
% ...\NRLK_FND\Experiments\Original Results\18 - BioGen\02_MRI
% for each participant (PXXX codes).
%
% Preprocessed data (& connectomes) are located in
% ./media/[!#path hidden]
%
% Contact: nicolas.gninenko@gmail.com
%
%%%
%
%   Plots of outputs:
%       1 - Mean connectomes for HC & FND populations (9 subplots in 1 Fig):
%           - Mean number of fibers crossing each ROI (column 1)
%           - Mean length of each bundle between each ROI (column 2)
%           - Mean FA between each ROI (column 3)
%
%   Usage:
%       If you want to see plots, modify the following vector 'skip_plots'
%       with the numbers above being 1 to see the corresponding plot
%       or 0 to hide it. Example : skip_plots = [1 0 0 ... 0] to see only
%       the first plot (defined above as the 9 subplots Figure).
%            [Fig1 Fig2 Fig3 ...]
% show_plots = [   1    1    1   1]; % feature dropped
%
%%%

% Init
tmp_getLocalHost = char(java.net.InetAddress.getLocalHost.getHostName);
if ismac && strcmp(getenv('USER'),'nicogn')
    basepath_data = '/Volumes/Data/Nico/BioGen/02_MRI/'; % [mounted MIP:Lab NAS - there is no participant data locally]
    proc_dir = '/Volumes/Data/Nico/BioGen/DTIproc/';
    fs_reconall_dir = [proc_dir '_fs_reconall/'];
    %dcm2nii_tmppath = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';
    %tpm_spm12path = '/Users/nicogn/Documents/MATLAB/spm12/tpm/TPM.nii';
    addpath('/Users/nicogn/Documents/MATLAB/NIfTI_20140122/');
    addpath('/Users/nicogn/Documents/MATLAB/exportfig');
    %addpath(genpath('/Users/nicogn/Documents/Bern buffer/Nicolas/Code/toolboxes/ENIGMA-2.0.0'));
    addpath('/Users/nicogn/Documents/Bern buffer/Nicolas/Code/toolboxes/biChordChart');
    %addpath('/Users/nicogn/Documents/MATLAB/');
    DWIcode_path = '/Users/nicogn/Documents/Bern buffer/Nicolas/Code/DWI/';
    addpath(genpath(DWIcode_path));
    cd(DWIcode_path);
    export_figs_path = [DWIcode_path 'Figures/'];
else
    % you may add your own path definitions here
end

% Saved [final] patients list, if changes need to use
% FND_DTI_proc.m to pre-process any additional/modified data
if exist([proc_dir 'patients_list.mat'],'file')
    patients_list = load([proc_dir 'patients_list.mat']); patients_list = patients_list.patients_list;
    patients_nb = length(patients_list);
end



%% Load all available connectomes
% You may reimport them as well by unskiping this step below,
% but this step can take a while...

skip_step = true;
if ~skip_step
    connectomes = cell(patients_nb,2); % second dimension is for 01_T1 & 02_FUP
    connectomes_ml = cell(patients_nb,2); connectomes_fa = cell(patients_nb,2);
    for p = 1:patients_nb
        tmp_path1 = [proc_dir patients_list{p,1} filesep '01_T1' filesep 'DTI' filesep ...
            'mrtrix_pipe' filesep 'vlm' filesep 'connectome_fa.csv'];
        if exist(tmp_path1,'file') % check last step of pre-proc
            connectomes{p,1} = importdata(strrep(tmp_path1,'_fa.csv','.csv'));
            connectomes_ml{p,1} = importdata(strrep(tmp_path1,'_fa.csv','_ml.csv'));
            connectomes_fa{p,1} = importdata(tmp_path1);
        end
        if exist(strrep(tmp_path1,'01_T1','02_FUP'),'file')
            connectomes{p,2} = importdata(strrep(strrep(tmp_path1,'_fa.csv','.csv'),'01_T1','02_FUP'));
            connectomes_ml{p,2} = importdata(strrep(strrep(tmp_path1,'_fa.csv','_ml.csv'),'01_T1','02_FUP'));
            connectomes_fa{p,2} = importdata(strrep(tmp_path1,'01_T1','02_FUP'));
        end
    end
else % local for faster loading
    load([DWIcode_path 'FND_DTI_connectomes.mat']);
    % contains variables above (connectomes, connectomes_ml, connectomes_fa)
    if ~exist('patients_nb','var'), patients_nb = length(patients_list); end
end



%% Load covariates [and reshuffle patients list]

skip_step = false;
if ~skip_step
    patients_datatable = readtable([DWIcode_path '230912 BioGen_Covariates_Nico SW.xlsx']);
    patients_datatable = sortrows(patients_datatable,1); % sort according to patient ID (increasing)
    [~,j] = sort(cell2mat(cellfun(@(x) str2double(x(2:end)), patients_list(:,1), 'UniformOutput', false)));
    patients_list = patients_list(j,:); % re-sort the list according to the Excel list above
    connectomes = connectomes(j,:); connectomes_ml = connectomes_ml(j,:); connectomes_fa = connectomes_fa(j,:);
    if ~isequal(cell2mat(cellfun(@(x) str2double(x(2:end)), patients_list(:,1), 'UniformOutput', false)),...
            table2array(patients_datatable(:,1)))
        error('Imported Excel patients table and processed patient data do not match...');
    end % check that we have the same patients list
end



%% Mean connectomes for HCs and FNDs at T01 and FUP

connectomes_HC_ind  = cell(2,2); % the second column is to store the contributing number of connectomes
connectomes_FND_ind = cell(2,2); % the first column is to store the connectomes' indices
                                  % i.e., first row is 01_T1 and second row 02_FUP
connectomes_HC_ind(:,:) = {0}; connectomes_FND_ind(:,:) = {0};

for j = 1:length(patients_list)
    if strcmp(patients_datatable{j,2}{1,1},'HC') && patients_list{j,2} % patient is HC and has T01 data
        connectomes_HC_ind{1,1} = [connectomes_HC_ind{1,1};j];
    end
    if strcmp(patients_datatable{j,2}{1,1},'HC') && patients_list{j,3} % patient is HC and has FUP data
        connectomes_HC_ind{2,1} = [connectomes_HC_ind{2,1};j]; % [we do not have any in fact]
    end
    if strcmp(patients_datatable{j,2}{1,1},'FND') && patients_list{j,2} % patient is FND and has T01 data
        connectomes_FND_ind{1,1} = [connectomes_FND_ind{1,1};j];
    end
    if strcmp(patients_datatable{j,2}{1,1},'FND') && patients_list{j,3} % patient is FND and has FUP data
        connectomes_FND_ind{2,1} = [connectomes_FND_ind{2,1};j];
    end
end

connectomes_HC_ind{1,1} = connectomes_HC_ind{1,1}(2:end); connectomes_HC_ind{2,1} = connectomes_HC_ind{2,1}(2:end);
connectomes_FND_ind{1,1} = connectomes_FND_ind{1,1}(2:end); connectomes_FND_ind{2,1} = connectomes_FND_ind{2,1}(2:end);
connectomes_HC_ind{1,2} = length(connectomes_HC_ind{1,1}); connectomes_FND_ind{1,2} = length(connectomes_FND_ind{1,1});
connectomes_HC_ind{2,2} = length(connectomes_HC_ind{2,1}); connectomes_FND_ind{2,2} = length(connectomes_FND_ind{2,1});

mean_connectomes_HC  = repmat({zeros(size(connectomes{1,1}))},1,3); % to store mean connectome, _ml, and _fa
mean_connectomes_FND = repmat(mean_connectomes_HC,2,1);

for j = connectomes_HC_ind{1,1}'
    mean_connectomes_HC{1,1} = mean_connectomes_HC{1,1} + connectomes{j,1};
    mean_connectomes_HC{1,2} = mean_connectomes_HC{1,2} + connectomes_ml{j,1};
    mean_connectomes_HC{1,3} = mean_connectomes_HC{1,3} + connectomes_fa{j,1};
end
for j = connectomes_FND_ind{1,1}'
    mean_connectomes_FND{1,1} = mean_connectomes_FND{1,1} + connectomes{j,1};
    mean_connectomes_FND{1,2} = mean_connectomes_FND{1,2} + connectomes_ml{j,1};
    mean_connectomes_FND{1,3} = mean_connectomes_FND{1,3} + connectomes_fa{j,1};
end
for j = connectomes_FND_ind{2,1}'
    mean_connectomes_FND{2,1} = mean_connectomes_FND{2,1} + connectomes{j,1};
    mean_connectomes_FND{2,2} = mean_connectomes_FND{2,2} + connectomes_ml{j,1};
    mean_connectomes_FND{2,3} = mean_connectomes_FND{2,3} + connectomes_fa{j,1};
end

for j = 1:3 % re-divide by the number of contributing connectomes in each case
    mean_connectomes_HC{1,j} = mean_connectomes_HC{1,j}/connectomes_HC_ind{1,2};
    % symmetrize
    mean_connectomes_HC{1,j} = mean_connectomes_HC{1,j} + transpose(triu(mean_connectomes_HC{1,j},1));
    for k = 1:2 % because we have T01 and FUP for FND (none of FUP for HC yet)
        mean_connectomes_FND{k,j} = mean_connectomes_FND{k,j}/connectomes_FND_ind{k,2};
        % symmetrize
        mean_connectomes_FND{k,j} = mean_connectomes_FND{k,j} + transpose(triu(mean_connectomes_FND{k,j},1));
    end
end



%% Plot mean connectomes (number of streamlines, mean length, and FA)
%  for HC & FND at timepoints 01_T1 [T1] and 02_FUP [FUP]

if show_plots(1)
    close all;
    tmp_f = figure('Position',[60 60 1600 900]);
    subplot(331);
    imagesc(mean_connectomes_HC{1,1});
    axis square; title('Mean # of streamlines','FontWeight','normal','Interpreter','none');
    colorbar; caxis([0 1800]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    ylabel('HC [T1]','Interpreter','none','Rotation',0,'FontSize',20);
    tmp_f.CurrentAxes.YAxis.Label.Position(1) = tmp_f.CurrentAxes.YAxis.Label.Position(1)-25;
    tmp_f.CurrentAxes.Title.FontSize = 20;
    tmp_f.CurrentAxes.Title.Position(2) = tmp_f.CurrentAxes.Title.Position(2)-10;
    subplot(332);
    imagesc(mean_connectomes_HC{1,2});
    axis square; title('Mean fibers'' length','FontWeight','normal','Interpreter','none');
    colorbar; caxis([0 180]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    tmp_f.CurrentAxes.Title.FontSize = 20;
    tmp_f.CurrentAxes.Title.Position(2) = tmp_f.CurrentAxes.Title.Position(2)-10;
    subplot(333);
    imagesc(mean_connectomes_HC{1,3});
    axis square; title('Mean FA','FontWeight','normal','Interpreter','none');
    colorbar; caxis([0 .75]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    tmp_f.CurrentAxes.Title.FontSize = 20;
    tmp_f.CurrentAxes.Title.Position(2) = tmp_f.CurrentAxes.Title.Position(2)-10;
    subplot(334);
    imagesc(mean_connectomes_FND{1,1});
    axis square; colorbar; caxis([0 1800]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    ylabel('FND [T1]','Interpreter','none','Rotation',0,'FontSize',20);
    tmp_f.CurrentAxes.YAxis.Label.Position(1) = tmp_f.CurrentAxes.YAxis.Label.Position(1)-27.5;
    subplot(335);
    imagesc(mean_connectomes_FND{1,2});
    axis square; colorbar; caxis([0 180]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    subplot(336);
    imagesc(mean_connectomes_FND{1,3});
    axis square; colorbar; caxis([0 .75]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    subplot(337);
    imagesc(mean_connectomes_FND{2,1});
    axis square; colorbar; caxis([0 1800]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    ylabel('FND [FUP]','Interpreter','none','Rotation',0,'FontSize',20);
    tmp_f.CurrentAxes.YAxis.Label.Position(1) = tmp_f.CurrentAxes.YAxis.Label.Position(1)-31;
    subplot(338);
    imagesc(mean_connectomes_FND{2,2});
    axis square; colorbar; caxis([0 180]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
    subplot(339);
    imagesc(mean_connectomes_FND{2,3});
    axis square; colorbar; caxis([0 .75]); xticks(10:10:80); yticks(10:10:80); tmp_f.CurrentAxes.FontSize = 12;
end
if ~exist([export_figs_path 'Figure1.pdf'],'file')
    %export_fig([export_figs_path 'Figure1.svg'],tmp_f);
    %exportgraphics(tmp_f,[export_figs_path 'Figure1.pdf'],...
    %    'BackgroundColor','none','ContentType','vector');
end



%% Plot differences in all connectomes (& _ml, _fa) between timepoints
%  [FUP-T1] for FND patients and between HC and FND patients at 01_T1
%  (baseline)

for p = 1 % does nothing, just to wrap the code [+/- on the left in MATLAB]
    if show_plots(2)
        %% Differences in # of streamlines between FND [FUP-T1]
        tmp_f1 = figure('Position',[1800 100 900 800]);
        %tmp_f1 = gcf;
        subplot(211);
        imagesc(mean_connectomes_FND{2,1}-mean_connectomes_FND{1,1}); % difference between FUP and 01_T1 for mean # of streamlines
        axis square; colorbar; caxis([-500 500]); xticks(10:10:80); yticks(10:10:80); tmp_f1.CurrentAxes.FontSize = 12;
        title({'Differences in # of streamlines','FND [FUP-T1]'},'FontWeight','normal','Interpreter','none');
        tmp_f1.CurrentAxes.Title.FontSize = 16;
        tmp_f1.CurrentAxes.Title.Position(2) = tmp_f1.CurrentAxes.Title.Position(2)-2;
        xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
        subplot(212);
        plot(jUpperTriMatToVec(triu(mean_connectomes_FND{2,1}-mean_connectomes_FND{1,1},1)));
        ylim([-1150 1150]); xlim([0 3500]);
        xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
        ylabel('Difference in # of streamlines','FontSize',14);
        % compute additionally the inter-hemispheric and intra-hemispheric average differences
        tmp_M_in = mean_connectomes_FND{2,1}-mean_connectomes_FND{1,1};
        tmp_M_in(35,84) = NaN; % we mask the cerebellar difference because it is out of scale with the used atlas
        [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
        title({['Absolute average intra-hemispheric diff. between FND [FUP-T1]: ' num2str(tmp_intra)],...
            ['Absolute average inter-hemispheric diff. between FND [FUP-T1]: ' num2str(tmp_inter)]},...
            'FontWeight','normal','Interpreter','none');
        tmp_f1.CurrentAxes.Title.FontSize = 16;
        tmp_f1.CurrentAxes.Title.Position(2) = tmp_f1.CurrentAxes.Title.Position(2)+100;
        %hold on; tmp_xlim = xlim;
        %plot([tmp_xlim(1) tmp_xlim(2)],repmat(...
        %    mean(jUpperTriMatToVec(triu(tmp_M_in,1)),'omitnan'),1,2),'r-.');
        %hold off; clear tmp_M_in tmp_xlim;
        if ~exist([export_figs_path 'Figure2a.pdf'],'file')
            %export_fig([export_figs_path 'Figure2a.pdf'],tmp_f1);
            % exportgraphics(tmp_f1,[export_figs_path 'Figure2a.pdf'],...
            %     'BackgroundColor','none','ContentType','vector');
        end

        %% Differences of mean lengths between FND [FUP-T1]
        tmp_f2 = figure('Position',[1800 100 900 800]);
        subplot(211);
        imagesc(mean_connectomes_FND{2,2}-mean_connectomes_FND{1,2});
        axis square; colorbar; caxis([-20 20]); xticks(10:10:80); yticks(10:10:80); tmp_f2.CurrentAxes.FontSize = 12;
        title({'Differences in mean length',...
            'FND [FUP-T1]'},'FontWeight','normal','Interpreter','none');
        tmp_f2.CurrentAxes.Title.FontSize = 16;
        tmp_f2.CurrentAxes.Title.Position(2) = tmp_f2.CurrentAxes.Title.Position(2)-2;
        xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
        subplot(212);
        plot(jUpperTriMatToVec(triu(mean_connectomes_FND{2,2}-mean_connectomes_FND{1,2},1)));
        ylim([-30 30]); xlim([0 3500]);
        xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
        ylabel('Difference in mean length','FontSize',14);
        % compute additionally the inter-hemispheric and intra-hemispheric average differences
        tmp_M_in = mean_connectomes_FND{2,2}-mean_connectomes_FND{1,2};
        %tmp_M_in(35,84) = NaN; % we mask the cerebellar difference because it is out of scale with the used atlas
        [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
        title({['Absolute average intra-hemispheric diff. between FND [FUP-T1]: ' num2str(tmp_intra)],...
            ['Absolute average inter-hemispheric diff. between FND [FUP-T1]: ' num2str(tmp_inter)]},...
            'FontWeight','normal','Interpreter','none');
        tmp_f2.CurrentAxes.Title.FontSize = 16;
        tmp_f2.CurrentAxes.Title.Position(2) = tmp_f2.CurrentAxes.Title.Position(2)+2;
        hold on; tmp_xlim = xlim;
        plot([tmp_xlim(1) tmp_xlim(2)],repmat(...
            mean(jUpperTriMatToVec(triu(tmp_M_in,1))),1,2),'r-.');
        hold off; clear tmp_M_in tmp_xlim;
        if ~exist([export_figs_path 'Figure2b.pdf'],'file')
            %export_fig([export_figs_path 'Figure2b.pdf'],tmp_f2);
            % exportgraphics(tmp_f2,[export_figs_path 'Figure2b.pdf'],...
            %     'BackgroundColor','none','ContentType','vector');
        end

        %% Differences of FA between FND [FUP-T1]
        tmp_f3 = figure('Position',[1800 100 900 800]);
        subplot(211);
        imagesc(mean_connectomes_FND{2,3}-mean_connectomes_FND{1,3});
        axis square; colorbar; caxis([-.1 .1]); xticks(10:10:80); yticks(10:10:80); tmp_f3.CurrentAxes.FontSize = 12;
        title({'Differences in mean FA',...
            'FND [FUP-T1]'},'FontWeight','normal','Interpreter','none');
        tmp_f3.CurrentAxes.Title.FontSize = 16;
        tmp_f3.CurrentAxes.Title.Position(2) = tmp_f3.CurrentAxes.Title.Position(2)-2;
        xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
        subplot(212);
        plot(jUpperTriMatToVec(triu(mean_connectomes_FND{2,3}-mean_connectomes_FND{1,3},1)));
        ylim([-.15 .15]); xlim([0 3500]);
        xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
        ylabel('Difference in FA','FontSize',14);
        % compute additionally the inter-hemispheric and intra-hemispheric average differences
        tmp_M_in = mean_connectomes_FND{2,3}-mean_connectomes_FND{1,3};
        %tmp_M_in(35,84) = NaN; % we mask the cerebellar difference because it is out of scale with the used atlas
        [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
        title({['Absolute average intra-hemispheric diff. between FND [FUP-T1]: ' num2str(tmp_intra)],...
            ['Absolute average inter-hemispheric diff. between FND [FUP-T1]: ' num2str(tmp_inter)]},...
            'FontWeight','normal','Interpreter','none');
        tmp_f3.CurrentAxes.Title.FontSize = 16;
        tmp_f3.CurrentAxes.Title.Position(2) = tmp_f3.CurrentAxes.Title.Position(2)+.01;
        hold on; tmp_xlim = xlim;
        plot([tmp_xlim(1) tmp_xlim(2)],repmat(...
            mean(jUpperTriMatToVec(triu(tmp_M_in,1))),1,2),'r-.');
        hold off; clear tmp_M_in tmp_xlim;
        if ~exist([export_figs_path 'Figure2c.pdf'],'file')
            %export_fig([export_figs_path 'Figure2c.pdf'],tmp_f3);
            % exportgraphics(tmp_f3,[export_figs_path 'Figure2c.pdf'],...
            %     'BackgroundColor','none','ContentType','vector');
        end

        %% Differences of # of streamlines between HC-FND [T1]
        tmp_f4 = figure('Position',[1800 100 900 800]);
        subplot(211);
        imagesc(mean_connectomes_HC{1,1}-mean_connectomes_FND{1,1}); % difference between FUP and 01_T1 for mean # of streamlines
        axis square; colorbar; caxis([-500 500]); xticks(10:10:80); yticks(10:10:80); tmp_f4.CurrentAxes.FontSize = 12;
        title({'Differences in # of streamlines','HC-FND [T1]'},'FontWeight','normal','Interpreter','none');
        tmp_f4.CurrentAxes.Title.FontSize = 16;
        tmp_f4.CurrentAxes.Title.Position(2) = tmp_f4.CurrentAxes.Title.Position(2)-2;
        xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
        subplot(212);
        plot(jUpperTriMatToVec(triu(mean_connectomes_HC{1,1}-mean_connectomes_FND{1,1},1)));
        ylim([-1150 1150]); xlim([0 3500]); % ylim kept as for FND [FUP-T1] for visual comparability between plots
        xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
        ylabel('Difference in # of streamlines','FontSize',14);
        % compute additionally the inter-hemispheric and intra-hemispheric average differences
        tmp_M_in = mean_connectomes_HC{1,1}-mean_connectomes_FND{1,1};
        tmp_M_in(35,84) = NaN; % we mask the cerebellar difference because it is out of scale with the used atlas
        [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
        title({['Absolute average intra-hemispheric diff. between HC and FND [T1]: ' num2str(tmp_intra)],...
            ['Absolute average inter-hemispheric diff. between HC and FND [T1]: ' num2str(tmp_inter)]},...
            'FontWeight','normal','Interpreter','none');
        tmp_f4.CurrentAxes.Title.FontSize = 16;
        tmp_f4.CurrentAxes.Title.Position(2) = tmp_f4.CurrentAxes.Title.Position(2)+100;
        if ~exist([export_figs_path 'Figure2d.pdf'],'file')
            %export_fig([export_figs_path 'Figure2d.pdf'],tmp_f4);
            % exportgraphics(tmp_f4,[export_figs_path 'Figure2d.pdf'],...
            %     'BackgroundColor','none','ContentType','vector');
        end

        %% Differences of mean length between HC-FND [T1]
        tmp_f5 = figure('Position',[1800 100 900 800]);
        subplot(211);
        imagesc(mean_connectomes_HC{1,2}-mean_connectomes_FND{1,2});
        axis square; colorbar; caxis([-20 20]); xticks(10:10:80); yticks(10:10:80); tmp_f5.CurrentAxes.FontSize = 12;
        title({'Differences in mean length',...
            'HC-FND [T1]'},'FontWeight','normal','Interpreter','none');
        tmp_f5.CurrentAxes.Title.FontSize = 16;
        tmp_f5.CurrentAxes.Title.Position(2) = tmp_f5.CurrentAxes.Title.Position(2)-2;
        xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
        subplot(212);
        plot(jUpperTriMatToVec(triu(mean_connectomes_HC{1,2}-mean_connectomes_FND{1,2},1)));
        ylim([-30 30]); xlim([0 3500]);
        xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
        ylabel('Difference in mean length','FontSize',14);
        % compute additionally the inter-hemispheric and intra-hemispheric average differences
        tmp_M_in = mean_connectomes_HC{1,2}-mean_connectomes_FND{1,2};
        %tmp_M_in(35,84) = NaN; % we mask the cerebellar difference because it is out of scale with the used atlas
        [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
        title({['Absolute average intra-hemispheric diff. between HC and FND [T1]: ' num2str(tmp_intra)],...
            ['Absolute average inter-hemispheric diff. between HC and FND [T1]: ' num2str(tmp_inter)]},...
            'FontWeight','normal','Interpreter','none');
        tmp_f5.CurrentAxes.Title.FontSize = 16;
        tmp_f5.CurrentAxes.Title.Position(2) = tmp_f5.CurrentAxes.Title.Position(2)+2;
        hold on; tmp_xlim = xlim;
        plot([tmp_xlim(1) tmp_xlim(2)],repmat(...
            mean(jUpperTriMatToVec(triu(tmp_M_in,1))),1,2),'r-.');
        hold off; clear tmp_M_in tmp_xlim;
        if ~exist([export_figs_path 'Figure2e.pdf'],'file')
            %export_fig([export_figs_path 'Figure2e.pdf'],tmp_f5);
            % exportgraphics(tmp_f5,[export_figs_path 'Figure2e.pdf'],...
            %     'BackgroundColor','none','ContentType','vector');
        end

        %% Differences of FA between HC and FND at [T1]
        tmp_f6 = figure('Position',[1800 100 900 800]);
        subplot(211);
        imagesc(mean_connectomes_HC{1,3}-mean_connectomes_FND{1,3});
        axis square; colorbar; caxis([-.1 .1]); xticks(10:10:80); yticks(10:10:80); tmp_f6.CurrentAxes.FontSize = 12;
        title({'Differences in mean FA','HC-FND [T1]'},'FontWeight','normal','Interpreter','none');
        tmp_f6.CurrentAxes.Title.FontSize = 16;
        tmp_f6.CurrentAxes.Title.Position(2) = tmp_f6.CurrentAxes.Title.Position(2)-2;
        xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
        subplot(212);
        plot(jUpperTriMatToVec(triu(mean_connectomes_HC{1,3}-mean_connectomes_FND{1,3},1)));
        ylim([-.15 .15]); xlim([0 3500]);
        xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
        ylabel('Difference in FA','FontSize',14);
        % compute additionally the inter-hemispheric and intra-hemispheric average differences
        tmp_M_in = mean_connectomes_HC{1,3}-mean_connectomes_FND{1,3};
        %tmp_M_in(35,84) = NaN; % we mask the cerebellar difference because it is out of scale with the used atlas
        [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
        title({['Absolute average intra-hemispheric diff. between HC and FND [T1]: ' num2str(tmp_intra)],...
            ['Absolute average inter-hemispheric diff. between HC and FND [T1]: ' num2str(tmp_inter)]},...
            'FontWeight','normal','Interpreter','none');
        tmp_f6.CurrentAxes.Title.FontSize = 16;
        tmp_f6.CurrentAxes.Title.Position(2) = tmp_f6.CurrentAxes.Title.Position(2)+.01;
        hold on; tmp_xlim = xlim;
        plot([tmp_xlim(1) tmp_xlim(2)],repmat(...
            mean(jUpperTriMatToVec(triu(tmp_M_in,1))),1,2),'r-.');
        hold off; clear tmp_M_in tmp_xlim;
        if ~exist([export_figs_path 'Figure2f.pdf'],'file')
            %export_fig([export_figs_path 'Figure2f.pdf'],tmp_f6);
            %exportgraphics(tmp_f6,[export_figs_path 'Figure2f.pdf'],...
            %    'BackgroundColor','none','ContentType','vector');
        end
    end
end



%% Compute WD (weighted-degree graph metric) from FA connectomes
%  to reproduce approach from Ibai Diez et al. (2021), in which
%  first significantly different regions (ROIs) in terms of WD
%  among their 87 ROIs parcellation (ours is 84)
%  were identified and then all links were explored within this
%  subset of significantly altered regions (i.e., proxy for
%  assessing the integrity of white matter originating from the
%  given region(s))

WD_HC = cell(connectomes_HC_ind{1,2},1); % number of HC subjects
WD_FND{1,1} = cell(connectomes_FND_ind{1,2},1);
WD_FND{2,1} = cell(connectomes_FND_ind{2,2},1);
for j = 1:length(WD_HC) % for each HC subject
    WD_HC{j,1} = zeros(size(connectomes_fa{1,1},1),1);
    for n = 1:size(connectomes_fa{1,1},1) % for each node (i.e., region), EXCLUDING self-connections (i=j)
        WD_HC{j,1}(n,1) = sum([...
            connectomes_fa{connectomes_HC_ind{1,1}(j),1}(1:(n-1),n)' ...
            connectomes_fa{connectomes_HC_ind{1,1}(j),1}(n,(n+1):end)]); % because connectomes_* is not symmetrized
    end
end
for j = 1:length(WD_FND{1,1})
    WD_FND{1,1}{j,1} = zeros(size(connectomes_fa{1,1},1),1);
    for n = 1:size(connectomes_fa{1,1},1)
        WD_FND{1,1}{j,1}(n,1) = sum([...
            connectomes_fa{connectomes_FND_ind{1,1}(j),1}(1:(n-1),n)' ...
            connectomes_fa{connectomes_FND_ind{1,1}(j),1}(n,(n+1):end)]);
    end
end
for j = 1:length(WD_FND{2,1})
    WD_FND{2,1}{j,1} = zeros(size(connectomes_fa{1,1},1),1);
    for n = 1:size(connectomes_fa{1,1},1)
        WD_FND{2,1}{j,1}(n,1) = sum([...
            connectomes_fa{connectomes_FND_ind{2,1}(j),2}(1:(n-1),n)' ...
            connectomes_fa{connectomes_FND_ind{2,1}(j),2}(n,(n+1):end)]);
    end
end

% Now we can compare both groups (HC vs FND at [T1])
% in all ROIs [84 for now] accounting for covariates of interest
% (as per the previously loaded covariates Excel file) using a regression
% approach
Y_reg  = ...
    [cell2mat(WD_HC')'; ...
     cell2mat(WD_FND{1,1}')'];

% X may vary according to model used (i.e., adding, or removing covariates)

% X      = ...
%     [%cell2mat(WD_HC')' ... % these are all the WD for HCs (75x84)
%      double(table2array(patients_datatable(connectomes_HC_ind{1,1},4))) ... % age covariate
%      double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_HC_ind{1,1},3)))+1) ... % gender covariate
%      table2array(patients_datatable(connectomes_HC_ind{1,1},7)) ... % BDI
%      table2array(patients_datatable(connectomes_HC_ind{1,1},9)); ... % STAI Y-2 (trait)
%      %zeros(length(WD_HC),1) ... % include a set of 0 as symptom duration vector for HCs
%      %zeros(length(WD_HC),1); ... % include a set of 0 as psychmed vector for HCs
%      %cell2mat(WD_FND{1,1}')' ... % these are all the WD for FNDs (at [T1]) (85x84)
%      double(table2array(patients_datatable(connectomes_FND_ind{1,1},4))) ... % age covariate
%      double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_FND_ind{1,1},3)))+1) ... % gender covariate
%      table2array(patients_datatable(connectomes_FND_ind{1,1},7)) ... % BDI
%      table2array(patients_datatable(connectomes_FND_ind{1,1},9)) ... % STAI Y-2 (trait)
%      %table2array(patients_datatable(connectomes_FND_ind{1,1},6)) ... % symptom duration vector for FND patients
%      %table2array(patients_datatable(connectomes_FND_ind{1,1},5)) % psych med vector for FND patients (0 or 1)
%      ];

% X      = ...
%     [
%      double(table2array(patients_datatable(connectomes_HC_ind{1,1},4))) ... % age covariate
%      double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_HC_ind{1,1},3)))+1) ... % gender covariate
%      zeros(length(WD_HC),1); % include a set of 0 as psychmed vector for HCs
%      double(table2array(patients_datatable(connectomes_FND_ind{1,1},4))) ... % age covariate
%      double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_FND_ind{1,1},3)))+1) ... % gender covariate
%      table2array(patients_datatable(connectomes_FND_ind{1,1},5)) ... % psych med vector for FND patients (0 or 1)
%      ];

% only with antidepressant medication
% X      = ...
%     [
%      double(table2array(patients_datatable(connectomes_HC_ind{1,1},4))) ... % age covariate
%      double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_HC_ind{1,1},3)))+1) ... % gender covariate
%      zeros(length(WD_HC),1); % include a set of 0 as psychmed vector for HCs
%      double(table2array(patients_datatable(connectomes_FND_ind{1,1},4))) ... % age covariate
%      double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_FND_ind{1,1},3)))+1) ... % gender covariate
%      table2array(patients_datatable_meds(connectomes_FND_ind{1,1},22)) ... % antidepressant med vector for FND patients (0 or 1)
%      ];

X      = ...
    [double(table2array(patients_datatable(connectomes_HC_ind{1,1},4))) ... % age covariate
     double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_HC_ind{1,1},3)))); ...
     %table2array(patients_datatable(connectomes_HC_ind{1,1},7)) ... % BDI
     %table2array(patients_datatable(connectomes_HC_ind{1,1},8)) ... % STAI Y-1 (state)
     %table2array(patients_datatable(connectomes_HC_ind{1,1},9)); ... % STAI Y-2 (trait)
     %zeros(length(WD_HC),1) ... % include a set of 0 as symptom duration vector for HCs
     %zeros(length(WD_HC),1); ... % include a set of 0 as psychmed vector for HCs
     double(table2array(patients_datatable(connectomes_FND_ind{1,1},4))) ... % age covariate
     double(cellfun(@(x) strcmp(x,'female'),table2array(patients_datatable(connectomes_FND_ind{1,1},3)))) ...
     %table2array(patients_datatable(connectomes_FND_ind{1,1},7)) ... % BDI
     %table2array(patients_datatable(connectomes_FND_ind{1,1},8)) ... % STAI Y-1 (state)
     %table2array(patients_datatable(connectomes_FND_ind{1,1},9)) ... % STAI Y-2 (trait)
     %table2array(patients_datatable(connectomes_FND_ind{1,1},6)) ... % symptom duration vector for FND patients
     %table2array(patients_datatable(connectomes_FND_ind{1,1},5)) % psych med vector for FND patients (0 or 1)
     ];

X = [zscore(X) ones(size(X,1),1)];

for j = 1:size(Y_reg,2)
    %[~,tmp_residuals] = y_regress_ss(Y_reg(:,j),X);
    [~,~,tmp_residuals] = regress(Y_reg(:,j),X);
    Y_reg(:,j) = mean(Y_reg(:,j)) + tmp_residuals;
end

% Simple comparison with ttest2 to assess whether there
% are differences between certain nodes (between HCs and FNDs)
%   Here first without covariates regression
tmp_HC1 = cell2mat(WD_HC')'; tmp_FND1 = cell2mat(WD_FND{1,1}')';
h_ = zeros(size(tmp_HC1,2),2); p_ = h_; % init (first column without reg; second column with reg of covariates)
stats_ = struct;
tmp_p_val = 0.05; % uncorrected p-value for significance threshold at single comparison
%tmp_p_val_corr = tmp_p_val*(size(tmp_HC1,2)+1)/(2*size(tmp_HC1,2));
%tmp_p_val_corr = 0.05/size(tmp_HC1,2); % e.g. Bonferroni but very strict
tmp_q_val = 0.05; %*tmp_p_val; % adj q-value
for j = 1:size(tmp_HC1,2)
    [h_(j,1),p_(j,1)] = ttest2(tmp_HC1(:,j),tmp_FND1(:,j),...
        'alpha',tmp_p_val,'tail','right','vartype','unequal');
    [h_(j,2),p_(j,2),~,stats_.(['stat' num2str(j)])] = ttest2(Y_reg(1:length(WD_HC),j),Y_reg(length(WD_HC)+1:end,j),...
        'alpha',tmp_p_val,'tail','right','vartype','unequal');
    % p_(j,2) = kruskalwallis(Y_reg(:,j),[ones(75,1);2*ones(85,1)],'off');
    % [p_(j,3),h_(j,3)] = ranksum(Y_reg(1:length(WD_HC),j),Y_reg(length(WD_HC)+1:end,j),...
    %   'alpha',tmp_p_val,'tail','right'); % non-parametric
end

% Adjusting for FDR (Benjamini-Hochberg)
h__ = zeros(size(tmp_HC1,2),2);
%[h__(:,1), crit_p, ajd_ci_cvrg, adj_p] = fdr_bh(p_(:,1),tmp_q_val,'pdep','yes');
h__(:,1) = fdr_bh(p_(:,1),tmp_p_val,'pdep','yes');
h__(:,2) = fdr_bh(p_(:,2),tmp_q_val,'pdep','yes'); % regressed-out with the chosen covariates above


%%% Gather t-stats from significantly differing regions (as in h__(:,2))
% tmp_tvals = zeros(length(find(h__(:,2))),1); cc = 1;
% for j = find(h__(:,2))' % for all sign. regions, store corresp. t-val
%     tmp_tvals(cc) = stats_.(['stat' num2str(j)]).tstat; cc = cc+1;
% end
% clear cc j;

%%% Fill in the MNI .nii map for t-stats at specific region's locations (temporary code for Figure)
% cc = find(h__(:,2)); ccc = 1;
% for j = 1:size(h__(:,2),1) % 84 ROIs
%     if ismember(j,cc)
%         tmp_MNImod(tmp_MNImod==cc(ccc)) = tmp_tvals(ccc);
%         ccc = ccc+1;
%     else
%         tmp_MNImod(tmp_MNImod==j) = 0;
%     end
% end



%% Additional analyses: # of streamlines + mean fiber length (ML)

tmp_nb_streamlines_HC = zeros(length(connectomes{1,1})*(length(connectomes{1,1})-1)/2,length(connectomes_HC_ind{1,1}));
tmp_ml_HC = tmp_nb_streamlines_HC;
for j = 1:length(connectomes_HC_ind{1,1})
    tmp_nb_streamlines_HC(:,j) = jUpperTriMatToVec(connectomes{connectomes_HC_ind{1,1}(j),1},1);
    tmp_ml_HC(:,j) = jUpperTriMatToVec(connectomes_ml{connectomes_HC_ind{1,1}(j),1},1);
end
tmp_nb_streamlines_FND = zeros(length(connectomes{1,1})*(length(connectomes{1,1})-1)/2,length(connectomes_FND_ind{1,1}));
tmp_ml_FND = tmp_nb_streamlines_FND;
for j = 1:length(connectomes_FND_ind{1,1})
    tmp_nb_streamlines_FND(:,j) = jUpperTriMatToVec(connectomes{connectomes_FND_ind{1,1}(j),1},1);
    tmp_ml_FND(:,j) = jUpperTriMatToVec(connectomes_ml{connectomes_FND_ind{1,1}(j),1},1);
end
Y_reg_nbstreamlines = [tmp_nb_streamlines_HC tmp_nb_streamlines_FND]';
Y_reg_ml = [tmp_ml_HC tmp_ml_FND]';
clear tmp_nb_streamlines_* tmp_ml_HC tmp_ml_FND;

% note that X is defined in the cell above
for j = 1:size(Y_reg_nbstreamlines,2)
    %[~,tmp_residuals] = y_regress_ss(Y_reg(:,j),X);
    [~,~,tmp_residuals] = regress(Y_reg_nbstreamlines(:,j),X);
    Y_reg_nbstreamlines(:,j) = mean(Y_reg_nbstreamlines(:,j)) + tmp_residuals;
    [~,~,tmp_residuals] = regress(Y_reg_ml(:,j),X);
    Y_reg_ml(:,j) = mean(Y_reg_ml(:,j)) + tmp_residuals;
end
Y_reg_ml(Y_reg_ml<0) = 0;
Y_reg_nbstreamlines(Y_reg_nbstreamlines<0) = 0;

% Simple comparison with ttest2 to assess whether there
% are differences in # of streamlines or mean fiber length between HCs and FNDs
%   Here first without covariates regression
h_ = zeros(size(Y_reg_ml,2),1); p_ = h_; % init (first column without reg; second column with reg of covariates)
stats_ = struct;
tmp_p_val = 0.05; % uncorrected p-value for significance threshold at single comparison
%tmp_p_val_corr = 0.05/size(tmp_HC1,2); % e.g. Bonferroni but very strict
tmp_q_val = 0.05; %*tmp_p_val; % adj q-value
for j = 1:size(Y_reg_ml,2)
    [h_(j,1),p_(j,1),~,stats_.(['stat' num2str(j)])] = ttest2(...
        Y_reg_nbstreamlines(1:length(connectomes_HC_ind{1,1}),j),...
        Y_reg_nbstreamlines((length(connectomes_HC_ind{1,1})+1):end,j),...
        'alpha',tmp_p_val,'tail','right','vartype','unequal');
end

% Adjusting for FDR (Benjamini-Hochberg)
h__ = zeros(size(Y_reg_ml,2),1);
h__(:,1) = fdr_bh(p_(:,1),tmp_q_val,'pdep','yes');



%% Inspect WD distributions in HCs and FND patients
if show_plots(3)
    % First without the regressed-out covariates above
    if exist('tmp_f7','var') && isvalid(tmp_f7), close(tmp_f7); clear tmp_f7; pause(0.1); end
    tmp_f7 = figure('Position',[60 60 1600 900]);
    subplot(211);
    plot(mean(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),2),...
        'Color',[44 44 255]./256,'LineWidth',1.5); hold on;
    plot(mean(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),2),...
        'Color',[255 44 44]./256,'LineWidth',1.5);
    tmp_l7 = legend('HC','FND','FontSize',20,'Location','southeast','FontName','Basis Grotesque Pro');
    %
    plot(mean(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),2)+...
        std(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),[],2),...
        'Color',[255 123 123]./256,'LineWidth',0.7);
    plot(mean(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),2)-...
        std(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),[],2),...
        'Color',[255 123 123]./256,'LineWidth',0.7);
    %
    plot(mean(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),2)+...
        std(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),[],2),...
        'Color',[78 145 253]./256,'LineWidth',0.7);
    plot(mean(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),2)-...
        std(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),[],2),...
        'Color',[78 145 253]./256,'LineWidth',0.7);
    for j = 1:length(h__(:,1))
        if h__(j,1), text(j-.35,55,'*','Color','black','FontSize',20,'FontName','Basis Grotesque Pro'); % [j-.35 for graphical purpose only]
        end
    end
    title(['Average baseline WDs across ' num2str(length(WD_HC{1,1})) ' ROIs in ' num2str(length(WD_FND{1,1})) ' FND vs ' ...
        num2str(length(WD_HC)) ' HCs with q = .05 (FDR)'],'FontWeight','normal','Interpreter','none','FontName','Basis Grotesque Pro');
    axis([0 85 0 60]); xticks(4:4:84); yticks(0:5:60); tmp_f7.CurrentAxes.FontSize = 12;
    tmp_f7.CurrentAxes.YGrid = 'on';
    xlabel('ROIs','FontSize',18,'FontName','Basis Grotesque Pro');
    ylabel('Weighted degree','Interpreter','none','Rotation',90,'FontSize',18,'FontName','Basis Grotesque Pro');
    tmp_f7.CurrentAxes.YAxis.Label.Position(1) = tmp_f7.CurrentAxes.YAxis.Label.Position(1)-1;
    tmp_f7.CurrentAxes.Title.FontSize = 20;
    tmp_f7.CurrentAxes.Title.Position(2) = tmp_f7.CurrentAxes.Title.Position(2)+2;
    tmp_l7.String = tmp_l7.String(1:2); % remove other lines from legend
    subplot(212);
    %histogram(mean(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),2),'BinMethod','scott',...
    %    'Normalization','pdf','FaceColor',[44 44 255]./256); hold on;
    tmp_h1 = histfit(mean(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),2),10); hold on;
    %histogram(mean(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),2),...
    %    'BinMethod','scott','Normalization','pdf','FaceColor',[255 44 44]./256);
    tmp_h2 = histfit(mean(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),2),10);
    tmp_h1(1).FaceAlpha = 0.75; tmp_h2(1).FaceAlpha = 0.75;
    tmp_h1(1).FaceColor = [78 145 253]./256; tmp_h2(1).FaceColor = [255 123 123]./256;
    tmp_h1(2).Color = [44 44 255 256*.75]./256; tmp_h2(2).Color = [255 44 44 256*.75]./256;
    axis([10 60 0 24]); xticks(15:5:60); yticks(0:4:24); tmp_f7.CurrentAxes.FontSize = 12;
    tmp_f7.CurrentAxes.YGrid = 'on';
    xlabel('Weighted degree','FontSize',18,'FontName','Basis Grotesque Pro');
    ylabel('Count','Interpreter','none','Rotation',90,'FontSize',18,'FontName','Basis Grotesque Pro');
    tmp_f7.CurrentAxes.YAxis.Label.Position(1) = tmp_f7.CurrentAxes.YAxis.Label.Position(1)-.5;
    tmp_l7 = legend('HC','','FND','FontSize',20,'Location','northeast','FontName','Basis Grotesque Pro');
    title(['Distribution of average baseline WDs across ' num2str(length(WD_HC{1,1})) ' ROIs in ' ...
        num2str(length(WD_FND{1,1})) ' FND vs ' num2str(length(WD_HC)) ' HCs'],...
        'FontWeight','normal','Interpreter','none','FontName','Basis Grotesque Pro');
    tmp_f7.CurrentAxes.Title.FontSize = 20;
    tmp_f7.CurrentAxes.Title.Position(2) = tmp_f7.CurrentAxes.Title.Position(2)+1;
    if ~exist([export_figs_path 'Figure3a.pdf'],'file')
        %export_fig([export_figs_path 'Figure3a.pdf'],tmp_f7);
        %exportgraphics(tmp_f7,[export_figs_path 'Figure3a.pdf'],...
        %    'BackgroundColor','none','ContentType','vector');
    end

    %% With regressed-out covariates (q = .05; see above)
    if exist('tmp_f8','var') && isvalid(tmp_f8), close(tmp_f8); clear tmp_f8; pause(0.1); end
    tmp_f8 = figure('Position',[60 60 1600 900]);
    subplot(211);
    plot(mean(Y_reg(1:length(WD_HC),:),1),...
        'Color',[44 44 255]./256,'LineWidth',1.5); hold on;
    plot(mean(Y_reg(length(WD_HC)+1:end,:),1),...
        'Color',[255 44 44]./256,'LineWidth',1.5);
    tmp_l7 = legend('HC','FND','FontSize',20,'Location','southeast','FontName','Basis Grotesque Pro'); % [j-.35 for graphical purpose only]
    %
    plot(mean(Y_reg(length(WD_HC)+1:end,:),1)+...
        std(Y_reg(length(WD_HC)+1:end,:),[],1),...
        'Color',[255 123 123]./256,'LineWidth',0.7);
    plot(mean(Y_reg(length(WD_HC)+1:end,:),1)-...
        std(Y_reg(length(WD_HC)+1:end,:),[],1),...
        'Color',[255 123 123]./256,'LineWidth',0.7);
    %
    plot(mean(Y_reg(1:length(WD_HC),:),1)+...
        std(Y_reg(1:length(WD_HC),:),[],1),...
        'Color',[78 145 253]./256,'LineWidth',0.7);
    plot(mean(Y_reg(1:length(WD_HC),:),1)-...
        std(Y_reg(1:length(WD_HC),:),[],1),...
        'Color',[78 145 253]./256,'LineWidth',0.7);
    for j = 1:length(h__(:,2))
        if h__(j,2), text(j-.35,55,'*','Color','black','FontSize',20,'FontName','Basis Grotesque Pro');
        end
    end
    title(['Average baseline WDs across ' num2str(length(WD_HC{1,1})) ' ROIs in ' num2str(length(WD_FND{1,1})) ' FND vs ' ...
        num2str(length(WD_HC)) ' HCs with q = 0.05 (FDR), ' ...
        'age and sex regressed out'],'FontWeight','normal','Interpreter','none','FontName','Basis Grotesque Pro');
    axis([0 85 0 60]); xticks(4:4:84); yticks(0:5:60); tmp_f8.CurrentAxes.FontSize = 12;
    tmp_f8.CurrentAxes.YGrid = 'on';
    xlabel('ROIs','FontSize',18);
    ylabel('Weighted degree','Interpreter','none','Rotation',90,'FontSize',18,'FontName','Basis Grotesque Pro');
    tmp_f8.CurrentAxes.YAxis.Label.Position(1) = tmp_f8.CurrentAxes.YAxis.Label.Position(1)-1;
    tmp_f8.CurrentAxes.Title.FontSize = 20;
    tmp_f8.CurrentAxes.Title.Position(2) = tmp_f8.CurrentAxes.Title.Position(2)+2;
    tmp_l7.String = tmp_l7.String(1:2); % remove other lines from legend
    subplot(212);
    %histogram(mean(reshape(cell2mat(WD_HC),length(WD_HC{1,1}),length(WD_HC)),2),'BinMethod','scott',...
    %    'Normalization','pdf','FaceColor',[44 44 255]./256); hold on;
    tmp_h1 = histfit(mean(Y_reg(1:length(WD_HC),:),1),10); hold on;
    %histogram(mean(reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})),2),...
    %    'BinMethod','scott','Normalization','pdf','FaceColor',[255 44 44]./256);
    tmp_h2 = histfit(mean(Y_reg(length(WD_HC)+1:end,:),1),10);
    tmp_h1(1).FaceAlpha = 0.75; tmp_h2(1).FaceAlpha = 0.75;
    tmp_h1(1).FaceColor = [78 145 253]./256; tmp_h2(1).FaceColor = [255 123 123]./256;
    tmp_h1(2).Color = [44 44 255 256*.75]./256; tmp_h2(2).Color = [255 44 44 256*.75]./256;
    axis([10 60 0 24]); xticks(15:5:60); yticks(0:4:24); tmp_f8.CurrentAxes.FontSize = 12;
    tmp_f8.CurrentAxes.YGrid = 'on';
    xlabel('Weighted degree','FontSize',18,'FontName','Basis Grotesque Pro');
    ylabel('Count','Interpreter','none','Rotation',90,'FontSize',18,'FontName','Basis Grotesque Pro');
    tmp_f8.CurrentAxes.YAxis.Label.Position(1) = tmp_f8.CurrentAxes.YAxis.Label.Position(1)-.5;
    tmp_l7 = legend('HC','','FND','FontSize',20,'Location','northeast','FontName','Basis Grotesque Pro');
    title(['Distribution of average baseline WDs across ' num2str(length(WD_HC{1,1})) ' ROIs in ' ...
        num2str(length(WD_FND{1,1})) ' FND vs ' num2str(length(WD_HC)) ' HCs, age and sex regressed out'],...
        'FontWeight','normal','Interpreter','none','FontName','Basis Grotesque Pro');
    tmp_f8.CurrentAxes.Title.FontSize = 20;
    tmp_f8.CurrentAxes.Title.Position(2) = tmp_f8.CurrentAxes.Title.Position(2)+1;
    if ~exist([export_figs_path 'Figure3b.pdf'],'file')
        %exportgraphics(tmp_f8,[export_figs_path 'Figure3b.pdf'],...
        %    'BackgroundColor','none','ContentType','vector');
    end
%     if ~exist([export_figs_path 'Figure3c.pdf'],'file')
%         exportgraphics(tmp_f8,[export_figs_path 'Figure3c.pdf'],...
%             'BackgroundColor','none','ContentType','vector');
%     end
end



%% Inspect the links that are significantly different for
%  the ROIs found above using the WD (FND<HC) approach

links_HC_fa  = zeros(84,connectomes_HC_ind{1,2},length(find(h__(:,2)))); % 84 links x number of HCs [at T1] x number of sign. ROIs
links_FND_fa = zeros(84,connectomes_FND_ind{1,2},length(find(h__(:,2)))); % same for FND patients [at T1]
sign_ROIs = find(h__(:,2));
for j = 1:length(sign_ROIs) % for every significant ROI (from regressed-out approach above)
    for p = 1:connectomes_HC_ind{1,2}
        links_HC_fa(1:sign_ROIs(j),p,j) = connectomes_fa{connectomes_HC_ind{1,1}(p),1}(1:sign_ROIs(j),sign_ROIs(j));
        links_HC_fa((sign_ROIs(j)+1):end,p,j) = ...
            connectomes_fa{connectomes_HC_ind{1,1}(p),1}(sign_ROIs(j),(sign_ROIs(j)+1):end);
    end
    for p = 1:connectomes_FND_ind{1,2}
        links_FND_fa(1:sign_ROIs(j),p,j) = connectomes_fa{connectomes_FND_ind{1,1}(p),1}(1:sign_ROIs(j),sign_ROIs(j));
        links_FND_fa((sign_ROIs(j)+1):end,p,j) = ...
            connectomes_fa{connectomes_FND_ind{1,1}(p),1}(sign_ROIs(j),(sign_ROIs(j)+1):end);
    end
end

% Take average across HCs and FND patients to get the average link values
% (84 regions) x number of significant ROIs (i.e., length(sign_ROIs))
avglinks_HC_fa  = squeeze(mean(links_HC_fa,2)); % mean across number of HCs [at T1]
avglinks_FND_fa = squeeze(mean(links_FND_fa,2)); % same; both are now 84xlength(sign_ROIs)

links_HC_fa  = permute(links_HC_fa,[2 1 3]);
links_FND_fa = permute(links_FND_fa,[2 1 3]);
% compute Y_reg2 in case
Y_reg2 = zeros(size(links_HC_fa,1)+size(links_FND_fa,1),size(links_HC_fa,2),length(sign_ROIs));
for j = 1:length(sign_ROIs)
    Y_reg2(:,:,j) = [links_HC_fa(:,:,j);links_FND_fa(:,:,j)];
end
for j = 1:length(sign_ROIs)
    for l = 1:size(Y_reg2,2)
        %[~,tmp_residuals] = y_regress_ss(Y_reg(:,j),X);
        [~,~,tmp_residuals] = regress(Y_reg2(:,l,j),X);
        Y_reg2(:,l,j) = mean(Y_reg2(:,l,j)) + tmp_residuals;
    end
end

% [[
% Approach with whole brain links (irrespective of sign. ROIs from above)
links_HC_fa_wholebrain  = zeros(connectomes_HC_ind{1,2},size(connectomes_fa{1,1},1)*...
    (size(connectomes_fa{1,1},1)-1)/2); % store all FA values in a matrix, without diagonal
links_FND_fa_wholebrain = zeros(connectomes_FND_ind{1,2},size(links_HC_fa_wholebrain,2));

for p = 1:connectomes_HC_ind{1,2}, links_HC_fa_wholebrain(p,:) = jUpperTriMatToVec(connectomes_fa{connectomes_HC_ind{1,1}(p),1},1); end
for p = 1:connectomes_FND_ind{1,2}, links_FND_fa_wholebrain(p,:) = jUpperTriMatToVec(connectomes_fa{connectomes_FND_ind{1,1}(p),1},1); end

Y_reg3 = [links_HC_fa_wholebrain;links_FND_fa_wholebrain];
for l = 1:size(Y_reg3,2)
    [~,~,tmp_residuals] = regress(Y_reg3(:,l),X);
    Y_reg3(:,l) = mean(Y_reg3(:,l)) + tmp_residuals;
end

p4_  = zeros(1,length(Y_reg3));
for j = 1:length(Y_reg3)
    [~,p4_(:,j)] = ttest2(...
        Y_reg3(1:connectomes_HC_ind{1,2},j),...
        Y_reg3((connectomes_HC_ind{1,2}+1):end,j),...
        'tail','right','vartype','unequal');
end
h4__ = fdr_bh(p4_,0.05,'pdep','yes');
h4__ = jVecToUpperTriMat(h4__,size(connectomes_fa{1,1},1));
% ]]

%%% Workaround for ROIs, without testing the same pvals twice (in the fdr approach):
% First create mask in matrix (logical)
h_mask_tmp = false(size(connectomes_fa{1,1},1));
for j = 1:length(sign_ROIs)
    h_mask_tmp(1:(sign_ROIs(j)-1),sign_ROIs(j)) = true; % note the -1 for the removal of the diagonal elements (self ROI-to-ROI FA values not assessed here for links' connectivity)
    h_mask_tmp(sign_ROIs(j),(sign_ROIs(j)+1):end) = true;
end
%figure; imagesc(h_mask_tmp); axis('square'); % visual assessment of mask

tmp_p5__inds = find(jUpperTriMatToVec(h_mask_tmp));
p5_  = zeros(1,length(find(h_mask_tmp))); % p5_ will contain the pvalues of individual t-tests but limited to ROI-to-ROIs of interest (without repetitions)
for j = 1:length(find(h_mask_tmp))
    [~,p5_(1,j)] = ttest2(...
        Y_reg3(1:connectomes_HC_ind{1,2},tmp_p5__inds(j)),...
        Y_reg3((connectomes_HC_ind{1,2}+1):end,tmp_p5__inds(j)),...
        'tail','right','vartype','unequal');
end

h5__ = fdr_bh(p5_,0.05,'pdep','yes'); % assess FDR with q = 0.05
h5__fdr = zeros(1,length(Y_reg3));

tmp_h5_vec = zeros(1,length(Y_reg3));
tmp_h5_vec(tmp_p5__inds) = 1;

ccc = 1;
for j = 1:length(Y_reg3)
    if tmp_h5_vec(1,j)
        if h5__(1,ccc), h5__fdr(1,j) = 1; end
        ccc = ccc+1;
    end
end
clear ccc tmp_h5_vec;

h5__fdr = jVecToUpperTriMat(h5__fdr,size(connectomes_fa{1,1},1));

% for k = 1:size(connectomes_fa{1,1},1) % re-sort back by column first (according to jUpperTriMatToVec)
%     for j = 1:size(connectomes_fa{1,1},1) 
%         if h_mask_tmp(j,k) && 
%             h5__fdr(j,k) = 1;
%         end
%     end
% end


% Simple comparison with ttest2 to assess which links differ
% (between HCs and FNDs) within the selected ROIs [, directly
% done with the same regression of covariates as for the ROIs above]
h1_ = zeros(size(Y_reg2,2),length(sign_ROIs)); p1_ = h1_;
%tmp_p_val = 0.05; % uncorrected p-value for significance threshold at single comparison
%tmp_q_val = 0.05; %*tmp_p_val; % adj q-value
for j = 1:length(sign_ROIs)
    for l = 1:size(Y_reg2,2)
        [h1_(l,j),p1_(l,j)] = ttest2(Y_reg2(1:size(links_HC_fa,1),l,j),Y_reg2(size(links_HC_fa,1)+1:end,l,j),...
            'alpha',tmp_p_val,'tail','right','vartype','unequal');
    end
end

% Adjusting for FDR (Benjamini-Hochberg)
h2__ = zeros(size(Y_reg2,2),length(sign_ROIs));
for j = 1:length(sign_ROIs)
    h2__(:,j) = fdr_bh(p1_(:,j),tmp_q_val,'pdep','yes'); % regressed-out with the chosen covariates above
end



%%
%%% Plot ROIs on brain [simple line-art SVG out] using
%   (an adapted code for) plotBrain.m — Original git page:
%   https://github.com/dutchconnectomelab/Simple-Brain-Plot

DesikanKilliany_atlas_regionlist = load('DesikanKilliany_atlas_regionlist.mat');
aparc_aseg_to_DKilliany_reorder_vec = DesikanKilliany_atlas_regionlist.aparc_aseg_to_DKilliany_reorder_vec;
DesikanKilliany_atlas = DesikanKilliany_atlas_regionlist.DesikanKilliany_atlas;
DesikanKilliany_atlas_regionlist = DesikanKilliany_atlas_regionlist.DesikanKilliany_atlas_regionlist;
DesikanKilliany_atlas_sorted = DesikanKilliany_atlas([1:2:end 2:2:end],:);

% % remove Left-Cerebellum-Cortex and Right-Cerebellum-Cortex from the ROIs
% % (not plotted with the Simple-Brain-Plot)
% sign_ROIs2 = sign_ROIs;
% %if ~isempty(find(ismember(sign_ROIs2,35),1)), sign_ROIs2 = setxor(sign_ROIs2,35); end
% %if ~isempty(find(ismember(sign_ROIs2,84),1)), sign_ROIs2 = setxor(sign_ROIs2,84); end
% 
% find corresponding indices from DesikanKilliany_atlas
% tmp_indz = zeros(length(sign_ROIs2),1);
% for j = 1:length(sign_ROIs2)
%     tmp_ind = find(cellfun(@(jj) str2double(jj),DesikanKilliany_atlas(:,1))==sign_ROIs2(j));
%     if ~isempty(tmp_ind), tmp_indz(j,1) = tmp_ind; else, tmp_indz(j,1) = NaN;
%     end
% end
tmp_indz = zeros(length(sign_ROIs),1);
for j = 1:length(sign_ROIs)
    tmp_ind = find(cellfun(@(jj) str2double(jj),DesikanKilliany_atlas(:,1))==sign_ROIs(j));
    if ~isempty(tmp_ind), tmp_indz(j,1) = tmp_ind; else, tmp_indz(j,1) = NaN;
    end
end
% %DesikanKilliany_atlas(tmp_indz,3); % gives the names of the ROIs

% colormap (just red interp.)
%tmp_cm = [1.0000    1.0000    1.0000; ...
%          0.9102    0.2422    0.2266];
%tmp_cm = interp1(tmp_cm, 1:0.002:size(tmp_cm,1));

% show in the WebViewer of MATLAB
%tmp_indz2 = zeros(length(regionDescriptions.aparc_aseg),1); tmp_indz2(tmp_indz) = 1;
%plotBrain(regionDescriptions.aparc_aseg, tmp_indz2, tmp_cm, 'atlas', 'aparc_aseg', 'scaling', 0.1);



%% Plot all links (metric: FA) with ROIs' names between both groups
% (HCs and FND [at T1])

if exist('tmp_f9','var') && ishandle(tmp_f9), close(tmp_f9); end
tmp_f9 = figure('Position',[30 5 1680 980]);
for j = 1:length(sign_ROIs)
    subplot(7,3,j);
    plot(avglinks_HC_fa(:,j),...
        'Color',[44 44 255]./256,'LineWidth',1.5); hold on;
    plot(avglinks_FND_fa(:,j),...
        'Color',[255 44 44]./256,'LineWidth',1.5);
    tmp_l9 = legend('HC','FND','FontSize',10,'Location','southwest','FontName','Basis Grotesque Pro');
    %
    plot(avglinks_HC_fa(:,j)+...
        std(Y_reg2(1:size(links_HC_fa,1),:,j))',...
        'Color',[255 123 123]./256,'LineWidth',0.7);
    plot(avglinks_HC_fa(:,j)-...
        std(Y_reg2(1:size(links_HC_fa,1),:,j))',...
        'Color',[255 123 123]./256,'LineWidth',0.7);
    %
    plot(avglinks_FND_fa(:,j)+...
        std(Y_reg2((size(links_HC_fa,1)+1):end,:,j))',...
        'Color',[78 145 253]./256,'LineWidth',0.7);
    plot(avglinks_FND_fa(:,j)-...
        std(Y_reg2((size(links_HC_fa,1)+1):end,:,j))',...
        'Color',[78 145 253]./256,'LineWidth',0.7);
    for l = 1:length(h5__fdr(:,sign_ROIs(j))) % here replace by h2__ if you do the fdr_bh test separately (not as all-in-one in matrix)
        if h5__fdr(l,sign_ROIs(j)), text(l-.1,.8,'*','Color','black',...
                'FontSize',16,'FontName','Basis Grotesque Pro'); % [l-.1 only for graphical purposes]  
        end
    end
    tmp_l9.String = tmp_l9.String(1:2); % remove other lines from legend
    axis([0 size(Y_reg2,2)+1 -.2 .9]);
    title(DesikanKilliany_atlas{tmp_indz(j),3},'FontWeight','bold','FontName','Basis Grotesque Pro');
    tmp_f9.CurrentAxes.Title.Position(2) = tmp_f9.CurrentAxes.Title.Position(2) + .05;
end



%% Connectogram of altered links (FA) between HCs and FND [at T1]
%  Makes use of the biChordChart MATLAB addon available on
%  https://ch.mathworks.com/matlabcentral/fileexchange/121043-digraph-chord-chart

%if exist('tmp_f10','var') && ishandle(tmp_f10), close(tmp_f10); end
tmp_f10 = figure('Position',[325 6 1000 1000]);
if ~exist('DesikanKilliany_atlas_connectogram_labels','var')
    load('DesikanKilliany_atlas_connectogram_labels.mat');
end

% significant links into data matrix (dataMat)
dataMat = zeros(size(avglinks_FND_fa,1),size(avglinks_FND_fa,1));
for j = 1:length(sign_ROIs)
%     dataMat(cell2mat(DesikanKilliany_atlas_connectogram_labels(:,2))==...
%       sign_ROIs(j),:) = h5__fdr(:,sign_ROIs(j)); % h3__(:,j);
    tmp_indzz = find(h5__fdr(:,sign_ROIs(j)));
    for k = tmp_indzz'
        dataMat(cell2mat(DesikanKilliany_atlas_connectogram_labels(:,2))==...
            sign_ROIs(j),cell2mat(DesikanKilliany_atlas_connectogram_labels(:,2))==k) = 1;
    end
end
clear tmp_indzz k j;

% appropriate labels (re-ordered manually to match Diez et al., 2021) (NameList)
%NameList = transpose(DesikanKilliany_atlas_sorted(:,3));
NameList = DesikanKilliany_atlas_connectogram_labels(:,1);

BCC = biChordChart(dataMat, 'Label', NameList, ...
    'Arrow', 'off', 'CData', repmat([.2 .2 .2],size(dataMat,1),1)); % hsv(84) / bone(size(dataMat,1))
BCC = BCC.draw(); 
BCC.labelRotate('on');
BCC.setLabelRadius(1.2);
%BCC.tickState('off');
%BCC.tickLabelState('off');
%BCC.setTickFont('FontName','Cambria','Color',[0,0,.6])
%BCC.setFont('FontName','Cambria','FontSize',17,'Color',[.2,.2,.2])
BCC.setFont('FontName', 'Basis Grotesque Arabic Pro', ...
    'FontSize', 16, 'Color', [.2 .2 .2]);

% customized coloring according to Diez et al., 2021
tmp_alllabels = findobj(tmp_f10.CurrentAxes,'Tag','BiChordLabel');
tmp_alllabels = flipud(tmp_alllabels); % [because those were retrieved in reverse order]
%tmp_colorpalette_ChordN = interp1([...
%    [255 186 186]./256; [0 0 167]./256], 1:.0625:(2-.0625));
%tmp_colorpalette_ChordN = repmat(['r';'g';'b';'c';'m';'y';'k';'w'],2,1);
% tmp_colorpalette_ChordN = ...
%     ['#fabebe'; '#3cb44b'; '#ffe119'; ...
%      '#4363d8'; '#f58231'; '#911eb4'; ...
%      '#46f0f0'; '#f032e6'; '#bcf60c'; ...
%      '#e6194b'; '#008080'; '#e6beff'; ...
%      '#9a6324'; '#42d4f4'; '#800000'; ...
%      '#808000'];% '#808000'; '#ffd8b1'; ...
%      %'#000075'; '#808080'; '#ffffff'; '#000000']; % use this if h2__ is used as a correction (fdr_bh)
tmp_colorpalette_ChordN = ...
    ['#fabebe'; '#3cb44b'; '#ffe119'; ... % InfPar InfTemp LatOcc (left)
     '#4363d8'; '#f58231'; '#911eb4'; ...
     '#46f0f0'; '#f032e6'; '#bcf60c'; ...
     '#e6194b'; '#008080'; '#e6beff'; ...
     '#9a6324'; '#42d4f4'; '#800000'; ...
     '#808000'];% '#808000'; '#ffd8b1'; ...
     %'#000075'; '#808080'; '#ffffff'; '#000000'];

tmp_cc = 1;
for j = 1:length(sign_ROIs)
    tmp_ind10 = find(cell2mat(...
        DesikanKilliany_atlas_connectogram_labels(:,2))==sign_ROIs(j));
    tmp_alllabels(tmp_ind10).FontWeight = 'Bold';
    tmp_alllabels(tmp_ind10).FontSize = 18;
    if sum(dataMat(tmp_ind10,:))>0
        tmp_alllabels(tmp_ind10).Color = tmp_colorpalette_ChordN(tmp_cc,:);
        BCC.setSquareN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:));
        %BCC.setSquareN(tmp_ind10,'EdgeColor',[.8 .2 .2],'LineWidth',2.5);
        BCC.setChordN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:)); % sanity check
        tmp_cc = tmp_cc+1;
    else
        tmp_alllabels(tmp_ind10).Color = '#a9a9a9';
        BCC.setSquareN(tmp_ind10,'FaceColor','#a9a9a9');
    end
end
clear tmp_cc tmp_ind10; 

%BCC.setSquareN(2,'FaceColor','black')
%BCC.setChordN(1,'FaceColor','red')
% [m,n]=find(dataMat==max(max(dataMat)));
% for i=1:length(m)
%     BCC.setChordMN(m(i),n(i),'EdgeColor',[.8,0,0],'LineWidth',2)
% end

%exportgraphics(tmp_f10,[pwd filesep 'Figures' filesep 'Figure10.pdf'],...
%            'BackgroundColor','none','ContentType','vector');
%export_fig([pwd filesep 'Figures' filesep 'Figure10a.pdf'],tmp_f10);



%% Add clinical data (∆S-FMDRS and CGI-1 scores) to the patients_datatable
%  Need the updated clinical files called
%  231201 BioGen_FUP_clean SW.xlsx
%  231212 FinalResults_BioGen.xlsx

if exist('patients_datatable','var')
    % [S-FMDRS [T1, FUP] ∆S-FMDRS CGI-1 [T1, FUP] ∆CGI-1 FW (yes:true/no:false) Improved (yes:true/no:false)]
    tmp_vartypes = {'double','double','double','double','double','double','double','double','double'};
    tmp_additional_tabledata = table('Size',[size(patients_datatable,1) 9],'VariableTypes',tmp_vartypes,...
        'VariableNames',{'p_code','S-FMDRS [T1]','S-FMDRS [FUP]','∆S-FMDRS','CGI-1 [T1]','CGI-1 [FUP]',...
        '∆CGI-1','FW','Improved'}); %clear tmp_vartypes; % [before was 'logical' but now all 'double' ...]
    tmp_additional_tabledata.p_code = patients_datatable.p_code;
    tmp_additional_tabledata{:,2:9} = nan; % for all 'double' columns
    %tmp_additional_tabledata{:,8:9} = false;
    tmp_newdata1  = readtable([DWIcode_path '231201 BioGen_FUP_clean SW.xlsx']);
    tmp_FWnFW     = readtable([DWIcode_path '231006 Covariates Update SW' filesep ...
        '231006 FW_Groups SW.xlsx']);
    for j = 1:size(tmp_FWnFW,1)
        tmp_FWnFW{j,1} = {str2double(tmp_FWnFW{j,1}{1}(2:end))}; % [convert to double for the P codes...]
    end
    tmp_improved1 = readtable([DWIcode_path '231006 Covariates Update SW' filesep ...
        '231006 ClinicalOutcomeGroups SW.xlsx']);
    for j = 1:size(tmp_newdata1,1)
        tmp_j = find(tmp_additional_tabledata{:,1}==tmp_newdata1.p_code(j,1),1);
        if ~isempty(tmp_j)
            tmp_additional_tabledata{tmp_j,2:7} = [tmp_newdata1.sfmdrs_t1(j,1) tmp_newdata1.sfmdrs_fup(j,1) ...
                tmp_newdata1.sfmdrs_fup(j,1)-tmp_newdata1.sfmdrs_t1(j,1) tmp_newdata1.cgi_1_t1(j,1) ...
                tmp_newdata1.cgi_1_fup(j,1) tmp_newdata1.cgi_1_fup(j,1)-tmp_newdata1.cgi_1_t1(j,1)];
        end
    end
    for j = 1:size(tmp_FWnFW,1)
        tmp_j = find(tmp_additional_tabledata{:,1}==tmp_FWnFW{j,1}{1});
        if ~isempty(tmp_j) && strcmp(tmp_FWnFW{j,2}{1},'FW')
            tmp_additional_tabledata{tmp_j,8} = true;
        elseif ~isempty(tmp_j) && strcmp(tmp_FWnFW{j,2}{1},'no-FW')
            tmp_additional_tabledata{tmp_j,8} = false;
        end
    end
    for j = 1:size(tmp_improved1,1)
        tmp_j = find(tmp_additional_tabledata{:,1}==tmp_improved1.p_code(j,1));
        if ~isempty(tmp_j) && strcmp(tmp_improved1.group(j,1),'improved')
            tmp_additional_tabledata{tmp_j,9} = true;
        elseif ~isempty(tmp_j) && strcmp(tmp_improved1.group(j,1),'worse')
            tmp_additional_tabledata{tmp_j,9} = false;
        end
    end

    % Merge into existing patients_datatable
    patients_datatable{:,10:17} = tmp_additional_tabledata{:,2:end};
    patients_datatable{:,18:21} = NaN; % add SF36-PhysHealth and SF36-MentalHealth etc.
    patients_datatable.Properties.VariableNames(10:21) = {'S-FMDRS [T1]','S-FMDRS [FUP]',...
        '∆S-FMDRS','CGI-1 [T1]','CGI-1 [FUP]','∆CGI-1','FW','Improved','SF36-PhysHealth',...
        'SF36-MentalHealth','SF36-GenHealth','SF36-PhysFunc'};

    tmp_newdata2  = readtable([DWIcode_path '231212 FinalResults_BioGen.xlsx']);
    tmp_newdata2.sfmdrs(88:end,1) = NaN; % keep NaN for HCs here not to confuse
    for j = 1:size(tmp_newdata2,1)
        tmp_j = find(patients_datatable.p_code==tmp_newdata2.p_code(j,1),1);
        if ~isempty(tmp_j)
            patients_datatable.("S-FMDRS [T1]")(tmp_j,1) = tmp_newdata2.sfmdrs(j,1);
            patients_datatable.("CGI-1 [T1]")(tmp_j,1) = tmp_newdata2.cgi_1(j,1);
            patients_datatable.("SF36-PhysHealth")(tmp_j,1) = tmp_newdata2.sf36_phyhealth(j,1);
            patients_datatable.("SF36-MentalHealth")(tmp_j,1) = tmp_newdata2.sf36_mentalhealth(j,1);
            patients_datatable.("SF36-GenHealth")(tmp_j,1) = tmp_newdata2.sf36_genheal(j,1);
            patients_datatable.("SF36-PhysFunc")(tmp_j,1) = tmp_newdata2.sf36_phyfun(j,1);
        end
    end
    clear j tmp_j tmp_vartypes tmp_FWnFW tmp_newdata1 tmp_newdata2 tmp_improved1 tmp_additional_tabledata;
end



%% Within-group analysis of FND patients for FW vs non-FW [at T1 & FUP when available]
%  for all metrics (FA, fiber density, and mean length)

patients_inds_FW = find(patients_datatable.FW==1); patients_inds_nFW = find(patients_datatable.FW==0);
patients_inds_impr = find(patients_datatable.Improved==1); patients_inds_nimpr = find(patients_datatable.Improved==0);

% manual adjustment for P48: no DTI data is available despite the available
% clinical classification -> remove its corresponding ind(ex) from list
if any(ismember(patients_inds_nFW,5)), patients_inds_nFW = setxor(patients_inds_nFW,5); end

mean_connectomes_FND_FW  = repmat({zeros(size(connectomes{1,1}))},2,3); % to store mean connectome, _ml, and _fa
mean_connectomes_FND_nFW = mean_connectomes_FND_FW; % initialize the same
mean_connectomes_FND_impr = mean_connectomes_FND_FW; % same
mean_connectomes_FND_nimpr= mean_connectomes_FND_FW; % same

for j = 1:length(patients_inds_FW)
    mean_connectomes_FND_FW{1,1} = mean_connectomes_FND_FW{1,1} + connectomes{patients_inds_FW(j),1};
    mean_connectomes_FND_FW{1,2} = mean_connectomes_FND_FW{1,2} + connectomes_ml{patients_inds_FW(j),1};
    mean_connectomes_FND_FW{1,3} = mean_connectomes_FND_FW{1,3} + connectomes_fa{patients_inds_FW(j),1};
    if ismember(patients_inds_FW(j),union(patients_inds_impr,patients_inds_nimpr)) % only if FUP data exists
        mean_connectomes_FND_FW{2,1} = mean_connectomes_FND_FW{2,1} + connectomes{patients_inds_FW(j),2};
        mean_connectomes_FND_FW{2,2} = mean_connectomes_FND_FW{2,2} + connectomes_ml{patients_inds_FW(j),2};
        mean_connectomes_FND_FW{2,3} = mean_connectomes_FND_FW{2,3} + connectomes_fa{patients_inds_FW(j),2};
    end
end
for j = 1:length(patients_inds_nFW)
    mean_connectomes_FND_nFW{1,1} = mean_connectomes_FND_nFW{1,1} + connectomes{patients_inds_nFW(j),1};
    mean_connectomes_FND_nFW{1,2} = mean_connectomes_FND_nFW{1,2} + connectomes_ml{patients_inds_nFW(j),1};
    mean_connectomes_FND_nFW{1,3} = mean_connectomes_FND_nFW{1,3} + connectomes_fa{patients_inds_nFW(j),1};
    if ismember(patients_inds_nFW(j),union(patients_inds_impr,patients_inds_nimpr)) % only if FUP data exists
        mean_connectomes_FND_nFW{2,1} = mean_connectomes_FND_nFW{2,1} + connectomes{patients_inds_nFW(j),2};
        mean_connectomes_FND_nFW{2,2} = mean_connectomes_FND_nFW{2,2} + connectomes_ml{patients_inds_nFW(j),2};
        mean_connectomes_FND_nFW{2,3} = mean_connectomes_FND_nFW{2,3} + connectomes_fa{patients_inds_nFW(j),2};
    end
end
for j = 1:length(patients_inds_impr)
    mean_connectomes_FND_impr{1,1} = mean_connectomes_FND_impr{1,1} + connectomes{patients_inds_impr(j),1};
    mean_connectomes_FND_impr{1,2} = mean_connectomes_FND_impr{1,2} + connectomes_ml{patients_inds_impr(j),1};
    mean_connectomes_FND_impr{1,3} = mean_connectomes_FND_impr{1,3} + connectomes_fa{patients_inds_impr(j),1};
    mean_connectomes_FND_impr{2,1} = mean_connectomes_FND_impr{2,1} + connectomes{patients_inds_impr(j),2};
    mean_connectomes_FND_impr{2,2} = mean_connectomes_FND_impr{2,2} + connectomes_ml{patients_inds_impr(j),2};
    mean_connectomes_FND_impr{2,3} = mean_connectomes_FND_impr{2,3} + connectomes_fa{patients_inds_impr(j),2};
end
for j = 1:length(patients_inds_nimpr)
    mean_connectomes_FND_nimpr{1,1} = mean_connectomes_FND_nimpr{1,1} + connectomes{patients_inds_nimpr(j),1};
    mean_connectomes_FND_nimpr{1,2} = mean_connectomes_FND_nimpr{1,2} + connectomes_ml{patients_inds_nimpr(j),1};
    mean_connectomes_FND_nimpr{1,3} = mean_connectomes_FND_nimpr{1,3} + connectomes_fa{patients_inds_nimpr(j),1};
    mean_connectomes_FND_nimpr{2,1} = mean_connectomes_FND_nimpr{2,1} + connectomes{patients_inds_nimpr(j),2};
    mean_connectomes_FND_nimpr{2,2} = mean_connectomes_FND_nimpr{2,2} + connectomes_ml{patients_inds_nimpr(j),2};
    mean_connectomes_FND_nimpr{2,3} = mean_connectomes_FND_nimpr{2,3} + connectomes_fa{patients_inds_nimpr(j),2};
end
for j = 1:3
    mean_connectomes_FND_FW{1,j} = mean_connectomes_FND_FW{1,j}./length(patients_inds_FW);
    mean_connectomes_FND_FW{1,j} = mean_connectomes_FND_FW{1,j} + transpose(triu(mean_connectomes_FND_FW{1,j},1));
    mean_connectomes_FND_FW{2,j} = mean_connectomes_FND_FW{2,j}./length(ismember(...
        patients_inds_FW,union(patients_inds_impr,patients_inds_nimpr)));
    mean_connectomes_FND_FW{2,j} = mean_connectomes_FND_FW{2,j} + transpose(triu(mean_connectomes_FND_FW{2,j},1));
    mean_connectomes_FND_nFW{1,j} = mean_connectomes_FND_nFW{1,j}./length(patients_inds_nFW);
    mean_connectomes_FND_nFW{1,j} = mean_connectomes_FND_nFW{1,j} + transpose(triu(mean_connectomes_FND_nFW{1,j},1));
    mean_connectomes_FND_nFW{2,j} = mean_connectomes_FND_nFW{2,j}./length(ismember(...
        patients_inds_nFW,union(patients_inds_impr,patients_inds_nimpr)));
    mean_connectomes_FND_nFW{2,j} = mean_connectomes_FND_nFW{2,j} + transpose(triu(mean_connectomes_FND_nFW{2,j},1));
    mean_connectomes_FND_impr{1,j} = mean_connectomes_FND_impr{1,j}./length(patients_inds_impr);
    mean_connectomes_FND_impr{1,j} = mean_connectomes_FND_impr{1,j} + transpose(triu(mean_connectomes_FND_impr{1,j},1));
    mean_connectomes_FND_impr{2,j} = mean_connectomes_FND_impr{2,j}./length(patients_inds_impr);
    mean_connectomes_FND_impr{2,j} = mean_connectomes_FND_impr{2,j} + transpose(triu(mean_connectomes_FND_impr{2,j},1));
    mean_connectomes_FND_nimpr{1,j} = mean_connectomes_FND_nimpr{1,j}./length(patients_inds_nimpr);
    mean_connectomes_FND_nimpr{1,j} = mean_connectomes_FND_nimpr{1,j} + transpose(triu(mean_connectomes_FND_nimpr{1,j},1));
    mean_connectomes_FND_nimpr{2,j} = mean_connectomes_FND_nimpr{2,j}./length(patients_inds_nimpr);
    mean_connectomes_FND_nimpr{2,j} = mean_connectomes_FND_nimpr{2,j} + transpose(triu(mean_connectomes_FND_nimpr{2,j},1));
end

%%
for p = 1 % to wrap plots
    tmp_f11 = figure('Position',[1800 100 900 800]);
    subplot(211);
    imagesc(mean_connectomes_FND_FW{1,3}-mean_connectomes_FND_nFW{1,3});
    axis square; colorbar; caxis([-.1 .1]); xticks(10:10:80); yticks(10:10:80); tmp_f11.CurrentAxes.FontSize = 12;
    title({'Differences in average FA between FND FW and nFW',...
        'patients [T1]'},'FontWeight','normal','Interpreter','none');
    tmp_f11.CurrentAxes.Title.FontSize = 16;
    tmp_f11.CurrentAxes.Title.Position(2) = tmp_f11.CurrentAxes.Title.Position(2)-2;
    xlabel('ROIs','FontSize',14); ylabel('ROIs','FontSize',14);
    subplot(212);
    plot(jUpperTriMatToVec(triu(mean_connectomes_FND_FW{1,3}-mean_connectomes_FND_nFW{1,3},1)));
    ylim([-.15 .15]); xlim([0 3500]);
    xlabel('Vectorized upper triangular part of connectome (84*83/2 ROI-pairs)','FontSize',14,'Interpreter','none');
    ylabel('Difference in FA','FontSize',14);
    % compute additionally the inter-hemispheric and intra-hemispheric average differences
    tmp_M_in = mean_connectomes_FND_FW{1,3}-mean_connectomes_FND_nFW{1,3};
    [tmp_intra, tmp_inter] = compute_hemispheric_metric(tmp_M_in,1);
    title({['Absolute average intra-hemispheric diff. between FW and nFW [T1]: ' num2str(tmp_intra)],...
        ['Absolute average inter-hemispheric diff. between FW and nFW [T1]: ' num2str(tmp_inter)]},...
        'FontWeight','normal','Interpreter','none');
    tmp_f11.CurrentAxes.Title.FontSize = 16;
    tmp_f11.CurrentAxes.Title.Position(2) = tmp_f11.CurrentAxes.Title.Position(2)+.01;
    hold on; tmp_xlim = xlim;
    plot([tmp_xlim(1) tmp_xlim(2)],repmat(...
        mean(jUpperTriMatToVec(triu(tmp_M_in,1))),1,2),'r-.');
    hold off; clear tmp_M_in tmp_xlim;
    if ~exist([export_figs_path 'Figure11_tmp.pdf'],'file')
        %export_fig([export_figs_path 'Figure2c.pdf'],tmp_f3);
        %exportgraphics(tmp_f11,[export_figs_path 'Figure11_tmp.pdf'],...
        %    'BackgroundColor','none','ContentType','vector');
    end
end



%% Compute WD (weighted-degree graph metric) from FA connectomes

tmp_FNDindzonly = find(cell2mat(cellfun(@(x) strcmp(x,'FND'),patients_datatable{:,2},'UniformOutput',false)));
tmp_FNDindzonly = tmp_FNDindzonly(2:end); % remove P48 (no DTI data)

patients_inds_FW_withinFND = zeros(length(patients_inds_FW),1);
patients_inds_nFW_withinFND= zeros(length(patients_inds_nFW),1);
patients_inds_mFND_withinFND = zeros(length(patients_inds_mFND),1);
patients_inds_PNES_withinFND = zeros(length(patients_inds_PNES),1);
patients_inds_PPPD_withinFND = zeros(length(patients_inds_PPPD),1);
for j = 1:length(patients_inds_FW)
    patients_inds_FW_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_FW(j));
end
for j = 1:length(patients_inds_nFW)
    patients_inds_nFW_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_nFW(j));
end
for j = 1:length(patients_inds_mFND)
    patients_inds_mFND_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_mFND(j));
end
for j = 1:length(patients_inds_PNES)
    patients_inds_PNES_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_PNES(j));
end
for j = 1:length(patients_inds_PPPD)
    patients_inds_PPPD_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_PPPD(j));
end

patients_inds_impr_withinFND = zeros(length(patients_inds_impr),1);
patients_inds_nimpr_withinFND = zeros(length(patients_inds_nimpr),1);
for j = 1:length(patients_inds_impr)
    patients_inds_impr_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_impr(j));
end
for j = 1:length(patients_inds_nimpr)
    patients_inds_nimpr_withinFND(j,1) = find(tmp_FNDindzonly==patients_inds_nimpr(j));
end


WD_FND_FW_T1 = cell(length(patients_inds_FW),1);
for j = 1:length(patients_inds_FW)
    WD_FND_FW_T1{j,1} = WD_FND{1,1}{patients_inds_FW_withinFND(j),1};
end
WD_FND_nFW_T1 = cell(length(patients_inds_nFW),1);
for j = 1:length(patients_inds_nFW)
    WD_FND_nFW_T1{j,1} = WD_FND{1,1}{patients_inds_nFW_withinFND(j),1};
end

WD_FND_impr_T1 = cell(length(patients_inds_impr),1);
WD_FND_nimpr_T1= cell(length(patients_inds_nimpr),1);
for j = 1:length(patients_inds_impr)
    WD_FND_impr_T1{j,1} = WD_FND{1,1}{patients_inds_impr_withinFND(j),1};
end
for j = 1:length(patients_inds_nimpr)
    WD_FND_nimpr_T1{j,1} = WD_FND{1,1}{patients_inds_nimpr_withinFND(j),1};
end

WD_FND_impr_T1_mat = reshape(cell2mat(WD_FND_impr_T1),84,31)';
WD_FND_nimpr_T1_mat= reshape(cell2mat(WD_FND_nimpr_T1),84,22)';

%meanWD_FND_FW_T1 = mean(reshape(cell2mat(WD_FND_FW_T1),84,29),2);
%meanWD_FND_nFW_T1= mean(reshape(cell2mat(WD_FND_nFW_T1),84,32),2);

WD_FND_FW_T1_mat = reshape(cell2mat(WD_FND_FW_T1),84,29)';
WD_FND_nFW_T1_mat= reshape(cell2mat(WD_FND_nFW_T1),84,32)';



%% Plot DurationSymptoms X S-FMDRS [T1] & CGI-1 [T1] covariates

for p = 1 % to wrap plots
    tmp_f12 = figure('Position',[1800 300 1100 600]);
    subplot(121);
    tmp_indsofint = find(~isnan(patients_datatable.DurationSymptoms));
    scatter(...
        patients_datatable.DurationSymptoms(tmp_indsofint,1),...
        patients_datatable.("S-FMDRS [T1]")(tmp_indsofint,1),'r*','LineWidth',1.5);
    axis square; %colorbar; %caxis([-.1 .1]); xticks(10:10:80); yticks(10:10:80); tmp_f11.CurrentAxes.FontSize = 12;
    title({'Symptoms Duration \times S-FMDRS','for FND patients [T1]'},'FontWeight','normal','Interpreter','tex');
    tmp_f12.CurrentAxes.Title.FontSize = 16;
    tmp_f12.CurrentAxes.Title.Position(2) = tmp_f12.CurrentAxes.Title.Position(2)+6;
    xlabel('Symptoms Duration (days)','FontSize',14); ylabel('S-FMDRS [T1]','FontSize',14);
    tmp_f12.CurrentAxes.YLabel.Position(1) = tmp_f12.CurrentAxes.YLabel.Position(1)-60;
    tmp_f12.CurrentAxes.XLabel.Position(2) = tmp_f12.CurrentAxes.XLabel.Position(2)-5;
    ylim([-4 50]); xlim([-40 450]);
    tmp_f12.CurrentAxes.Box = 'on'; 
    subplot(122);
    scatter(patients_datatable.DurationSymptoms(tmp_indsofint,1),...
        patients_datatable.("CGI-1 [T1]")(tmp_indsofint,1),'black','Marker','*','LineWidth',1.5);
    axis square; %colorbar; %caxis([-.1 .1]); xticks(10:10:80); yticks(10:10:80); tmp_f11.CurrentAxes.FontSize = 12;
    title({'Symptoms Duration \times CGI-1','for FND patients [T1]'},'FontWeight','normal','Interpreter','tex');
    tmp_f12.CurrentAxes.Title.FontSize = 16;
    tmp_f12.CurrentAxes.Title.Position(2) = tmp_f12.CurrentAxes.Title.Position(2)+1.15;
    xlabel('Symptoms Duration (days)','FontSize',14); ylabel('CGI-1 [T1]','FontSize',14);
    tmp_f12.CurrentAxes.YLabel.Position(1) = tmp_f12.CurrentAxes.YLabel.Position(1)-50;
    tmp_f12.CurrentAxes.XLabel.Position(2) = tmp_f12.CurrentAxes.XLabel.Position(2)-1.95;
    ylim([-.75 7]); xlim([-40 450]);
    tmp_f12.CurrentAxes.Box = 'on'; 
    if ~exist([export_figs_path 'Figure12.pdf'],'file')
        %export_fig([export_figs_path 'Figure2c.pdf'],tmp_f3);
        %exportgraphics(tmp_f12,[export_figs_path 'Figure12.pdf'],...
        %    'BackgroundColor','none','ContentType','vector');
    end
end



%% Probe for FA correlation with DurationSymptoms or SF36 metrics
%  according to the results of Diez et al. (2019):
%  their Stria Terminalis / Fornix bundle consists of the following
%  FA links:
%    - amyR-amyR            [not a sign. hub in our results] | 48
%    - hippR-hippR          [not a sign. hub in our results] | 47
%    - amyL-putamenL**      [putL is a sign. hub]            | (41,38)
%    - hippR-cerebcortexR*  [sign. hub but no links survive] | (47,84)
%    - amyR-entorehinalR    [not sign. hubs in our results]  | (48,54)

for p = 1 % to wrap plots
    %close all;
    tmp_f13 = figure('Position',[1800 50 920 920]);
    %tmp_xlabel = 'SF36-PhysFunc';
    tmp_xlabel = 'Illness duration (months)';
    subplot(2,6,[2 3]);
    tmp_what_to_scatter = patients_datatable.DurationSymptoms(connectomes_FND_ind{1,1},1);
    %tmp_what_to_scatter = patients_datatable.("S-FMDRS [T1]")(connectomes_FND_ind{1,1},1);
    %tmp_what_to_scatter = 100-tmp_what_to_scatter; % Inverted SF36-PhysHealth
    %tmp_what_to_scatter = patients_datatable.("SF36-PhysFunc")(connectomes_FND_ind{1,1},1);
    tmp_probe_FA_stria_fornix = zeros(connectomes_FND_ind{1,2},5); % 5 FA links as defined above / comparing to Diez et al. (2019)
    for j = 1:connectomes_FND_ind{1,2}
        tmp_probe_FA_stria_fornix(j,1) = connectomes_fa{connectomes_FND_ind{1,1}(j,1),1}(10,22);
        tmp_probe_FA_stria_fornix(j,2) = connectomes_fa{connectomes_FND_ind{1,1}(j,1),1}(14,73);
        tmp_probe_FA_stria_fornix(j,3) = connectomes_fa{connectomes_FND_ind{1,1}(j,1),1}(31,64); % because connectomes_FND_ind{1,1} is not symmetrized
        tmp_probe_FA_stria_fornix(j,4) = connectomes_fa{connectomes_FND_ind{1,1}(j,1),1}(47,84);
        tmp_probe_FA_stria_fornix(j,5) = connectomes_fa{connectomes_FND_ind{1,1}(j,1),1}(6,33);
    end
    scatter(tmp_what_to_scatter,tmp_probe_FA_stria_fornix(:,1),'MarkerEdgeColor','bla',...
        'MarkerFaceColor',[141/255 222/255 222/255],'Marker','sq','LineWidth',.75); %[255/255 201/255 135/255] % light orange
    tmp_f13.CurrentAxes.YGrid = 'on'; tmp_f13.CurrentAxes.XGrid = 'on';
    axis square; xticks(0:50:400); %yticks(10:10:80);
    tmp_f13.CurrentAxes.FontSize = 11;
    tmp_f13.CurrentAxes.FontName = 'Basis Grotesque Pro';
    title('Lat Occip L — Post Cing L','FontWeight','normal',...
        'Interpreter','none','FontSize',14,'FontName','Basis Grotesque Pro');
    %tmp_f13.CurrentAxes.Title.Position(2) = tmp_f13.CurrentAxes.Title.Position(2)+.01;
    xlabel(tmp_xlabel,'FontSize',14,'FontName','Basis Grotesque Pro');
    %tmp_f13.CurrentAxes.XLabel.Position(2) = tmp_f13.CurrentAxes.XLabel.Position(2)-.02;
    %tmp_f13.CurrentAxes.YLabel.Position(1) = tmp_f13.CurrentAxes.YLabel.Position(1)-20;
    ylabel('FA','FontSize',14,'FontName','Basis Grotesque Pro');
    tmp_f13.CurrentAxes.Box = 'on'; xlim([-20 420]); ylim([.2 .8]);
    subplot(2,6,[4 5]);
    scatter(tmp_what_to_scatter,tmp_probe_FA_stria_fornix(:,2),'MarkerEdgeColor','bla',...
        'MarkerFaceColor',[141/255 222/255 222/255],'Marker','sq','LineWidth',.75);
    tmp_f13.CurrentAxes.YGrid = 'on';
    axis square; xticks(0:10:50); %yticks(10:10:80);
    tmp_f13.CurrentAxes.FontSize = 11;
    tmp_f13.CurrentAxes.FontName = 'Basis Grotesque Pro';
    title('Precuneus R — Middle temporal L','FontWeight','normal',...
        'Interpreter','none','FontSize',16,'FontName','Basis Grotesque Pro');
    %tmp_f13.CurrentAxes.Title.Position(2) = tmp_f13.CurrentAxes.Title.Position(2)+.11;
    xlabel(tmp_xlabel,'FontSize',14,'FontName','Basis Grotesque Pro');
    %tmp_f13.CurrentAxes.XLabel.Position(2) = tmp_f13.CurrentAxes.XLabel.Position(2)-.19;
    ylabel('FA','FontSize',14,'FontName','Basis Grotesque Pro'); tmp_f13.CurrentAxes.Box = 'on'; xlim([-10 50]); ylim([0.3 0.8]);
    subplot(2,6,[7 8]);
    scatter(tmp_what_to_scatter,tmp_probe_FA_stria_fornix(:,3),'MarkerEdgeColor','bla',...
        'MarkerFaceColor',[141/255 222/255 222/255],'Marker','sq','LineWidth',.75);
    tmp_f13.CurrentAxes.YGrid = 'on';
    axis square; %xticks(0:20:100); %yticks(10:10:80);
    tmp_f13.CurrentAxes.FontSize = 11;
    tmp_f13.CurrentAxes.FontName = 'Basis Grotesque Pro';
    title('Amy L - Putamen L','FontWeight','normal','Interpreter','none','FontSize',16);
    tmp_f13.CurrentAxes.Title.Position(2) = tmp_f13.CurrentAxes.Title.Position(2)+.07;
    xlabel(tmp_xlabel,'FontSize',14);
    ylabel('FA','FontSize',14); tmp_f13.CurrentAxes.Box = 'on'; xlim([-2 8]); ylim([0.1 .5]);
    tmp_f13.CurrentAxes.YLabel.Position(1) = tmp_f13.CurrentAxes.YLabel.Position(1)-8;
    tmp_f13.CurrentAxes.XLabel.Position(2) = tmp_f13.CurrentAxes.XLabel.Position(2)-.025;
    subplot(2,6,[9 10]);
    scatter(tmp_what_to_scatter,tmp_probe_FA_stria_fornix(:,4),'MarkerEdgeColor','bla',...
        'MarkerFaceColor',[141/255 222/255 222/255],'Marker','sq','LineWidth',.75);
    tmp_f13.CurrentAxes.YGrid = 'on';
    axis square; %xticks(0:20:100); %yticks(10:10:80);
    tmp_f13.CurrentAxes.FontSize = 11;
    tmp_f13.CurrentAxes.FontName = 'Basis Grotesque Pro';
    title('Hipp R - CerebCortex R','FontWeight','normal','Interpreter','none','FontSize',16);
    tmp_f13.CurrentAxes.Title.Position(2) = tmp_f13.CurrentAxes.Title.Position(2)+.0425;
    xlabel(tmp_xlabel,'FontSize',14);
    tmp_f13.CurrentAxes.XLabel.Position(2) = tmp_f13.CurrentAxes.XLabel.Position(2)-.083;
    ylabel('','FontSize',14); tmp_f13.CurrentAxes.Box = 'on'; xlim([-5 45]); ylim([0.4 .65]); %xlim([-10 110]);
    subplot(2,6,[11 12]);
    scatter(tmp_what_to_scatter,tmp_probe_FA_stria_fornix(:,5),'MarkerEdgeColor','bla',...
        'MarkerFaceColor',[141/255 222/255 222/255],'Marker','sq','LineWidth',.75);
    tmp_f13.CurrentAxes.YGrid = 'on';
    axis square; %xticks(0:20:100); %yticks(10:10:80);
    tmp_f13.CurrentAxes.FontSize = 11;
    tmp_f13.CurrentAxes.FontName = 'Basis Grotesque Pro';
    title('Amy R - Entorehinal R','FontWeight','normal','Interpreter','none','FontSize',16);
    tmp_f13.CurrentAxes.Title.Position(2) = tmp_f13.CurrentAxes.Title.Position(2)+.066;
    xlabel(tmp_xlabel,'FontSize',14);
    tmp_f13.CurrentAxes.XLabel.Position(2) = tmp_f13.CurrentAxes.XLabel.Position(2)-.075;
    ylabel('','FontSize',14); tmp_f13.CurrentAxes.Box = 'on';
    xlim([-5 105]); ylim([0.3 .7]);
    if ~exist([export_figs_path 'Figure13_e.pdf'],'file')
        %export_fig([export_figs_path 'Figure13_tmp.pdf'],tmp_f3);
        %exportgraphics(tmp_f13,[export_figs_path 'Figure13_e.pdf'],...
        %    'BackgroundColor','none','ContentType','vector');
    end
end



%% Widesearch best correlates [best linear fit]
%   across available WDs / FDs for 85 FND patients
%   update: as well as across all 75 HCs and
%   across all HCs + FND patients taken together for clin. variables
%   that are available for both groups (e.g., BDI, STAI, etc.)

%widesearch_output_corr_fa = zeros(length(connectomes_fa{1,1})*(length(connectomes_fa{1,1})+1)/2,9);

%%% [at T1]
widesearch_vars_FND(1:connectomes_FND_ind{1,2},1:12) = nan; % init
widesearch_vars_FND(:,1) = patients_datatable.DurationSymptoms(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,2) = patients_datatable.("S-FMDRS [T1]")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,3) = patients_datatable.("CGI-1 [T1]")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,4) = patients_datatable.("SF36-PhysHealth")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,5) = patients_datatable.("SF36-MentalHealth")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,6) = patients_datatable.("SF36-GenHealth")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,7) = patients_datatable.("SF36-PhysFunc")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,8) = patients_datatable.("AUCi")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,9) = patients_datatable.("CAR_AUCi")(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,10) = patients_datatable.bdi(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,11) = patients_datatable.stai1(connectomes_FND_ind{1,1},1);
widesearch_vars_FND(:,12) = patients_datatable.stai2(connectomes_FND_ind{1,1},1);

%%% [at FUP/between FUP and T1]
widesearch_vars2(:,1) = patients_datatable.("∆S-FMDRS")(connectomes_FND_ind{2,1},1);

%%% update: add widesearch for HCs and for both groups together for available clin. variables
widesearch_vars_HC(1:connectomes_HC_ind{1,2},1:9) = nan;
widesearch_vars_HC(:,1) = patients_datatable.("SF36-PhysHealth")(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,2) = patients_datatable.("SF36-MentalHealth")(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,3) = patients_datatable.("SF36-GenHealth")(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,4) = patients_datatable.("SF36-PhysFunc")(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,5) = patients_datatable.("AUCi")(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,6) = patients_datatable.("CAR_AUCi")(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,7) = patients_datatable.bdi(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,8) = patients_datatable.stai1(connectomes_HC_ind{1,1},1);
widesearch_vars_HC(:,9) = patients_datatable.stai2(connectomes_HC_ind{1,1},1);

%%% update: add widesearch for both combined on those clin. vars available
widesearch_vars_combined(1:length([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}]),1:9) = nan;
widesearch_vars_combined(:,1) = patients_datatable.("SF36-PhysHealth")([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,2) = patients_datatable.("SF36-MentalHealth")([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,3) = patients_datatable.("SF36-GenHealth")([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,4) = patients_datatable.("SF36-PhysFunc")([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,5) = patients_datatable.("AUCi")([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,6) = patients_datatable.("CAR_AUCi")([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,7) = patients_datatable.bdi([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,8) = patients_datatable.stai1([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);
widesearch_vars_combined(:,9) = patients_datatable.stai2([connectomes_HC_ind{1,1};connectomes_FND_ind{1,1}],1);

%links_tmp_FND = reshape(links_FND_fa,connectomes_FND_ind{1,2},[]); % reshapes only for significant ROIs (85 FND subj x all relevant links, 19x84)
links_tmp_FND = reshape(Y_reg2((connectomes_HC_ind{1,2}+1):end,:,:),connectomes_FND_ind{1,2},[]);
widesearch_output_FND_corr_fa      = zeros(length(links_tmp_FND),size(widesearch_vars_FND,2));
widesearch_output_FND_corr_fa_pval = zeros(length(links_tmp_FND),size(widesearch_vars_FND,2));
%links_tmp_HC = reshape(links_HC_fa,connectomes_HC_ind{1,2},[]); % same for HCs only (75 HCs x all relevant links)
links_tmp_HC = reshape(Y_reg2(1:connectomes_HC_ind{1,2},:,:),connectomes_HC_ind{1,2},[]);
widesearch_output_HC_corr_fa      = zeros(length(links_tmp_HC),size(widesearch_vars_HC,2));
widesearch_output_HC_corr_fa_pval = zeros(length(links_tmp_HC),size(widesearch_vars_HC,2));
links_tmp_combined = [links_tmp_HC;links_tmp_FND]; % same for combined
widesearch_output_combined_corr_fa      = zeros(length(links_tmp_combined),size(widesearch_vars_combined,2));
widesearch_output_combined_corr_fa_pval = zeros(length(links_tmp_combined),size(widesearch_vars_combined,2));

for q = 1:size(widesearch_vars_FND,2)
    for k = 1:length(links_tmp_FND)
        [widesearch_output_FND_corr_fa(k,q),...
            widesearch_output_FND_corr_fa_pval(k,q)] = corr(links_tmp_FND(:,k),widesearch_vars_FND(:,q),'rows','pairwise');
    end
end
for q = 1:size(widesearch_vars_HC,2)
    for k = 1:length(links_tmp_HC)
        [widesearch_output_HC_corr_fa(k,q),...
            widesearch_output_HC_corr_fa_pval(k,q)] = corr(links_tmp_HC(:,k),widesearch_vars_HC(:,q),'rows','pairwise');
    end
end
for q = 1:size(widesearch_vars_combined,2)
    for k = 1:length(links_tmp_combined)
        [widesearch_output_combined_corr_fa(k,q),...
            widesearch_output_combined_corr_fa_pval(k,q)] = corr(links_tmp_combined(:,k),widesearch_vars_combined(:,q),'rows','pairwise');
    end
end

tmp_sorted_widesearch_output_FND_corr_fa_pval = widesearch_output_FND_corr_fa_pval;
tmp_sorted_widesearch_output_FND_corr_fa_indz = widesearch_output_FND_corr_fa_pval;
for q = 1:size(widesearch_output_FND_corr_fa_pval,2)
    [tmp_sorted_widesearch_output_FND_corr_fa_pval(:,q),...
        tmp_sorted_widesearch_output_FND_corr_fa_indz(:,q)] = sort(widesearch_output_FND_corr_fa_pval(:,q));
end
tmp_sorted_widesearch_output_HC_corr_fa_pval = widesearch_output_HC_corr_fa_pval;
tmp_sorted_widesearch_output_HC_corr_fa_indz = widesearch_output_HC_corr_fa_pval;
for q = 1:size(widesearch_output_HC_corr_fa_pval,2)
    [tmp_sorted_widesearch_output_HC_corr_fa_pval(:,q),...
        tmp_sorted_widesearch_output_HC_corr_fa_indz(:,q)] = sort(widesearch_output_HC_corr_fa_pval(:,q));
end
tmp_sorted_widesearch_output_combined_corr_fa_pval = widesearch_output_combined_corr_fa_pval;
tmp_sorted_widesearch_output_combined_corr_fa_indz = widesearch_output_combined_corr_fa_pval;
for q = 1:size(widesearch_output_combined_corr_fa_pval,2)
    [tmp_sorted_widesearch_output_combined_corr_fa_pval(:,q),...
        tmp_sorted_widesearch_output_combined_corr_fa_indz(:,q)] = sort(widesearch_output_combined_corr_fa_pval(:,q));
end

%%% add R2 computation for all cases (distribution is plotted then as well)
widesearch_output_FND_corr_fa_R2 = widesearch_output_FND_corr_fa_pval;
widesearch_output_HC_corr_fa_R2 = widesearch_output_HC_corr_fa_pval;
widesearch_output_combined_corr_fa_R2 = widesearch_output_combined_corr_fa_pval;
for q = 1:size(widesearch_output_FND_corr_fa_R2,2) % 12 clin vars
    for l = 1:size(links_tmp_FND,2)
        tmp_PP = polyfit(widesearch_vars_FND(:,q),links_tmp_FND(:,l),1);
        tmp_yfit = polyval(tmp_PP,widesearch_vars_FND(:,q));
        widesearch_output_FND_corr_fa_R2(l,q) = ...
            1 - sum((links_tmp_FND(:,l)-tmp_yfit).^2)/(...
        (length(links_tmp_FND(:,l))-1) * var(links_tmp_FND(:,l)));
    end
end
for q = 1:size(widesearch_output_HC_corr_fa_R2,2) % 9 clin vars
    for l = 1:size(links_tmp_HC,2)
        tmp_PP = polyfit(widesearch_vars_HC(:,q),links_tmp_HC(:,l),1);
        tmp_yfit = polyval(tmp_PP,widesearch_vars_HC(:,q));
        widesearch_output_HC_corr_fa_R2(l,q) = ...
            1 - sum((links_tmp_HC(:,l)-tmp_yfit).^2)/(...
        (length(links_tmp_HC(:,l))-1) * var(links_tmp_HC(:,l)));
        tmp_PP = polyfit(widesearch_vars_combined(:,q),links_tmp_combined(:,l),1);
        tmp_yfit = polyval(tmp_PP,widesearch_vars_combined(:,q));
        widesearch_output_combined_corr_fa_R2(l,q) = ...
            1 - sum((links_tmp_combined(:,l)-tmp_yfit).^2)/(...
        (length(links_tmp_combined(:,l))-1) * var(links_tmp_combined(:,l)));
    end
end

q_val = 0.05; % [if not defined already]
%q_val_strict = 0.01;
sign_links_nb_fa_all = zeros(size(widesearch_output_FND_corr_fa_pval,2),3); % store the # of sign. links surviving each fitting (HC, FND, HC+FND (combined))
for q = 1:size(widesearch_output_HC_corr_fa_pval,2)
    [~,tmp_pval_1] = fdr_bh(widesearch_output_HC_corr_fa_pval(:,q),q_val,'pdep','no');
    sign_links_nb_fa_all(q,1) = length(find(widesearch_output_HC_corr_fa_pval(:,q)<=tmp_pval_1)); % find # of links
    [~,tmp_pval_1] = fdr_bh(widesearch_output_combined_corr_fa_pval(:,q),q_val,'pdep','no');
    sign_links_nb_fa_all(q,3) = length(find(widesearch_output_combined_corr_fa_pval(:,q)<=tmp_pval_1)); % find # of links
end
for q = 1:size(widesearch_output_FND_corr_fa_pval,2)
    [~,tmp_pval_1] = fdr_bh(widesearch_output_FND_corr_fa_pval(:,q),q_val,'pdep','no');
    sign_links_nb_fa_all(q,2) = length(find(widesearch_output_FND_corr_fa_pval(:,q)<=tmp_pval_1)); % find # of links
end
clear q tmp_pval_1;



%% Plots
tmp_f14 = figure('Position',[-1250 100 960 520]);
tmp_clrmap = colormap(turbo(size(widesearch_output_FND_corr_fa_pval,2)));
plot(tmp_sorted_widesearch_output_FND_corr_fa_pval(:,1),'LineWidth',1.3,'Color',tmp_clrmap(1,:)); % DurationSymptoms (for FND)
hold(tmp_f14.CurrentAxes,'on');
for q = 2:size(widesearch_output_FND_corr_fa_pval,2)
    plot(tmp_sorted_widesearch_output_FND_corr_fa_pval(:,q),'LineWidth',1.3,'Color',tmp_clrmap(q,:));
    %hold(tmp_f14.CurrentAxes,'on');
end
tmp_f14.CurrentAxes.FontSize = 12;
tmp_f14.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-50 1650 -.1 1.1]);
xticks(0:200:1600); yticks(0:0.1:1); tmp_f14.CurrentAxes.YGrid = 'on';
plot([0 1600],[.05 .05],'bla-','LineWidth',1.3);
text(1500,.075,'\itp\rm < .05','FontSize',12,'FontName','Basis Grotesque Pro');
legend([' Illness Duration {' num2str(sign_links_nb_fa_all(1,2)) '}'],...
    [' S-FMDRS {' num2str(sign_links_nb_fa_all(2,2)) '}'],...
    [' CGI {' num2str(sign_links_nb_fa_all(3,2)) '}'],...
    [' SF36-PhysHealth {' num2str(sign_links_nb_fa_all(4,2)) '}'],...
    [' SF36-MentalHealth {' num2str(sign_links_nb_fa_all(5,2)) '}'],...
    [' SF36-GenHealth {' num2str(sign_links_nb_fa_all(6,2)) '}'],...
    [' SF36-PhysFunc {' num2str(sign_links_nb_fa_all(7,2)) '}'],[' AUCi [cortisol] {' num2str(sign_links_nb_fa_all(8,2)) '}'],...
    [' CAR_AUCi [cortisol] {' num2str(sign_links_nb_fa_all(9,2)) '}'],...
    [' BDI {' num2str(sign_links_nb_fa_all(10,2)) '}'],[' STAI-S {' num2str(sign_links_nb_fa_all(11,2)) '}'],...
    [' STAI-T {' num2str(sign_links_nb_fa_all(12,2)) '}'],...
    'FontSize',12,'FontName','Basis Grotesque Pro',...
    'Interpreter','None','Location','Best');
text(1100,.26,'CGI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(3,:));
text(1200,.165,'S-FMDRS','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(2,:));
text(880,.2,'SF36-PhysFunc','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(7,:));
text(600,.33,'BDI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(10,:));
text(600,.41,'STAI-T','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(12,:));
tmp_f14.CurrentAxes.YAxis.Label.String = '\itp\rm-value';
tmp_f14.CurrentAxes.YAxis.Label.FontSize = 14;
tmp_f14.CurrentAxes.YAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f14.CurrentAxes.XAxis.Label.String = 'Sorted links (FA)';
tmp_f14.CurrentAxes.XAxis.Label.FontSize = 14;
tmp_f14.CurrentAxes.XAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f14.CurrentAxes.Title.String = ...
    ['Sorted links (FA) between all significantly altered gray matter ROIs in FND patients only (q_{FDR} = ' num2str(q_val) ')'];
%tmp_f14.CurrentAxes.Title.String = ...
%    'Sorted links (FA) between all significantly altered gray matter ROIs in FND \itvs\rm HCs';
tmp_f14.CurrentAxes.Title.FontSize = 16;
tmp_f14.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f14.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';


tmp_f14_a = figure('Position',[-1250 100 960 520]);
%tmp_clrmap = colormap(turbo(size(widesearch_output_HC_corr_fa_pval,2))); % below: use same as for FND for consistency
plot(tmp_sorted_widesearch_output_HC_corr_fa_pval(:,1),'LineWidth',1.3,'Color',tmp_clrmap(4,:)); % "SF36-PhysHealth" (for HCs)
hold(tmp_f14_a.CurrentAxes,'on');
for q = 2:size(widesearch_output_HC_corr_fa_pval,2)
    plot(tmp_sorted_widesearch_output_HC_corr_fa_pval(:,q),'LineWidth',1.3,'Color',tmp_clrmap(q+3,:));
    %hold(tmp_f14_a.CurrentAxes,'on');
end
tmp_f14_a.CurrentAxes.FontSize = 12;
tmp_f14_a.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-50 1650 -.1 1.1]);
xticks(0:200:1600); yticks(0:0.1:1); tmp_f14_a.CurrentAxes.YGrid = 'on';
plot([0 1600],[.05 .05],'bla-','LineWidth',1.3);
text(1500,.075,'\itp\rm < .05','FontSize',12,'FontName','Basis Grotesque Pro');
legend([' SF36-PhysHealth {' num2str(sign_links_nb_fa_all(1,1)) '}'],...
    [' SF36-MentalHealth {' num2str(sign_links_nb_fa_all(2,1)) '}'],...
    [' SF36-GenHealth {' num2str(sign_links_nb_fa_all(3,1)) '}'],...
    [' SF36-PhysFunc {' num2str(sign_links_nb_fa_all(4,1)) '}'],...
    [' AUCi [cortisol] {' num2str(sign_links_nb_fa_all(5,1)) '}'],...
    [' CAR_AUCi [cortisol] {' num2str(sign_links_nb_fa_all(6,1)) '}'],...
    [' BDI {' num2str(sign_links_nb_fa_all(7,1)) '}'],...
    [' STAI-S {' num2str(sign_links_nb_fa_all(8,1)) '}'],[' STAI-T {' num2str(sign_links_nb_fa_all(9,1)) '}'],...
    'FontSize',12,'FontName','Basis Grotesque Pro',...
    'Interpreter','None','Location','Best');
%text(1100,.26,'CGI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(3,:));
%text(1200,.165,'S-FMDRS','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(2,:));
%text(880,.2,'SF36-PhysFunc','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(7,:));
%text(600,.33,'BDI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(10,:));
%text(600,.41,'STAI-T','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(12,:));
tmp_f14_a.CurrentAxes.YAxis.Label.String = '\itp\rm-value';
tmp_f14_a.CurrentAxes.YAxis.Label.FontSize = 14;
tmp_f14_a.CurrentAxes.YAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f14_a.CurrentAxes.XAxis.Label.String = 'Sorted links (FA)';
tmp_f14_a.CurrentAxes.XAxis.Label.FontSize = 14;
tmp_f14_a.CurrentAxes.XAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f14_a.CurrentAxes.Title.String = ...
    ['Sorted links (FA) between all significantly altered gray matter ROIs in HCs only (q_{FDR} = ' num2str(q_val) ')'];
tmp_f14_a.CurrentAxes.Title.FontSize = 16;
tmp_f14_a.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f14_a.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';


tmp_f14_b = figure('Position',[-1250 100 960 520]);
%tmp_clrmap = colormap(turbo(size(widesearch_output_HC_corr_fa_pval,2))); % below: use same as for FND for consistency
plot(tmp_sorted_widesearch_output_combined_corr_fa_pval(:,1),'LineWidth',1.3,'Color',tmp_clrmap(4,:)); % "SF36-PhysHealth" (for HCs)
hold(tmp_f14_b.CurrentAxes,'on');
for q = 2:size(widesearch_output_combined_corr_fa_pval,2)
    plot(tmp_sorted_widesearch_output_combined_corr_fa_pval(:,q),'LineWidth',1.3,'Color',tmp_clrmap(q+3,:));
    %hold(tmp_f14_b.CurrentAxes,'on');
end
tmp_f14_b.CurrentAxes.FontSize = 12;
tmp_f14_b.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-50 1650 -.1 1.1]);
xticks(0:200:1600); yticks(0:0.1:1); tmp_f14_b.CurrentAxes.YGrid = 'on';
plot([0 1600],[.05 .05],'bla-','LineWidth',1.3);
text(1500,.075,'\itp\rm < .05','FontSize',12,'FontName','Basis Grotesque Pro');
legend([' SF36-PhysHealth {' num2str(sign_links_nb_fa_all(1,3)) '}'],...
    [' SF36-MentalHealth {' num2str(sign_links_nb_fa_all(2,3)) '}'],...
    [' SF36-GenHealth {' num2str(sign_links_nb_fa_all(3,3)) '}'], ...
    [' SF36-PhysFunc {' num2str(sign_links_nb_fa_all(4,3)) '}'], ...
    [' AUCi [cortisol] {' num2str(sign_links_nb_fa_all(5,3)) '}'],...
    [' CAR_AUCi [cortisol] {' num2str(sign_links_nb_fa_all(6,3)) '}'], ...
    [' BDI {' num2str(sign_links_nb_fa_all(7,3)) '}'], ...
    [' STAI-S {' num2str(sign_links_nb_fa_all(8,3)) '}'], ...
    [' STAI-T {' num2str(sign_links_nb_fa_all(9,3)) '}'],...
    'FontSize',12,'FontName','Basis Grotesque Pro',...
    'Interpreter','None','Location','Best');
%text(1100,.26,'CGI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(3,:));
text(1370,.4,'SF36-GenHealth','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(6,:));
text(1370,.2,'SF36-PhysFunc','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(7,:));
text(1390,.5,'BDI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(10,:));
%text(600,.41,'STAI-T','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(12,:));
tmp_f14_b.CurrentAxes.YAxis.Label.String = '\itp\rm-value';
tmp_f14_b.CurrentAxes.YAxis.Label.FontSize = 14;
tmp_f14_b.CurrentAxes.YAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f14_b.CurrentAxes.XAxis.Label.String = 'Sorted links (FA)';
tmp_f14_b.CurrentAxes.XAxis.Label.FontSize = 14;
tmp_f14_b.CurrentAxes.XAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f14_b.CurrentAxes.Title.String = ...
    ['Sorted links (FA) between all significantly altered gray matter ROIs in HCs and FND patients (q_{FDR} = ' num2str(q_val) ')'];
tmp_f14_b.CurrentAxes.Title.FontSize = 16;
tmp_f14_b.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f14_b.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';

%
%exportgraphics(tmp_f14, 'Figures/Figure14_q001_age_gender_BDI_STAI-T_reg.png', 'Resolution', '300');
%exportgraphics(tmp_f14_a, 'Figures/Figure14a_q001_age_gender_BDI_STAI-T_reg.png', 'Resolution', '300');
%exportgraphics(tmp_f14_b, 'Figures/Figure14b_q001_age_gender_BDI_STAI-T_reg.png', 'Resolution', '300');


%%
%close all;
tmp_f16 = figure('Position',[-1250 100 960 520]);
scatter(widesearch_vars_FND(:,1),links_tmp_FND(:,190),30,tmp_clrmap(1,:),'filled');
hold(tmp_f16.CurrentAxes,'on');
tmp_f16.CurrentAxes.FontSize = 12;
tmp_f16.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-10 410 .3 .7]);
xticks(0:40:400); yticks(.35:.05:.65); tmp_f16.CurrentAxes.YGrid = 'on';
title('Lat Occip L — Post Cing L connectivity (FA) \itvs\rm illness duration',...
    'FontSize',16,'FontName','Basis Grotesque Pro','FontWeight','Normal');
xlabel('Illness duration (months)','FontSize',14,'FontName','Basis Grotesque Pro');
ylabel('FA','FontSize',14,'FontName','Basis Grotesque Pro');
tmp_PP = polyfit(widesearch_vars_FND(:,1),links_tmp_FND(:,190),1);
%tmp_yfit = P(1)*x+P(2);
plot(widesearch_vars_FND(:,1),tmp_PP(1)*widesearch_vars_FND(:,1)+tmp_PP(2),...
    '-','Color',tmp_clrmap(1,:),'LineWidth',1.3);
legend(' FND patient',' Linear interp.','FontSize',14,'FontName','Basis Grotesque Pro','Location','Best');
text(300,.35,{'\itp\rm < 0.0002','(FDR-corrected)'},'FontSize',10,...
    'FontName','Basis Grotesque Pro','Color',tmp_clrmap(1,:));
clear tmp_PP;



%% Filter the significantly correlated links with CGI or other variables
% that are also significantly different between FND and HCs (as per h5__fdr
% above)

q_val = 0.05;
%dataMat_SFMDRS = zeros(size(dataMat)); dataMat_CGI = zeros(size(dataMat)); dataMat_SF36PhysFunc = zeros(size(dataMat));
dataMat_2{2,1} = zeros(size(dataMat)); dataMat_2{3,1} = zeros(size(dataMat)); dataMat_2{7,1} = zeros(size(dataMat));
tmp_refmat = reshape(1:length(tmp_sorted_widesearch_output_FND_corr_fa_indz),size(dataMat,1),length(sign_ROIs));
for tmp_IOI = [2 3 7] % for S-FMDRS, CGI, and SF-36 PhysFunc which are the only ones with lots of sign. corr. with FA
    [~,tmp_TEST2] = fdr_bh(widesearch_output_FND_corr_fa_pval(:,tmp_IOI),q_val,'pdep','no');
    tmp_qmax = length(find(widesearch_output_FND_corr_fa_pval(:,tmp_IOI)<=tmp_TEST2));
    for q = 1:tmp_qmax
        [tmpi,tmpj] = ind2sub(size(tmp_refmat),tmp_sorted_widesearch_output_FND_corr_fa_indz(q,tmp_IOI));
        if h5__fdr(tmpi,sign_ROIs(tmpj))
            dataMat_2{tmp_IOI,1}(cell2mat(DesikanKilliany_atlas_connectogram_labels(:,2))==...
                sign_ROIs(tmpj),cell2mat(DesikanKilliany_atlas_connectogram_labels(:,2))==tmpi) = 1;
        end
    end
end



%% Plot the connectograms with clinical correlations of FA
%  as 'filtered' versions of the connectogram between HCs and FND

if exist('tmp_f18','var') && ishandle(tmp_f18), close(tmp_f18); end
tmp_f18 = figure('Position',[325 6 1000 1000]);
if ~exist('DesikanKilliany_atlas_connectogram_labels','var')
    load('DesikanKilliany_atlas_connectogram_labels.mat');
end

subplot(221);
BCC = biChordChart(dataMat, 'Label', NameList, ...
    'Arrow', 'on', 'CData', repmat([.2 .2 .2],size(dataMat,1),1)); % plot original dataMat
BCC = BCC.draw(); 
BCC.labelRotate('on');
BCC.setLabelRadius(1.2);
BCC.setFont('FontName', 'Basis Grotesque Pro', ...
    'FontSize', 8, 'Color', [.2 .2 .2]);

% customized coloring [according to Diez et al., 2021?]
tmp_alllabels = findobj(tmp_f18.CurrentAxes,'Tag','BiChordLabel');
tmp_alllabels = flipud(tmp_alllabels); % [because those are retrieved in reverse order]
tmp_colorpalette_ChordN = ...
    ['#fabebe'; '#3cb44b'; '#ffe119'; ... % InfPar InfTemp LatOcc (left)
     '#4363d8'; '#f58231'; '#911eb4'; ...
     '#46f0f0'; '#f032e6'; '#bcf60c'; ...
     '#e6194b'; '#008080'; '#e6beff'; ...
     '#9a6324'; '#42d4f4'; '#800000'; ...
     '#808000'];% '#808000'; '#ffd8b1'; ...
     %'#000075'; '#808080'; '#ffffff'; '#000000'];

tmp_hardcodedcoloringindices = [4 5 7:19]; % regions with j==1 2 3 or 6 are w/out sign. links
tmp_cc = 1;
for j = 1:length(sign_ROIs)
    tmp_ind10 = find(cell2mat(...
        DesikanKilliany_atlas_connectogram_labels(:,2))==sign_ROIs(j));
    tmp_alllabels(tmp_ind10).FontWeight = 'Bold';
    tmp_alllabels(tmp_ind10).FontSize = 12;
    if sum(dataMat(tmp_ind10,:))>0
        tmp_alllabels(tmp_ind10).Color = tmp_colorpalette_ChordN(tmp_cc,:);
        BCC.setSquareN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:));
        BCC.setChordN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:)); % sanity check
        tmp_cc = tmp_cc + 1;
    else
        tmp_alllabels(tmp_ind10).Color = '#a9a9a9';
        BCC.setSquareN(tmp_ind10,'FaceColor','#a9a9a9');
    end
end

subplot(222);
BCC = biChordChart(dataMat_2{2,1}, 'Label', NameList, ...
    'Arrow', 'on', 'CData', repmat([.2 .2 .2],size(dataMat_2{2,1},1),1));
BCC = BCC.draw(); 
BCC.labelRotate('on');
BCC.setLabelRadius(1.2);
BCC.setFont('FontName', 'Basis Grotesque Pro', ...
    'FontSize', 8, 'Color', [.2 .2 .2]);

tmp_alllabels = findobj(tmp_f18.CurrentAxes,'Tag','BiChordLabel');
tmp_alllabels = flipud(tmp_alllabels); % [because those are retrieved in reverse order]

tmp_f18.CurrentAxes.Title.String = 'S-FMDRS';
tmp_f18.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f18.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';
tmp_f18.CurrentAxes.Title.FontSize = 16;
tmp_f18.CurrentAxes.Title.Position(2) = tmp_f18.CurrentAxes.Title.Position(2)+.4;

tmp_cc = 1;
for j = 1:length(sign_ROIs)
    tmp_ind10 = find(cell2mat(...
        DesikanKilliany_atlas_connectogram_labels(:,2))==sign_ROIs(j));
    tmp_alllabels(tmp_ind10).FontWeight = 'Bold';
    tmp_alllabels(tmp_ind10).FontSize = 12;
    if sum(dataMat_2{2,1}(tmp_ind10,:))>0 && ismember(j,tmp_hardcodedcoloringindices)
        tmp_alllabels(tmp_ind10).Color = tmp_colorpalette_ChordN(tmp_cc,:);
        BCC.setSquareN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:));
        BCC.setChordN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:)); % sanity check
        %tmp_cc = tmp_cc+1; % see comment below
    else
        tmp_alllabels(tmp_ind10).Color = '#a9a9a9';
        BCC.setSquareN(tmp_ind10,'FaceColor','#a9a9a9');
    end
    if sum(dataMat(tmp_ind10,:))>0
        tmp_cc = tmp_cc+1; % to keep the same coloring scheme as for dataMat
    end
end

subplot(223);
BCC = biChordChart(dataMat_2{3,1}, 'Label', NameList, ...
    'Arrow', 'on', 'CData', repmat([.2 .2 .2],size(dataMat_2{3,1},1),1));
BCC = BCC.draw(); 
BCC.labelRotate('on');
BCC.setLabelRadius(1.2);
BCC.setFont('FontName', 'Basis Grotesque Arabic Pro', ...
    'FontSize', 8, 'Color', [.2 .2 .2]);

tmp_alllabels = findobj(tmp_f18.CurrentAxes,'Tag','BiChordLabel');
tmp_alllabels = flipud(tmp_alllabels); % [because those are retrieved in reverse order]

tmp_f18.CurrentAxes.Title.String = 'CGI';
tmp_f18.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f18.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';
tmp_f18.CurrentAxes.Title.FontSize = 16;
tmp_f18.CurrentAxes.Title.Position(2) = tmp_f18.CurrentAxes.Title.Position(2)+.4;

tmp_cc = 1;
for j = 1:length(sign_ROIs)
    tmp_ind10 = find(cell2mat(...
        DesikanKilliany_atlas_connectogram_labels(:,2))==sign_ROIs(j));
    tmp_alllabels(tmp_ind10).FontWeight = 'Bold';
    tmp_alllabels(tmp_ind10).FontSize = 12;
    if sum(dataMat_2{3,1}(tmp_ind10,:))>0 && ismember(j,tmp_hardcodedcoloringindices)
        tmp_alllabels(tmp_ind10).Color = tmp_colorpalette_ChordN(tmp_cc,:);
        BCC.setSquareN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:));
        BCC.setChordN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:)); % sanity check
        %tmp_cc = tmp_cc+1;
    else
        tmp_alllabels(tmp_ind10).Color = '#a9a9a9';
        BCC.setSquareN(tmp_ind10,'FaceColor','#a9a9a9');
    end
    if sum(dataMat(tmp_ind10,:))>0
        tmp_cc = tmp_cc+1; % to keep the same coloring scheme as for dataMat
    end
end

subplot(224);
BCC = biChordChart(dataMat_2{7,1}, 'Label', NameList, ...
    'Arrow', 'on', 'CData', repmat([.2 .2 .2],size(dataMat_2{7,1},1),1));
BCC = BCC.draw(); 
BCC.labelRotate('on');
BCC.setLabelRadius(1.2);
BCC.setFont('FontName', 'Basis Grotesque Pro', ...
    'FontSize', 8, 'Color', [.2 .2 .2]);

tmp_alllabels = findobj(tmp_f18.CurrentAxes,'Tag','BiChordLabel');
tmp_alllabels = flipud(tmp_alllabels); % [because those are retrieved in reverse order]

title('SF36-PhysFunc','FontSize',16,'FontName','Basis Grotesque Pro','FontWeight','Normal');
tmp_f18.CurrentAxes.Title.Position(2) = tmp_f18.CurrentAxes.Title.Position(2)+.4;

tmp_cc = 1;
for j = 1:length(sign_ROIs)
    tmp_ind10 = find(cell2mat(...
        DesikanKilliany_atlas_connectogram_labels(:,2))==sign_ROIs(j));
    tmp_alllabels(tmp_ind10).FontWeight = 'Bold';
    tmp_alllabels(tmp_ind10).FontSize = 12;
    if sum(dataMat_2{7,1}(tmp_ind10,:))>0 && ismember(j,tmp_hardcodedcoloringindices)
        tmp_alllabels(tmp_ind10).Color = tmp_colorpalette_ChordN(tmp_cc,:);
        BCC.setSquareN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:));
        %BCC.setSquareN(tmp_ind10,'EdgeColor',[.8 .2 .2],'LineWidth',2.5);
        BCC.setChordN(tmp_ind10,'FaceColor',tmp_colorpalette_ChordN(tmp_cc,:)); % sanity check
        %tmp_cc = tmp_cc+1;
    else
        tmp_alllabels(tmp_ind10).Color = '#a9a9a9';
        BCC.setSquareN(tmp_ind10,'FaceColor','#a9a9a9');
    end
    if sum(dataMat(tmp_ind10,:))>0
        tmp_cc = tmp_cc+1; % to keep the same coloring scheme as for dataMat
    end
end
clear tmp_cc tmp_ind10; 



%%% Widesearch with WDs

%%% recompute WDs with regressed out data (Y_reg2) which is age, gender,
%%% BDI, and STAI-T corrected (depending on the X definition earlier)

WD_HC_reg3 = cell(connectomes_HC_ind{1,2},1); % number of HC subjects
WD_FND_reg3{1,1} = cell(connectomes_FND_ind{1,2},1);
%WD_FND_reg3{2,1} = cell(connectomes_FND_ind{2,2},1);
for j = 1:length(WD_HC_reg3) % for each HC subject
    WD_HC_reg3{j,1} = zeros(size(connectomes_fa{1,1},1),1);
    tmp__connectome = jVecToUpperTriMat(Y_reg3(j,:),size(connectomes_fa{1,1},2));
    for n = 1:size(connectomes_fa{1,1},1) % for each node (i.e., region)
        WD_HC_reg3{j,1}(n,1) = sum([...
            tmp__connectome(1:(n-1),n)' ...
            tmp__connectome(n,(n+1):end)]); % because connectomes_* is not symmetrized
    end
end
for j = 1:length(WD_FND_reg3{1,1})
    WD_FND_reg3{1,1}{j,1} = zeros(size(connectomes_fa{1,1},1),1);
    tmp__connectome = jVecToUpperTriMat(Y_reg3(j+connectomes_HC_ind{1,2},:),size(connectomes_fa{1,1},2));
    for n = 1:size(connectomes_fa{1,1},1)
        WD_FND_reg3{1,1}{j,1}(n,1) = sum([...
            tmp__connectome(1:(n-1),n)' ...
            tmp__connectome(n,(n+1):end)]);
    end
end
clear tmp__connectome j n;

%WD_FND_mat = cell2mat(WD_FND{1,1}')'; % 85 FND patients x 84 ROIs (WDs)
WD_FND_mat = cell2mat(WD_FND_reg3{1,1}')'; % [regressed out version]
%WD_HC_mat  = cell2mat(WD_HC')';  % 75 HCs x 84 ROIs (WDs)
WD_HC_mat  = cell2mat(WD_HC_reg3')';
WD_combined_mat = [WD_HC_mat;WD_FND_mat]; % 160 x 84 ROIs (WDs)

widesearch_output_FND_corr_wd      = zeros(size(WD_FND_mat,2),size(widesearch_vars_FND,2)); % 84 ROIs x 12 clin vars
widesearch_output_FND_corr_wd_pval = zeros(size(WD_FND_mat,2),size(widesearch_vars_FND,2));
widesearch_output_HC_corr_wd       = zeros(size(WD_HC_mat,2),size(widesearch_vars_HC,2)); % 84 ROIs x 9 clin vars
widesearch_output_HC_corr_wd_pval  = zeros(size(WD_HC_mat,2),size(widesearch_vars_HC,2));
widesearch_output_combined_corr_wd       = zeros(size(WD_combined_mat,2),size(widesearch_vars_combined,2)); % 84 ROIs x 9 clin vars
widesearch_output_combined_corr_wd_pval  = zeros(size(WD_combined_mat,2),size(widesearch_vars_combined,2));

for q = 1:size(widesearch_output_FND_corr_wd,2)
    for k = 1:size(WD_FND_mat,2)
        [widesearch_output_FND_corr_wd(k,q),...
            widesearch_output_FND_corr_wd_pval(k,q)] = corr(WD_FND_mat(:,k),widesearch_vars_FND(:,q),'rows','pairwise');
    end
end
for q = 1:size(widesearch_output_HC_corr_wd,2)
    for k = 1:size(WD_HC_mat,2)
        [widesearch_output_HC_corr_wd(k,q),...
            widesearch_output_HC_corr_wd_pval(k,q)] = corr(WD_HC_mat(:,k),widesearch_vars_HC(:,q),'rows','pairwise');
    end
end
for q = 1:size(widesearch_output_combined_corr_wd,2)
    for k = 1:size(WD_combined_mat,2)
        [widesearch_output_combined_corr_wd(k,q),...
            widesearch_output_combined_corr_wd_pval(k,q)] = corr(WD_combined_mat(:,k),widesearch_vars_combined(:,q),'rows','pairwise');
    end
end

tmp_sorted_widesearch_output_FND_corr_wd_pval = widesearch_output_FND_corr_wd_pval;
tmp_sorted_widesearch_output_FND_corr_wd_indz = widesearch_output_FND_corr_wd_pval;
for q = 1:size(widesearch_output_FND_corr_wd_pval,2)
    [tmp_sorted_widesearch_output_FND_corr_wd_pval(:,q),...
        tmp_sorted_widesearch_output_FND_corr_wd_indz(:,q)] = sort(widesearch_output_FND_corr_wd_pval(:,q));
end
tmp_sorted_widesearch_output_HC_corr_wd_pval = widesearch_output_HC_corr_wd_pval;
tmp_sorted_widesearch_output_HC_corr_wd_indz = widesearch_output_HC_corr_wd_pval;
for q = 1:size(widesearch_output_HC_corr_wd_pval,2)
    [tmp_sorted_widesearch_output_HC_corr_wd_pval(:,q),...
        tmp_sorted_widesearch_output_HC_corr_wd_indz(:,q)] = sort(widesearch_output_HC_corr_wd_pval(:,q));
end
tmp_sorted_widesearch_output_combined_corr_wd_pval = widesearch_output_combined_corr_wd_pval;
tmp_sorted_widesearch_output_combined_corr_wd_indz = widesearch_output_combined_corr_wd_pval;
for q = 1:size(widesearch_output_combined_corr_wd_pval,2)
    [tmp_sorted_widesearch_output_combined_corr_wd_pval(:,q),...
        tmp_sorted_widesearch_output_combined_corr_wd_indz(:,q)] = sort(widesearch_output_combined_corr_wd_pval(:,q));
end

%%% add R2 computation for all cases (distribution is plotted then as well)
widesearch_output_FND_corr_wd_R2 = zeros(size(widesearch_output_FND_corr_wd_pval));
widesearch_output_HC_corr_wd_R2 = zeros(size(widesearch_output_HC_corr_wd_pval));
widesearch_output_combined_corr_wd_R2 = zeros(size(widesearch_output_combined_corr_wd_pval));
for q = 1:size(widesearch_output_FND_corr_wd_R2,2) % 12 clin vars
    for l = 1:size(WD_FND_mat,2)
        tmp_PP = polyfit(widesearch_vars_FND(:,q),WD_FND_mat(:,l),1);
        tmp_yfit = polyval(tmp_PP,widesearch_vars_FND(:,q));
        widesearch_output_FND_corr_wd_R2(l,q) = ...
            1 - sum((WD_FND_mat(:,l)-tmp_yfit).^2)/(...
        (length(WD_FND_mat(:,l))-1) * var(WD_FND_mat(:,l)));
    end
end
for q = 1:size(widesearch_output_HC_corr_wd_R2,2) % 9 clin vars
    for l = 1:size(WD_combined_mat,2)
        tmp_PP = polyfit(widesearch_vars_HC(:,q),WD_HC_mat(:,l),1);
        tmp_yfit = polyval(tmp_PP,widesearch_vars_HC(:,q));
        widesearch_output_HC_corr_wd_R2(l,q) = ...
            1 - sum((WD_HC_mat(:,l)-tmp_yfit).^2)/(...
        (length(WD_HC_mat(:,l))-1) * var(WD_HC_mat(:,l)));
        tmp_PP = polyfit(widesearch_vars_combined(:,q),WD_combined_mat(:,l),1);
        tmp_yfit = polyval(tmp_PP,widesearch_vars_combined(:,q));
        widesearch_output_combined_corr_wd_R2(l,q) = ...
            1 - sum((WD_combined_mat(:,l)-tmp_yfit).^2)/(...
        (length(WD_combined_mat(:,l))-1) * var(WD_combined_mat(:,l)));
    end
end

q_val = 0.05; % [if not defined already]
%q_val_strict = 0.01;
sign_links_nb_wd_all = zeros(size(widesearch_output_FND_corr_wd_pval,2),3); % store the # of sign. links surviving each fitting (HC, FND, HC+FND (combined))
for q = 1:size(widesearch_output_HC_corr_wd_pval,2)
    [~,tmp_pval_1] = fdr_bh(widesearch_output_HC_corr_wd_pval(:,q),q_val,'pdep','no');
    sign_links_nb_wd_all(q,1) = length(find(widesearch_output_HC_corr_wd_pval(:,q)<=tmp_pval_1)); % find # of links
    [~,tmp_pval_1] = fdr_bh(widesearch_output_combined_corr_wd_pval(:,q),q_val,'pdep','no');
    sign_links_nb_wd_all(q,3) = length(find(widesearch_output_combined_corr_wd_pval(:,q)<=tmp_pval_1)); % find # of links
end
for q = 1:size(widesearch_output_FND_corr_wd_pval,2)
    [~,tmp_pval_1] = fdr_bh(widesearch_output_FND_corr_wd_pval(:,q),q_val,'pdep','no');
    sign_links_nb_wd_all(q,2) = length(find(widesearch_output_FND_corr_wd_pval(:,q)<=tmp_pval_1)); % find # of links
end
clear q tmp_pval_1;



%% Find overlapping regions between signROIs (significant regions differing
%   between HCs and FND patients in terms of WD) and regions popping up
%   from the above analysis with the widesearch_output_FND_corr_wd_pval
%   variable, for clinical variables: S-FMDRS, CGI, and SF36-PhysFunc

signROIs_WD_FND_clinvars = zeros(84,3); % for the three clin. metrics
k = 1;
for q = [2 3 7]
    [~,tmp_pval_1] = fdr_bh(widesearch_output_FND_corr_wd_pval(:,q),q_val,'pdep','no'); % here we assess the FDR-corr. p-val (with q = 0.01) for the diff. clin. vars
    tmp_indz1 = find(widesearch_output_FND_corr_wd_pval(:,q)<=tmp_pval_1);
    %[~,tmp_pval_1] = fdr_bh(widesearch_output_combined_corr_wd_pval(:,q),q_val,'pdep','no');
    %tmp_indz1 = find(widesearch_output_combined_corr_wd_pval(:,q)<=tmp_pval_1);
    signROIs_WD_FND_clinvars(1:length(tmp_indz1),k) = tmp_indz1;
    k = k+1;
end
clear k tmp_indz1 tmp_pval_1;



%%
close all;
tmp_f20 = figure('Position',[-1250 100 960 520]);
tmp_clrmap = colormap(turbo(size(widesearch_output_FND_corr_wd_pval,2)));
plot(tmp_sorted_widesearch_output_FND_corr_wd_pval(:,1),'LineWidth',1.3,'Color',tmp_clrmap(1,:)); % DurationSymptoms (for FND)
hold(tmp_f20.CurrentAxes,'on');
for q = 2:size(widesearch_output_FND_corr_wd_pval,2)
    plot(tmp_sorted_widesearch_output_FND_corr_wd_pval(:,q),'LineWidth',1.3,'Color',tmp_clrmap(q,:));
    %hold(tmp_f20.CurrentAxes,'on');
end
tmp_f20.CurrentAxes.FontSize = 12;
tmp_f20.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-5 90 -.1 1.1]);
xticks(0:5:85); yticks(0:0.1:1); tmp_f20.CurrentAxes.YGrid = 'on';
plot([0 85],[.05 .05],'bla-','LineWidth',1.3);
text(78,.09,'\itp\rm < .05','FontSize',14,'FontName','Basis Grotesque Pro');
legend([' Illness Duration {' num2str(sign_links_nb_wd_all(1,2)) '}'],...
    [' S-FMDRS {' num2str(sign_links_nb_wd_all(2,2)) '}'],...
    [' CGI {' num2str(sign_links_nb_wd_all(3,2)) '}'],...
    [' SF36-PhysHealth {' num2str(sign_links_nb_wd_all(4,2)) '}'],...
    [' SF36-MentalHealth {' num2str(sign_links_nb_wd_all(5,2)) '}'],...
    [' SF36-GenHealth {' num2str(sign_links_nb_wd_all(6,2)) '}'],...
    [' SF36-PhysFunc {' num2str(sign_links_nb_wd_all(7,2)) '}'],[' AUCi [cortisol] {' num2str(sign_links_nb_wd_all(8,2)) '}'],...
    [' CAR_AUCi [cortisol] {' num2str(sign_links_nb_wd_all(9,2)) '}'],...
    [' BDI {' num2str(sign_links_nb_wd_all(10,2)) '}'],[' STAI-S {' num2str(sign_links_nb_wd_all(11,2)) '}'],...
    [' STAI-T {' num2str(sign_links_nb_wd_all(12,2)) '}'],...
    'FontSize',14,'FontName','Basis Grotesque Pro',...
    'Interpreter','None','Location','Best');
%text(1100,.26,'CGI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(3,:));
%text(1200,.165,'S-FMDRS','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(2,:));
%text(880,.2,'SF36-PhysFunc','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(7,:));
%text(600,.33,'BDI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(10,:));
%text(600,.41,'STAI-T','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(12,:));
tmp_f20.CurrentAxes.YAxis.Label.String = '\itp\rm-value';
tmp_f20.CurrentAxes.YAxis.Label.FontSize = 18;
tmp_f20.CurrentAxes.YAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f20.CurrentAxes.XAxis.Label.String = 'Sorted WD';
tmp_f20.CurrentAxes.XAxis.Label.FontSize = 18;
tmp_f20.CurrentAxes.XAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f20.CurrentAxes.Title.String = ...
    ['Sorted WD in all 84 ROIs in FND patients (q_{FDR} = ' num2str(q_val) ')'];
tmp_f20.CurrentAxes.Title.FontSize = 20;
tmp_f20.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f20.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';

tmp_f20_a = figure('Position',[-1250 100 960 520]);
%tmp_clrmap = colormap(turbo(size(widesearch_output_HC_corr_fa_pval,2)));
plot(tmp_sorted_widesearch_output_HC_corr_wd_pval(:,1),'LineWidth',1.3,'Color',tmp_clrmap(4,:)); % "SF36-PhysHealth" (for combined)
hold(tmp_f20_a.CurrentAxes,'on');
for q = 2:size(widesearch_output_HC_corr_wd_pval,2)
    plot(tmp_sorted_widesearch_output_HC_corr_wd_pval(:,q),'LineWidth',1.3,'Color',tmp_clrmap(q+3,:));
    %hold(tmp_f20_a.CurrentAxes,'on');
end
tmp_f20_a.CurrentAxes.FontSize = 12;
tmp_f20_a.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-5 90 -0.1 1.1]);
xticks(0:5:85); yticks(0:.1:1); tmp_f20_a.CurrentAxes.YGrid = 'on';
plot([0 85],[.05 .05],'bla-','LineWidth',1.3);
text(78,.09,'\itp\rm < .05','FontSize',14,'FontName','Basis Grotesque Pro');
legend([' SF36-PhysHealth {' num2str(sign_links_nb_wd_all(1,1)) '}'],...
    [' SF36-MentalHealth {' num2str(sign_links_nb_wd_all(2,1)) '}'],...
    [' SF36-GenHealth {' num2str(sign_links_nb_wd_all(3,1)) '}'],...
    [' SF36-PhysFunc {' num2str(sign_links_nb_wd_all(4,1)) '}'],...
    [' AUCi [cortisol] {' num2str(sign_links_nb_wd_all(5,1)) '}'],...
    [' CAR_AUCi [cortisol] {' num2str(sign_links_nb_wd_all(6,1)) '}'],...
    [' BDI {' num2str(sign_links_nb_wd_all(7,1)) '}'],...
    [' STAI-S {' num2str(sign_links_nb_wd_all(8,1)) '}'],[' STAI-T {' num2str(sign_links_nb_wd_all(9,1)) '}'],...
    'FontSize',14,'FontName','Basis Grotesque Pro',...
    'Interpreter','None','Location','Best');
%text(1100,.26,'CGI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(3,:));
%text(1370,.4,'SF36-GenHealth','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(6,:));
%text(1370,.2,'SF36-PhysFunc','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(7,:));
%text(1390,.5,'BDI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(10,:));
%text(600,.41,'STAI-T','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(12,:));
tmp_f20_a.CurrentAxes.YAxis.Label.String = '\itp\rm-value';
tmp_f20_a.CurrentAxes.YAxis.Label.FontSize = 18;
tmp_f20_a.CurrentAxes.YAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f20_a.CurrentAxes.XAxis.Label.String = 'Sorted WD';
tmp_f20_a.CurrentAxes.XAxis.Label.FontSize = 18;
tmp_f20_a.CurrentAxes.XAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f20_a.CurrentAxes.Title.String = ...
    ['Sorted WDs in all 84 ROIs in HCs (q_{FDR} = ' num2str(q_val) ')'];
tmp_f20_a.CurrentAxes.Title.FontSize = 20;
tmp_f20_a.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f20_a.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';

tmp_f20_b = figure('Position',[-1250 100 960 520]);
%tmp_clrmap = colormap(turbo(size(widesearch_output_HC_corr_fa_pval,2)));
plot(tmp_sorted_widesearch_output_combined_corr_wd_pval(:,1),'LineWidth',1.3,'Color',tmp_clrmap(4,:)); % "SF36-PhysHealth" (for combined)
hold(tmp_f20_b.CurrentAxes,'on');
for q = 2:size(widesearch_output_combined_corr_wd_pval,2)
    plot(tmp_sorted_widesearch_output_combined_corr_wd_pval(:,q),'LineWidth',1.3,'Color',tmp_clrmap(q+3,:));
    %hold(tmp_f20_b.CurrentAxes,'on');
end
tmp_f20_b.CurrentAxes.FontSize = 12;
tmp_f20_b.CurrentAxes.FontName = 'Basis Grotesque Pro';
axis([-5 90 -0.1 1.1]);
xticks(0:5:85); yticks(0:.1:1); tmp_f20_b.CurrentAxes.YGrid = 'on';
plot([0 85],[.05 .05],'bla-','LineWidth',1.3);
text(78,.09,'\itp\rm < .05','FontSize',14,'FontName','Basis Grotesque Pro');
legend([' SF36-PhysHealth {' num2str(sign_links_nb_wd_all(1,3)) '}'],...
    [' SF36-MentalHealth {' num2str(sign_links_nb_wd_all(2,3)) '}'],...
    [' SF36-GenHealth {' num2str(sign_links_nb_wd_all(3,3)) '}'], ...
    [' SF36-PhysFunc {' num2str(sign_links_nb_wd_all(4,3)) '}'], ...
    [' AUCi [cortisol] {' num2str(sign_links_nb_wd_all(5,3)) '}'],...
    [' CAR_AUCi [cortisol] {' num2str(sign_links_nb_wd_all(6,3)) '}'], ...
    [' BDI {' num2str(sign_links_nb_wd_all(7,3)) '}'], ...
    [' STAI-S {' num2str(sign_links_nb_wd_all(8,3)) '}'], ...
    [' STAI-T {' num2str(sign_links_nb_wd_all(9,3)) '}'],...
    'FontSize',14,'FontName','Basis Grotesque Pro',...
    'Interpreter','None','Location','Best');
%text(1100,.26,'CGI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(3,:));
%text(1370,.4,'SF36-GenHealth','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(6,:));
%text(1370,.2,'SF36-PhysFunc','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(7,:));
%text(1390,.5,'BDI','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(10,:));
%text(600,.41,'STAI-T','FontSize',12,'FontName','Basis Grotesque Pro','Color',tmp_clrmap(12,:));
tmp_f20_b.CurrentAxes.YAxis.Label.String = '\itp\rm-value';
tmp_f20_b.CurrentAxes.YAxis.Label.FontSize = 18;
tmp_f20_b.CurrentAxes.YAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f20_b.CurrentAxes.XAxis.Label.String = 'Sorted WD';
tmp_f20_b.CurrentAxes.XAxis.Label.FontSize = 18;
tmp_f20_b.CurrentAxes.XAxis.Label.FontName = 'Basis Grotesque Pro';
tmp_f20_b.CurrentAxes.Title.String = ...
    ['Sorted WDs in all ROIs in HC and FND patients (q_{FDR} = ' num2str(q_val) ')'];
tmp_f20_b.CurrentAxes.Title.FontSize = 20;
tmp_f20_b.CurrentAxes.Title.FontWeight = 'Normal';
tmp_f20_b.CurrentAxes.Title.FontName = 'Basis Grotesque Pro';

%
% exportgraphics(tmp_f20, 'Figures/Figure20_q001_forInkScape2.png', 'Resolution', '300');
% exportgraphics(tmp_f20_a, 'Figures/Figure20a_q001_forInkScape2.png', 'Resolution', '300');
% exportgraphics(tmp_f20_b, 'Figures/Figure20b_q001.png', 'Resolution', '300');

% exportgraphics(tmp_f20, 'Figures/Figure20_q005_age_gender_BDI_STAI-T_reg.png', 'Resolution', '300');
% exportgraphics(tmp_f20_a, 'Figures/Figure20a_q005_age_gender_BDI_STAI-T_reg.png', 'Resolution', '300');
% exportgraphics(tmp_f20_b, 'Figures/Figure20b_q005_age_gender_BDI_STAI-T_reg.png', 'Resolution', '300');



%% Explore the 21 ROIs (significant ones from WD first analysis)
%   for differences in subgroups (mFND, PNES, PPPD, etc.)

boxplots_PPPD = zeros(length(patients_inds_PPPD_withinFND),length(sign_ROIs));
boxplots_PNES = zeros(length(patients_inds_PNES_withinFND),length(sign_ROIs));
boxplots_mFND = zeros(length(patients_inds_mFND_withinFND),length(sign_ROIs));
boxplots_HCs  = zeros(length(connectomes_HC_ind{1,2}),length(sign_ROIs));
for j = 1:length(sign_ROIs)
    for i = 1:length(patients_inds_PPPD_withinFND)
        boxplots_PPPD(i,j) = WD_FND{1,1}{patients_inds_PPPD_withinFND(i),1}(j,1)/length(WD_FND{1,1}{1,1});
    end
    for i = 1:length(patients_inds_PNES_withinFND)
        boxplots_PNES(i,j) = WD_FND{1,1}{patients_inds_PNES_withinFND(i),1}(j,1)/length(WD_FND{1,1}{1,1});
    end
    for i = 1:length(patients_inds_mFND_withinFND)
        boxplots_mFND(i,j) = WD_FND{1,1}{patients_inds_mFND_withinFND(i),1}(j,1)/length(WD_FND{1,1}{1,1});
    end
    for i = 1:length(WD_HC)
        boxplots_HCs(i,j) = WD_HC{i,1}(j,1)/length(WD_HC{1,1});
    end
end

%
tmp_f14 = figure('Position',[1729 38 1920 976]);
for j = 1:length(sign_ROIs)
    tmp_axx = subplot(3,7,j);
    tmp_mrksize = 5;
    swarmchart(1,boxplots_PPPD(:,j),tmp_mrksize); hold(tmp_axx,'on');
    swarmchart(2,boxplots_PNES(:,j),tmp_mrksize);
    swarmchart(3,boxplots_mFND(:,j),tmp_mrksize);
    swarmchart(4,boxplots_HCs(:,j),tmp_mrksize);
    axis([0.5 4.5 0 .75]);
    xticks(1:4);
    xticklabels({'PPPD','PNES','mFND','HC'}); tmp_axx.FontSize = 9;
    tmp_axx.FontName = 'Basis Grotesque Arabic Pro';
    if mod(j,7)==1, ylabel('Mean FA','FontSize',12); end
    if j>14, xlabel('Subgroups','FontSize',12); end     
    title(DesikanKilliany_atlas_sorted{sign_ROIs(j),3},'FontSize',14,...
        'FontName','Basis Grotesque Arabic Pro','FontWeight','Normal');
    box('on');
    [~,ptmp_] = ttest2(boxplots_PPPD(:,j),boxplots_HCs(:,j),'tail','left');
    if ptmp_<=(0.05/length(sign_ROIs))
        fprintf(['\nSign. diff. for j=' num2str(j) ' (p-val = ' num2str(ptmp_) ').\n']);
        plot([1 4],[.65 .65],'black'); text((4-1)/2+1,.7,'*','FontSize',16);
    end
    [~,ptmp_] = ttest2(boxplots_PNES(:,j),boxplots_HCs(:,j),'tail','left');
    if ptmp_<=(0.05/length(sign_ROIs))
        fprintf(['\nSign. diff. for j=' num2str(j) ' (p-val = ' num2str(ptmp_) ').\n']);
        plot([2 4],[.65 .65],'black'); text((4-2)/2+2,.7,'*','FontSize',16);
    end
    [~,ptmp_] = ttest2(boxplots_mFND(:,j),boxplots_HCs(:,j),'tail','left');
    if ptmp_<=(0.05/length(sign_ROIs))
        fprintf(['\nSign. diff. for j=' num2str(j) ' (p-val = ' num2str(ptmp_) ').\n']);
        plot([3 4],[.65 .65],'black'); text((4-3)/2+3,.7,'*','FontSize',16);
    end
    hold(tmp_axx,'off');
end
if ~exist([export_figs_path 'Figure14_test.pdf'],'file')
    %export_fig([export_figs_path 'Figure13_tmp.pdf'],tmp_f3);
    %exportgraphics(tmp_f14,[export_figs_path 'Figure14_test.pdf'],...
    %    'BackgroundColor','none','ContentType','vector');
end
clear tmp_mrksize;



%% Within-subgroups analyses for T1 vs FUP data
%   Looking at PPPD, PNES, mFND

patients_inds_FUP = find(~isnan(patients_datatable.Improved));
patients_inds_FND = tmp_FNDindzonly;

patients_inds_PPPD_FUP_withinFND_T1 = ...
    patients_inds_PPPD(ismember(patients_inds_PPPD_withinFND,patients_inds_impr_withinFND) | ...
    ismember(patients_inds_PPPD_withinFND,patients_inds_nimpr_withinFND));
patients_inds_PPPD_FUP_withinFND_FUP = find(ismember(patients_inds_FUP,patients_inds_PPPD_FUP_withinFND_T1));
patients_inds_PPPD_FUP_withinFND_T1 = find(ismember(patients_inds_FND,patients_inds_PPPD_FUP_withinFND_T1));

patients_inds_PNES_FUP_withinFND_T1 = ...
    patients_inds_PNES(ismember(patients_inds_PNES_withinFND,patients_inds_impr_withinFND) | ...
    ismember(patients_inds_PNES_withinFND,patients_inds_nimpr_withinFND));
patients_inds_PNES_FUP_withinFND_FUP = find(ismember(patients_inds_FUP,patients_inds_PNES_FUP_withinFND_T1));
patients_inds_PNES_FUP_withinFND_T1 = find(ismember(patients_inds_FND,patients_inds_PNES_FUP_withinFND_T1));

patients_inds_mFND_FUP_withinFND_T1 = ...
    patients_inds_mFND(ismember(patients_inds_mFND_withinFND,patients_inds_impr_withinFND) | ...
    ismember(patients_inds_mFND_withinFND,patients_inds_nimpr_withinFND));
patients_inds_mFND_FUP_withinFND_FUP = find(ismember(patients_inds_FUP,patients_inds_mFND_FUP_withinFND_T1));
patients_inds_mFND_FUP_withinFND_T1 = find(ismember(patients_inds_FND,patients_inds_mFND_FUP_withinFND_T1));

%   Differences within the PPPD subgroup between T1 and FUP,
%   in terms of microstructural integrity?

WD_FND_T1_mat  = reshape(cell2mat(WD_FND{1,1}),length(WD_FND{1,1}{1,1}),length(WD_FND{1,1})); % per column for each subject
WD_FND_FUP_mat = reshape(cell2mat(WD_FND{2,1}),length(WD_FND{2,1}{1,1}),length(WD_FND{2,1}));
tmp_diffz = zeros(84,2);
for j = 1:length(tmp_diffz)
    [tmp_diffz(j,1), tmp_diffz(j,2)] = ttest2(...
        WD_FND_T1_mat(j,patients_inds_PPPD_FUP_withinFND_T1),...
        WD_FND_FUP_mat(j,patients_inds_PPPD_FUP_withinFND_FUP),'tail','right');
    if tmp_diffz(j,1)
        fprintf(['\nFound a sign. change (btw T1 and FUP) for PPPD (j = ' ...
            num2str(j) '), p (unc.) = ' num2str(tmp_diffz(j,2)) '.\n']);
    end
end

%   Same for PNES and mFND?
tmp_diffz = zeros(84,2);
for j = 1:length(tmp_diffz)
    [tmp_diffz(j,1), tmp_diffz(j,2)] = ttest2(...
        WD_FND_T1_mat(j,patients_inds_PNES_FUP_withinFND_T1),...
        WD_FND_FUP_mat(j,patients_inds_PNES_FUP_withinFND_FUP),'tail','right');
    if tmp_diffz(j,1)
        fprintf(['\nFound a sign. change (btw T1 and FUP) for PNES (j = ' ...
            num2str(j) '), p (unc.) = ' num2str(tmp_diffz(j,2)) '.\n']);
    end
end

tmp_diffz = zeros(84,2);
for j = 1:length(tmp_diffz)
    [tmp_diffz(j,1), tmp_diffz(j,2)] = ttest2(...
        WD_FND_T1_mat(j,patients_inds_mFND_FUP_withinFND_T1),...
        WD_FND_FUP_mat(j,patients_inds_mFND_FUP_withinFND_FUP),'tail','right');
    if tmp_diffz(j,1)
        fprintf(['\nFound a sign. change (btw T1 and FUP) for mFND (j = ' ...
            num2str(j) '), p (unc.) = ' num2str(tmp_diffz(j,2)) '.\n']);
    end
end

tmp_diffz = zeros(84,2);
for j = 1:length(tmp_diffz)
    [tmp_diffz(j,1), tmp_diffz(j,2)] = ttest(...
        WD_FND_T1_mat(j,ismember(patients_inds_FND,patients_inds_FUP)),...
        WD_FND_FUP_mat(j,:),'tail','left');
    if tmp_diffz(j,1)
        fprintf(['\nFound a sign. change (btw T1 and FUP) for all FND taken together (j = ' ...
            num2str(j) '), p (unc.) = ' num2str(tmp_diffz(j,2)) '.\n']);
    end
end



%% Misc

function [avg_intra, avg_inter] = compute_hemispheric_metric(M_in,flag_abs)
    tmp_midbrain_th = size(M_in,1)/2; % should be dividible though (usually most atlases are)
    M_in_up = M_in(1:tmp_midbrain_th,1:tmp_midbrain_th); % split matrix in two subparts
    M_in_lo = M_in(tmp_midbrain_th+1:end,tmp_midbrain_th+1:end);
    if ~flag_abs
        avg_intra = mean([jUpperTriMatToVec(triu(M_in_up,1));jUpperTriMatToVec(triu(M_in_lo,1))],'omitnan');
        avg_inter = mean(mean(M_in(1:tmp_midbrain_th,(tmp_midbrain_th+1):end),'omitnan'),'omitnan');
    else
        avg_intra = mean(abs([jUpperTriMatToVec(triu(M_in_up,1));jUpperTriMatToVec(triu(M_in_lo,1))]),'omitnan');
        avg_inter = mean(mean(abs(M_in(1:tmp_midbrain_th,(tmp_midbrain_th+1):end)),'omitnan'),'omitnan');
    end
end


