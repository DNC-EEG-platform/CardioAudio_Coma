%% Description
% The following script performs the cluster permutation statistical analysis
% for comparing condition-specific EEG data as part of the cardio-audio
% experiment in comatose patients

% Requires Fieldtrip toolbox (https://www.fieldtriptoolbox.org/)

% called by Run_ClusterAnalyses.m script to perform statistical analysis step
%
% %
%     Copyright (C) 2024, Andria Pelentritou
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% %
%
%% Code
function stat = run_clustperm(data,time,fs,channames,elc_file,stat,compute_stat)
%input data = 720tf x 62chan x 2cond x 31subj %600 ms trials x 62 EEG channels x 2 conditions x 31 subjects


[ntf,nchan,ncond,nsubj] = size(data);

% make electrode structure, needs the standard .elc file
elec = ft_read_sens(elc_file);
idx = nan(size(channames));
for i=1:length(idx)
    tmp = find(strcmpi(channames{i},elec.label));
    if ~isempty(tmp)
        idx(i) = tmp;
    end
end

if sum(isnan(idx))==0
    elec.chanpos=elec.chanpos(idx,:);
    elec.elecpos=elec.elecpos(idx,:);
    elec.label=elec.label(idx,:);
else
    sprintf('Error! 3D electrode positions missing for %d channels.',sum(isnan(idx)))
    return
end

% make a template ft structure
data_avg              = [];
data_avg.fsample      = fs; %chan x time
data_avg.avg          = nan(nchan,ntf);% 1 x tf
data_avg.time         = time; %1 x time
data_avg.label        = channames;%chan x 1;
data_avg.dimord       = 'chan_time';
data_avg.elec         = elec;
data_avg.cfg          = [];

% make cell array with all data in separate ft structures
allsubj = cell(nsubj,ncond);
allsubj(:,:) = {data_avg};
clear data_avg chanlocs;

for cond = 1:ncond
    for subj = 1:nsubj    
        allsubj{subj,cond}.avg = squeeze(data(:,:,cond,subj))'; %chan x time  
    end
end

%following steps from fieldtrip tutorial sit: http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock

% prepare neighbours for significant cluster definition
cfg_neighb        = [];
cfg_neighb.method = 'distance';  
neighbours        = ft_prepare_neighbours(cfg_neighb, allsubj{1,1});

% Cluster permutation statistics 
if compute_stat == 1
    %cluster permutation parameters                           
    cfg = [];
    cfg.channel         = 'all';                      
    cfg.latency         = 'all';
    cfg.method          = 'montecarlo';
    cfg.statistic       = 'depsamplesT';
    cfg.correctm        = 'cluster';
    cfg.clusteralpha    = 0.05;
    cfg.alpha           = 0.05;
    cfg.avgoverchan     = 'no';
    cfg.avgovertime     = 'no';
    cfg.correcttail     = 'prob';
    cfg.numrandomization = 5000;
    cfg.uvar            = 1;
    cfg.ivar            = 2;
    cfg.neighbours      = neighbours;  
    
    %statistical design
    design = zeros(2,2*nsubj);
    for i = 1:nsubj
        design(1,i) = i;
    end
    for i = 1:nsubj
    design(1,nsubj+i) = i;
    end
    design(2,1:nsubj)        = 1;
    design(2,nsubj+1:2*nsubj) = 2;

    cfg.design = design;

    [stat] = ft_timelockstatistics(cfg, allsubj{:,1}, allsubj{:,2}); %cluster permutation fieldtrip command
end
end