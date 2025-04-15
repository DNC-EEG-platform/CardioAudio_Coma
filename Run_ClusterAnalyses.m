%% Description
% The following script performs the cluster permutation statistical analysis
% for comparing condition-specific EEG data as part of the cardio-audio
% experiment in comatose patients
%
% %
% Requires modification of data paths according to user
%
% %
% Dependencies 
% Fieldtrip toolbox (https://www.fieldtriptoolbox.org/download)
% run_clustperm.m script to perform statistical analysis step
% 'standard_1020.elc' electrode layout file from fieldtrip
%
% % 
% Inputs 
%
% dataall = 4-D matrix as timeframesxchannelxconditionlabelxNsubjects
% contains subject-wise grand average ERPs for conditions of interest over
% the trial of interest for each channel in the EEG layout
% e.g. 720 timeframes x 65 channels x 2 condition labels x 31 subjects
% for condition labels in our case see below
%
% chanlabel = channelx1 cell
% each row contains the channel name in string
% dataall = 4-D matrix as timeframesxchannelsxconditionlabelxNsubjects
% contains subject-wise grand average ERPs for conditions of interest
% contrasted in the cluster permutation test
% e.g. 720 timeframes x 65 channels x 2 condition labels x 31 subjects
% for condition labels in our case see below
%
% loaded in this script as e.g. '../groupN31.mat'
%
%
% Outputs
%
% data2 = 4-D matrix
% stat = structure of the output of fieldtrip's ft_timelockstatistics
% containing the cluster permutation statistical analysis results
% for details visit:
% http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock
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
clear all;
OSflag=1; %1 = Linux; 2 = Windows

%clear out fieldtrip dependencies
restoredefaultpath

if OSflag==1
    addpath /mnt/HDD1/CardioAudio_coma/scripts/fieldtrip-20201205/;
    addpath /mnt/HDD1/CardioAudio_coma/scripts/ScriptsLinux/;
    addpath /mnt/HDD1/CardioAudio_coma/data/;
    
    mpath=['/mnt/HDD1/CardioAudio_coma/data/'];
    
    pt          = [];
    pt.root     = mpath; 
    pt.data     = [pt.root];
    pt.group    = [pt.data,'group/'      ];
elseif OSflag==2
    addpath D:\CardioAudio_coma\scripts\fieldtrip-20201205\;
    addpath D:\CardioAudio_coma\scripts\Scripts\;
    addpath D:\CardioAudio_coma\data\;
    
    mpath=['D:\CardioAudio_coma\data\'];
    
    pt          = [];
    pt.root     = mpath; 
    pt.data     = [pt.root];
    pt.group    = [pt.data,'group\'      ];
end

%restore fieldtrip defaults
ft_defaults

%valid subjects
subjs=[1:31]; %1:31 for FO included subjects; 1:17 for UO included subjects

%select conditions to contrast
%   CONDITION LABELS
%         first number
%             %1=synch
%             %2=asynch
%             %3=isoch
%         second number
%             %1=rpeaks during sound stimulation
%             %2=rpeaks during omission
%             %3=sound onset
%             %4=omission onset based on mean ISI

condoi=[12]; %deviant stimulus data label 
condoi_er=[22]; %standard stimulus data label
%the topoplot refers to condoi-condoi_er

clust_compute = 1; %run cluster permutation analysis in fieldtrip

%% ========================================================================            
%% ###11 cluster permutation test  ==================================
%% ========================================================================
% if find(cnt.steps==11)

    cd(pt.group);
    
    % load and organize data from event-related blocks
    filename = sprintf('groupN%d.mat',length(subjs));
    load(filename);
   
    alldata_er=alldata;
    clear alldata filename
    
    chanind = 1:62; %EEG channels to include
    condind = zeros(size(alldata_er,2));
    for i = 1:length(condoi_er)
        condind(conds==condoi_er(i))=1;
    end
    condind = logical(condind);
    
    channames = chanlabel(chanind); %EEG channel labels
    fs = 1200; %sampling rate
    elc_file = 'standard_1020.elc'; %electrode template from fieldtrip
    
    data_er = alldata_er(:,chanind,condind,:); %standard stimulus data
    
    % load and organize data from second condition
    filename = sprintf('groupN%d.mat',length(subjs));
    load(filename)
    
    chanind = 1:62; %EEG channels to include
    condind = zeros(size(alldata,2));
    for i = 1:length(condoi)
        condind(conds==condoi(i))=1;
    end
    condind = logical(condind);
    
    channames = chanlabel(chanind); %EEG channel labels
    fs = 1200; %sampling rate
    elc_file = 'standard_1020.elc'; %electrode template from fieldtrip
    
    data = alldata(:,chanind,condind,:); %deviant stimulus data
    
    data2=cat(3,data,data_er); %compile two condition data

     filename = sprintf('fstats_cond=%dvs%d_clust_n=%d.mat',condoi_er(1),condoi(1),length(subjs));
    
    %if computation wanted, use the 
        stat = [];
        clust_compute=1;
        stat = run_clustperm(data2,time,fs,channames,elc_file,stat,clust_compute); %cluster permutation function

       save(filename,'stat','data2');
    clear fs data chanind condind channames elc_file data filename stat; 
    
'done.'