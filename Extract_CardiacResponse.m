%% Description
%
% The following script extracts the RR intervals prior, during and after
% sound omissions to identify cardiac omission responses from ECG signals
% for the cardio audio experiment in comatose patients
%
% %
% Requires modification of data paths according to user
%
% %
% Inputs
%
% FOR EACH EXPERIMENTAL CONDITION AND PARTICIPANT
% T variable is a table array of dimension Yx3, where Y is the number of
% trials per experimental block and 3 columns as follows:
% T.trial = trial number per block
% T.cond = binary variable (1 for sound trial; 2 for omission trial)
% T.lsound = actual sound onset in samples
% loaded in this script as e.g. '../trg_S1.mat'
% example
% trial     cond    lsound
% 1         1       9173
% 2         1       10064
% 3         1       10984
% 4         1       11945
% 5         2       NaN
%
% r variable (Nx1 vector) containing all offline detected time stamps of
% R peaks (from raw ECG data timeseries)
% loaded in this script as e.g. '../ecg_S1.mat'
%
% %
% Output
%
% comaFO structure containing RR intervals for omission trials as follows:
% comaFO.RR_omis_m1 = average RR interval for trial prior to omission
% comaFO.RR_omis = average RR interval for omission trial
% comaFO.RR_omis_p1 = average RR interval for trial after omission
% comaFO.RR_omis_p2 = average RR interval for two trials after omission
% comaFO.RR_omis_XX variable (3xN matrix) with rows for each condition (synch,
% asynch, isoch) and columns for N subjects

% comaFO.RRall_omis_m1 = all RR intervals for trial prior to omission
% comaFO.RRall_omis = all RR intervals for omission trial
% comaFO.RRall_omis_p1 = all RR intervals for trial after omission
% comaFO.RRall_omis_p2 = all RR intervals for two trials after omission
% comaFO.RRall_omis_XX cell structure (3xN) with rows for each condition (synch,
% asynch, isoch) and columns for N sujects
% each cell is a Yx1 vector containing all available RR interval values
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
%set path
if OSflag==1
    mpath=['/mnt/HDD1/CardioAudio_coma/data/'];
    mpath_analysis=['/mnt/HDD1/CardioAudio_coma/analysis/CR/'];
    load '/mnt/HDD1/CardioAudio_coma/data/subj_list.mat' %variable containing list of valid fo and uo subjects
elseif OSflag==2
    mpath=['D:\CardioAudio_coma\data\'];
    mpath_analysis=['D:\CardioAudio_coma\analysis\CR\'];
    load 'D:\CardioAudio_coma\data\data\subj_list.mat' %variable containing list of valid fo and uo subjects
end

rng(1,'twister'); %set randomization

%valid subjects
sub_type=1; %1 for FO subjects; 2 for UO subjects
if sub_type==1
    subjs=fo_subjs; %list of FO subjects included in the study
elseif sub_type==2
    subjs=uo_subjs; %list of UO subjects included in the study
end
%valid conditions
condition{1}='S'; %synch
condition{2}='A'; %asynch
condition{3}='I'; %isoch
%sampling rate
Fs=1200;
%block number
block_num=[1:2]; %included blocks

compute_means=1; %1 for yes; 0 for no
exclude_outliers=1; %1 for yes; 0 for no

%% Check for changes across subjects
if compute_means==1
    for k=1:length(condition)
        for j=1:length(subjs)
            Tomis=[]; %R to R interval during omission
            Tomis_m1=[]; %R to R interval before omission
            Tomis_p1=[]; %R to R interval after omission
            Tomis_p2=[]; %2 R to R interval after omission

            for block=block_num

                if OSflag==1
                    load([mpath,'s',num2str(subjs(j)),'/process/trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
                    load([mpath,'s',num2str(subjs(j)),'/process/ecg_',condition{k},num2str(block),'.mat']) % ecg info before artefact rejection
                elseif OSflag==2
                    load([mpath,'s',num2str(subjs(j)),'\process\trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
                    load([mpath,'s',num2str(subjs(j)),'\process\ecg_',condition{k},num2str(block),'.mat']) % ecg info before artefact rejection
                end

                %prepare T variable
                s = T.lsound(T.cond==1);
                ss = diff(T.lsound);
                ss_avg = nanmean(ss);

                T.lomiss = nan(size(T.lsound));
                ind = find(T.cond==2);
                ind=ind(ind>1);
                T.lomiss(ind) = T.lsound(ind-1)+ss_avg;
                clear ind s ss ss_avg

                %remove SOSO trials (i.e. we do not consider trials where the omissions occur twice in a series of four trials such as sound, omission, sound, omission)
                T_temp = T;
                is = [];
                io=find(T.cond(2:end-3) == 2);
                io=io+1;
                for i=1:length(io)
                    if T.cond(io(i,1))==2 && T.cond((io(i,1))+2)==2
                        is = [is;io(i,1)-1;io(i,1);io(i,1)+1;io(i,1)+2];
                    end
                end

                T_temp.lomiss(is,:) = NaN;
                T = T_temp;
                clear T_temp io is;

                %start at 3 trial
                io=find(isnan(T.lomiss(3:end-3)) == 0);
                io=io+2;

                %identify omission R peaks
                for i=1:length(io)
                    valom = round(T.lomiss(io(i,1)));
                    [~,idx]=min(abs(r-valom));
                    minval = r(idx);
                    if minval > valom
                        Tomis = [Tomis; r(idx) - r(idx-1)];
                        Tomis_m1=[Tomis_m1; r(idx-1) - r(idx-2)];
                        Tomis_p1=[Tomis_p1; r(idx+1) - r(idx)];
                        Tomis_p2=[Tomis_p2; r(idx+2) - r(idx+1)];
                    else
                        %elseif minval < valom
                        Tomis = [Tomis; r(idx+1) - r(idx)];
                        Tomis_m1=[Tomis_m1; r(idx) - r(idx-1)];
                        Tomis_p1=[Tomis_p1; r(idx+2) - r(idx+1)];
                        Tomis_p2=[Tomis_p2; r(idx+3) - r(idx+2)];
                    end
                end
                clear io r T
            end

            %identify and remove outliers
            if exclude_outliers==1
                Tomisplot=Tomis;
                Tomisplot(find(isnan(Tomis) == 1))=[];
                [TomisClean,II1]=rmoutliers(Tomisplot);
                Tomism1plot=Tomis_m1;
                Tomism1plot(find(isnan(Tomis_m1) == 1))=[];
                [Tomism1Clean,II2]=rmoutliers(Tomism1plot);
                Tomisp1plot=Tomis_p1;
                Tomisp1plot(find(isnan(Tomis_p1) == 1))=[];
                [Tomisp1Clean,II3]=rmoutliers(Tomisp1plot);
                Tomisp2plot=Tomis_p2;
                Tomisp2plot(find(isnan(Tomis_p2) == 1))=[];
                [Tomisp2Clean,II4]=rmoutliers(Tomisp2plot);
            end


            if exclude_outliers==1
                %ensure equal numbers of distance to omission trials following outlier removal
                min_omistrl = min([length(TomisClean); length(Tomism1Clean); length(Tomisp1Clean); length(Tomisp2Clean)]);
                final_trl=randperm(length(TomisClean), min_omistrl);
                %Store subject-wise outlier free data
                RRall_omis{k,j}=TomisClean(final_trl,1);
                RRall_omis_m1{k,j}=Tomism1Clean(final_trl,1);
                RRall_omis_p1{k,j}=Tomisp1Clean(final_trl,1);
                RRall_omis_p2{k,j}=Tomisp2Clean(final_trl,1);

                %store group outlier free data
                RR_omis(k,j)=mean(RRall_omis{k,j},1);
                RR_omis_m1(k,j)=mean(RRall_omis_m1{k,j},1);
                RR_omis_p1(k,j)=mean(RRall_omis_p1{k,j},1);
                RR_omis_p2(k,j)=mean(RRall_omis_p2{k,j},1);
                clear TomisClean Tomism1Clean Tomisp1Clean Tomisp2Clean

            else %don't exclude outliers
                %ensure equal numbers of distance to omission trials following outlier removal
                min_omistrl = min([length(Tomis); length(Tomis_m1); length(Tomis_p1); length(Tomis_p2)]);
                final_trl=randperm(length(Tomis), min_omistrl);
                %Store subject-wise data
                RRall_omis{k,j}=Tomis(final_trl,1);
                RRall_omis_m1{k,j}=Tomis_m1(final_trl,1);
                RRall_omis_p1{k,j}=Tomis_p1(final_trl,1);
                RRall_omis_p2{k,j}=Tomis_p2(final_trl,1);

                %store group data
                RR_omis(k,j)=mean(RRall_omis{k,j},1);
                RR_omis_m1(k,j)=mean(RRall_omis_m1{k,j},1);
                RR_omis_p1(k,j)=mean(RRall_omis_p1{k,j},1);
                RR_omis_p2(k,j)=mean(RRall_omis_p2{k,j},1);
            end
            clear min_omistrl final_trl
        end
    end

    cd(mpath_analysis) %output path for saving
    if sub_type==1 %FO subjects
        %store relevant data for statistical analysis
        comaFO = [];
        comaFO.RR_omis_m1 = RR_omis_m1;
        comaFO.RR_omis = RR_omis;
        comaFO.RR_omis_p1 = RR_omis_p1;
        comaFO.RR_omis_p2 = RR_omis_p2;
        comaFO.RRall_omis_m1 = RRall_omis_m1;
        comaFO.RRall_omis = RRall_omis;
        comaFO.RRall_omis_p1 = RRall_omis_p1;
        comaFO.RRall_omis_p2 = RRall_omis_p2;
        if exclude_outliers==1
            save comaFO_CR_all comaFO
        else
            save comaFO_CR_withoutliers_all comaFO
        end
    elseif sub_type==2 %UO subjects
        %store relevant data for statistical analysis
        comaUO = [];
        comaUO.RR_omis_m1 = RR_omis_m1;
        comaUO.RR_omis = RR_omis;
        comaUO.RR_omis_p1 = RR_omis_p1;
        comaUO.RR_omis_p2 = RR_omis_p2;
        comaUO.RRall_omis_m1 = RRall_omis_m1;
        comaUO.RRall_omis = RRall_omis;
        comaUO.RRall_omis_p1 = RRall_omis_p1;
        comaUO.RRall_omis_p2 = RRall_omis_p2;
        if exclude_outliers==1
            save comaUO_CR_all comaUO
        else
            save comaUO_CR_withoutliers_all comaUO
        end
    end
end