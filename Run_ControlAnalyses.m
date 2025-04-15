%% Description
%
% The following script performs the control analyses for the cardio-audio
% experiment in comatose patients
% SR, RS, SS, RR intervals and variabilities are calculated to assess
% whether the experimental design was successfuly administered
%
% %
% Requires modification of data paths according to user
%
% %
% Inputs
%
% FOR EACH AUDITORY CONDITION AND PARTICIPANT
% T variable is a table array of dimension Yx4, where Y is the number of
% trials per experimental block and 4 columns as follows:
% T.trial = trial number per block
% T.cond = binary variable (1 for sound trial; 2 for omission trial)
% T.lsound = actual sound trigger onset in samples
% T.lr = closest r peak onset in samples (to stimulus onset)
% loaded in this script as e.g. '../trg_S1.mat'
% example
% trial     cond    lsound      lr
% 1         1       9173        9105
% 2         1       10064       9999
% 3         1       10984       10919
% 4         1       11945       11875
% 5         2       NaN         13895
%
% FOR BASELINE CONDITION AND PARTICIPANT
% rr variable (Nx1 vector) containing all RR intervals
% calculated as the difference between all offline detected time stamps of
% R peaks (from raw ECG data timeseries)
% loaded in this script as e.g. ..'/ecg_B1.mat'
%
% %
% Output
%
% comaFO.SR = mean SR interval
% comaFO.RS = mean RS interval
% comaFO.SS = mean SS interval
% comaFO.RR = mean RR interval
% comaFO.SR, comaFO.RS, comaFO.SS variable (3xN matrix) with rows for each
% condition (synch, asynch, isoch) and columns for N sujects
% comaFO.RR variable (4xN matrix) with rows for each
% condition (synch, asynch, isoch, baseline) and columns for N sujects
%
% comaFO.varSR = mean SR variability
% comaFO.varRS = mean RS variability
% comaFO.varSS = mean SS variability
% comaFO.varRR = mean RR variability
% variability refers to SEM of corresponding intervals
% comaFO.varSR, comaFO.varRS, comaFO.varSS variable (3xN matrix) with rows for each
% condition (synch, asynch, isoch) and columns for N sujects
% comaFO.varRR variable (4xN matrix) with rows for each
% condition (synch, asynch, isoch, baseline) and columns for N sujects
%
% comaFO.SRall = SRall; %all SR interval values
% comaFO.RSall = RSall; %all RS interval values
% comaFO.SSall = SSall; %all SS interval values
% comaFO.RRall = final_RRall; %all RR interval values
% comaFO.SRall, comaFO.RSall, comaFO.SSall cell structure (3xN) with rows for each
% condition (synch, asynch, isoch) and columns for N sujects
% comaFO.RRall cell structure (4xN) with rows for each
% condition (synch, asynch, isoch, baseline) and columns for N sujects
% each cell is a Yx1 vector containing all available interval values
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
%% Code
clear all;

OSflag=1; %1 = Linux; 2 = Windows
% set path
if OSflag==1
    mpath=['/mnt/HDD1/CardioAudio_coma/data/'];
    mpath_analysis=['/mnt/HDD1/CardioAudio_coma/analysis/CA/'];
    load '/mnt/HDD1/CardioAudio_coma/data/subj_list.mat' %variable containing list of valid fo and uo subjects
elseif OSflag==2
    mpath=['D:\CardioAudio_coma\data\'];
    mpath_analysis=['D:\CardioAudio_coma\analysis\CA\'];
    load 'D:\CardioAudio_coma\data\data\subj_list.mat' %variable containing list of valid fo and uo subjects
end

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

            Ts=[]; %sound onset
            SoundR=[]; %sound to R interval
            dTs=[]; % sound to sound interval
            RSound=[]; %R to sound interval
            RRint=[]; %%R to R interval
            condstage=[]; %block number

            for block=block_num

                if OSflag==1
                	load([mpath,'s',num2str(subjs(j)),'/process/trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
                elseif OSflag==2
                    load([mpath,'s',num2str(subjs(j)),'\process\trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
                end

                %prepare T variable
                s = T.lsound(T.cond==1);
                ss = diff(T.lsound);
                ss_avg = nanmean(ss);

                T.lomiss = nan(size(T.lsound));
                ind = find(T.cond==2);
                ind=ind(ind>1);
                T.lomiss(ind) = T.lsound(ind-1)+ss_avg;
                clear ind

                %compute different control analyses variables
                Ts=[Ts; T.lsound]; %S latencies

                io=find(T.cond == 2);
                io(1)=[];

                SoundR=[SoundR; T.lr(io)-T.lsound(io-1)]; %SR

                ir=find(T.cond == 1);
                ir(1)=[];

                RSound=[RSound; T.lsound(ir)-T.lr(ir)]; %RS

                RRint=[RRint; diff(T.lr)]; %RR

                dTs=[dTs; diff(T.lsound)]; %SS

                clear io ir

            end

            if exclude_outliers==1
                %identify and remove outliers
                dTsplot=dTs; %SS outliers - currently not used further
                dTsplot(find(isnan(dTs) == 1))=[];
                [dTsClean,II1]=rmoutliers(dTsplot);
                SoundRplot=SoundR; %SR outliers
                SoundRplot(find(isnan(SoundR) == 1))=[];
                [SoundRClean,II2]=rmoutliers(SoundRplot);
                RSoundplot=RSound; %RS outliers
                RSoundplot(find(isnan(RSound) == 1))=[];
                [RSoundClean,II3]=rmoutliers(RSoundplot);
                RRintplot=RRint; %RR outliers
                RRintplot(find(isnan(RRint) == 1))=[];
                [RRintClean,II4]=rmoutliers(RRintplot);

            end


            if exclude_outliers==1
                %keep all values
                SRall{k,j}=SoundRClean;
                RSall{k,j}=RSoundClean;
                SSall{k,j}=dTsClean;
                RRall{k,j}=RRintClean;

                SR(k,j)=nanmean(SoundRClean); %SR interval
                varSR(k,j)=nanstd(SoundRClean); %SR variability

                SS(k,j)=nanmean(dTsClean); %SS interval
                varSS(k,j)=nanstd(dTsClean); %% SS variability

                RS(k,j)=nanmean(RSoundClean); %RS interval
                varRS(k,j)=nanstd(RSoundClean); %RS variability

                RR(k,j)=nanmean(RRintClean); %RR interval
                varRR(k,j)=nanstd(RRintClean); %RR variability


            else %if exclude_outliers==0
                SRall{k,j}=SoundR;
                SSall{k,j}=dTs;
                RSall{k,j}=RSound;
                RRall{k,j}=RRint;

                SR(k,j)=nanmean(SoundR); %SR interval
                varSR(k,j)=nanstd(SoundR); %SR variability

                SS(k,j)=nanmean(dTs); %SS interval
                varSS(k,j)=nanstd(dTs); %SS variability

                RS(k,j)=nanmean(RSound); %RS interval
                varRS(k,j)=nanstd(RSound); %RS variability

                RR(k,j)=nanmean(RRint); %RR interval
                varRR(k,j)=nanstd(RRint); %RR variability

            end

            clear Ts Tr dTs SoundR RSound RRint dTsClean SoundRClean RSoundClean RRintClean

        end
    end

    %include baseline data for RR info
    for j=1:length(subjs)

        clear irr
        if OSflag==1
            load([mpath,'s',num2str(subjs(j)),'/process/ecg_B.mat']) % ECG info before artefact rejection
        elseif OSflag==2
            load([mpath,'s',num2str(subjs(j)),'\process\ecg_B.mat']) % ECG info before artefact rejection
        end

        if exclude_outliers==1
            %identify and remove outliers
         	rrbplot=rr;
         	rrbplot(find(isnan(rr) == 1))=[];
         	[rrbClean,II1]=rmoutliers(rrbplot);

            %store baseline info
            RRall{4,j}=rrbClean;
            RR(4,j)=nanmean(RRbClean); %RR interval
            varRR(4,j)=nanstd(RRbClean); %RR variability

        else %if exclude outliers==0
            %store baseline info
            RRall{4,j}=rr;
            RR(4,j)=nanmean(rr); %RR interval
            varRR(4,j)=nanstd(rr); %RR variability
        end

        clear rr
    end
    
    cd(mpath_analysis) %output path for saving
    if sub_type==1 %FO subjects
        %store relevant data for statistical analysis
        comaFO.SR = SR; %mean SR interval
        comaFO.RS = RS; %mean RS interval
        comaFO.SS = SS; %mean SS interval
        comaFO.RR = RR; %mean RR interval
        comaFO.SRall = SRall; %all SR interval values
        comaFO.RSall = RSall; %all RS interval values
        comaFO.SSall = SSall; %all SS interval values
        comaFO.RRall = RRall; %all RR interval values
        comaFO.varSR = varSR; %mean SR variability
        comaFO.varRS = varRS; %mean RS variability
        comaFO.varSS = varSS; %mean SS variability
        comaFO.varRR = varRR; %mean RR variability
        if exclude_outliers==1
            save comaFO_CA_all comaFO
        else
            save comaFO_CA_withoutliers_all comaFO
        end
    elseif sub_type==2 %UO subjects
        comaUO.SR = SR; %mean SR interval
        comaUO.RS = RS; %mean RS interval
        comaUO.SS = SS; %mean SS interval
        comaUO.RR = final_RR; %mean RR interval
        comaUO.SRall = SRall; %all SR interval values
        comaUO.RSall = RSall; %all RS interval values
        comaUO.SSall = SSall; %all SS interval values
        comaUO.RRall = final_RRall; %all RR interval values
        comaUO.varSR = varSR; %mean SR variability
        comaUO.varRS = varRS; %mean RS variability
        comaUO.varSS = varSS; %mean SS variability
        comaUO.varRR = final_varRR; %mean RR variability
        if exclude_outliers==1
            save comaUO_CA_all comaUO
        else
            save comaUO_CA_withoutliers_all comaUO
        end
    end %sub_type
end %compute_means