%% Description
%
% The following script performs a linear multivariate decoding analysis on
% sound omission and baseline EEG trials for the cardio audio experiment 
% in comatose patients
%
% %
% Requires modification of data paths according to user
%
% Dependencies
% roc2.mat = function for the computation of the ROC curve
% rocarea2.m = function for the computation of the AUC
%
% %
% Inputs
%
% FOR EACH EXPERIMENTAL CONDITION AND PARTICIPANT
% data variable is a fieldtrip data structure containing the cell array data.trial
% with n trials where each cell is a 2-D matrix as channels x timeframes
% where each trial represents the evoked response between -100 ms to 500 ms
% (sampled at 1200 Hz) locked to the R peak during the omission period in the synch
% or the R peaks during the baseline
% e.g. 65 channels x 720 timeframes
%
%
% %
% Output
% fo_acc & fo_auc or uo_acc & uo_ auc = 3-D matrices for accuracy and AUC of the linear SVM 
% classifier as permutations x subjects x timepoints 
% e.g. 501 permutations x 31 subjects x 720 timeframes
% where the 1st permutation is the true unpermutated case
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
    addpath /mnt/HDD1/CardioAudio_coma/scripts/ScriptsLinux/;
    addpath /mnt/HDD1/CardioAudio_coma/data/;
    
    mpath=['/mnt/HDD1/CardioAudio_coma/data/group_trials/'];
    mpath_analysis=['/mnt/HDD1/CardioAudio_coma/analysis/LC/'];

    load '/mnt/HDD1/CardioAudio_coma/data/subj_list.mat' %variable containing list of valid fo and uo subjects
elseif OSflag==2
    addpath D:\CardioAudio_coma\scripts\Scripts\;
    addpath D:\CardioAudio_coma\data\;
    
    mpath=['D:\CardioAudio_coma\data\group_trials\'];
    mpath_analysis=['D:\CardioAudio_coma\analysis\LC\'];

    load 'D:\CardioAudio_coma\data\data\subj_list.mat' %variable containing list of valid fo and uo subjects
end

%valid subjects
sub_type=1; %1 for FO subjects; 2 for UO subjects
if sub_type==1
    subjs=fo_subjs; %list of FO subjects included in the study
elseif sub_type==2
    subjs=uo_subjs; %list of UO subjects included in the study
end

%set linear classification parameters
epochlength=720; %epoch length
train_epochs=120; %number of training epochs
numperm=500; %number of permutations to run
nchannels=62; %number of EEG channels

for s=1:length(subjs)
    cd(mpath) %data path
    %filenames for each condition comparison
    filename1 = sprintf('synchOHEP-s%d.mat',subjs(s)); %synch condition data structure
    filename2 = sprintf('bslOHEP-s%d.mat',subjs(s)); %baseline condition data structure
    
    %load condition data
    load(filename1,'data');
    data1=data;
    clear data filename1
    load(filename2,'data');
    data2=data;
    clear data filename2
    
    %identify number of test epochs
    min_epochs=min([numel(data1.trial) numel(data2.trial)]);
    test_epochs = min_epochs-train_epochs;
    clear min_epochs
    
    %select training and test data from both datasets
    cnt1=1;
    for i=1:train_epochs
        train1(cnt1,:,:)=data1.trial{1,i}(1:nchannels,:);  %training dataset for the first condition
        train2(cnt1,:,:)=data2.trial{1,i}(1:nchannels,:);  %training dataset for the second condition
        cnt1=cnt1+1;
    end
    cnt2=1;
    for j=1+train_epochs:train_epochs+test_epochs
        test1(cnt2,:,:)=data1.trial{1,j}(1:nchannels,:);  %test dataset for the first condition
        test2(cnt2,:,:)=data2.trial{1,j}(1:nchannels,:);  %test dataset for the second condition
        cnt2=cnt2+1;
    end
    clear data1 data2 cnt1 cnt2
    
    %initiate counter for permutations
    cnt_perm=1;

    %loop over timepoints
    for tp=1:epochlength
        %permutation testing
        for perm=1:(numperm+1) %loop over permutations + the unpermuted case (1)
            %train data
            xx1train=squeeze(train1(:,:,tp));
            xx2train=squeeze(train2(:,:,tp));
            
            %test data
            xx1test=squeeze(test1(:,:,tp));
            xx2test=squeeze(test2(:,:,tp));
            
            %prepare final training data
            xxtrain=[xx1train; xx2train];
            if perm==1 %generate training labels for unpermutated case
                y=[zeros(1,train_epochs) zeros(1,train_epochs)+1];
            else %generate permuted training labels
                rng(cnt_perm,'twister')
                y=zeros(1,size(xxtrain,1));
                ix=randperm(numel(y));
                ix=ix(1:train_epochs);
                y(ix)=1;
                clear ix
                cnt_perm=cnt_perm+1;
            end
            
            %Generate linear svm model based on training data and labels
            Rsvm=fitclinear(xxtrain,y);
            %Predict labels for test data
            [Label,score] = predict(Rsvm,[xx1test',xx2test']');
            
            %Compute accuracy and auc based on test data prediction
            if sub_type==1
                fo_acc(perm,s,tp)=(length(find(Label(1:size(xx1test,1)) == 0))+length(find(Label(size(xx1test,1)+1:2*size(xx1test,1)) == 1)))/(2*size(xx1test,1));
            elseif sub_type==2
                uo_acc(perm,s,tp)=(length(find(Label(1:size(xx1test,1)) == 0))+length(find(Label(size(xx1test,1)+1:2*size(xx1test,1)) == 1)))/(2*size(xx1test,1));
            end
            dfcn=score(:,2)-score(:,1);
            [n1,n2,FP,FN] = roc2([dfcn(1:size(xx1test,1)); dfcn(size(xx1test,1)+1:2*size(xx1test,1))]',[zeros(1,size(xx1test,1))+2 zeros(1,size(xx1test,1))+1]);      
            TP = 1 - FN;
            if sub_type==1
                fo_auc(perm,s,tp)=rocarea2(FP, TP);
            elseif sub_type==2
                uo_auc(perm,s,tp)=rocarea2(FP, TP);
            end
            clear xx1train xx2train xx1test xx2test xxtrain y Rsvm Label score dfcn n1 n2 FP FN TP
        end %perm
	end %tp
    clear train1 train2 test1 test2 test_epochs
end %s

cd(mpath_analysis)
if sub_type==1
    save fo_synchVSbsl_linearSVM fo_acc fo_auc
elseif sub_type==2
    save uo_synchVSbsl_linearSVM uo_acc uo_auc
end