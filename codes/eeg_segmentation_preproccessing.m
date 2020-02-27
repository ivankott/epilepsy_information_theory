function [struct_w_seizures_for_cspmvarica]=eeg_segmentation_preproccessing(dir_in, time_before_baseline, time_after_baseline, time_before_preictal, time_after_postictal, do_windows)
    % dir_in = '/home/administrator/Downloads/MATLAB_R2019a/scripts/focal_seizures/';
    % folder with all "d" structures stored in .mat
    
    % time_before_baseline = 60; %seconds
    % time_after_baseline = 60; %seconds
    % time_before_preictal = 30; %seconds
    % time_after_postictal = 30; %seconds
   
    % do_windows = True
    
    if nargin < 2 || isempty(time_before_baseline) || isempty(time_after_baseline) || isempty(time_before_preictal) || isempty(time_after_postictal) 
        time_before_baseline = 60; %seconds
        time_after_baseline = 60; %seconds
        time_before_preictal = 30; %seconds
        time_after_postictal = 30; %seconds
    end
    if nargin < 7 || isempty(do_windows)
        do_windows = true;
    end
    
    len_classes = 30; % all classes are 30s length
    win_len = 5; %seconds
    %samples_win_len = win_len*Fs - 1;
    
    struct_w_seizures_for_cspmvarica = struct;

    fileStruct = dir(strcat(dir_in,'*.mat'));
    
    for i=1:length(fileStruct)% all mat files
        fileStruct(i).name
        load([dir_in,fileStruct(i).name]);

        if ~isfield(d, 'seizureStart')% no seizures in the signal
            disp('No seizures in the signal.')
            continue;
        elseif isempty(d.seizureStart)
            continue;
        elseif ~isfield(d, 'seizureEnd')% no seizures in the signal
            continue;
        elseif length(d.seizureStart)~=length(d.seizureEnd)
            error('Number of Starts and Ends of seizure are mismatched!!!');
            continue;
        else   
            S=length(d.seizureStart);% number of seizures
            
            sigs = struct;

            wheredot=strfind(ans, '.');
            a_name=ans(1:wheredot-1);
            
            %% Let's remove some worthless electrodes
            [n_eeg, ~] = regexp(d.labels, 'EEG', 'match', 'split');
            n_eeg_ind = find(~cellfun('isempty', n_eeg));

            % let's delete them
            [to_delete, ~] = regexp(d.labels, 'A1-N', 'match', 'split'); 
            ind_to_delete = find(~cellfun('isempty', to_delete));
            if ~isempty(ind_to_delete)
                n_eeg_ind(n_eeg_ind == ind_to_delete) = [];
            end

            [to_delete, ~] = regexp(d.labels, 'A1-A2', 'match', 'split');  
            ind_to_delete = find(~cellfun('isempty', to_delete));
            if ~isempty(ind_to_delete)
                n_eeg_ind(n_eeg_ind == ind_to_delete) = [];
            end

            [to_delete, ~] = regexp(d.labels, 'A2-A1', 'match', 'split'); 
            ind_to_delete = find(~cellfun('isempty', to_delete));
            if ~isempty(ind_to_delete)
                n_eeg_ind(n_eeg_ind == ind_to_delete) = [];
            end
            clear to_delete ind_to_delete
            
            data_eeg = [];
            new_labels = [];
            % and save to new data_filtfilted
            for n_i = 1:length(n_eeg_ind)
                data_eeg(:,n_i) = d.data(:, n_eeg_ind(n_i));
                new_labels{n_i} = d.labels{n_i};
            end
            
            Fs = d.Fs;
            
            samples_before_seizure_60_30 = time_before_baseline * Fs - 1;
            samples_after_seizure_30_60 = time_after_baseline * Fs - 1;

            samples_before_seizure_30_0 = time_before_preictal * Fs - 1;
            samples_after_seizure_0_30 = time_after_postictal * Fs - 1;

            for s=1:S % for every seizure in the signal
                disp(['Processing seizure ',num2str(s), ' of ',num2str(S)])
                % calculation of indicies for seizures
                %%
                %if contains(a_name, 'boyko_18') & 
                %%
                if S==1% single seizure in the record
                    s_tmp=find(d.time>=d.seizureStart(s));
                    s_st=s_tmp(1);% start of current seizure (index)
                    s_prev=1;% end of previous seizure (index) = start of record in this case
                    s_next=length(d.time);% start of next seizure (index) = end of record in this case
                    s_tmp=find(d.time>=d.seizureEnd(s));
                    s_end=s_tmp(1)-1;% end of the current seizure (index)

                    if samples_before_seizure_60_30<length(s_prev:s_st)
                        ind_preictal_60_30 = s_st - samples_before_seizure_60_30:s_st - samples_before_seizure_30_0 - 1;
                    else
                        disp('seizure oneset is in the first 1 min of recording or time between seizerus is less then 1 minutes')
                        ind_preictal_60_30 = [];
                    end

                    if samples_before_seizure_30_0<length(s_prev:s_st)
                        ind_preictal_30_0 = s_st - samples_before_seizure_30_0:s_st;
                    else
                        disp('seizure oneset is in the first 30 sec of recording or time between seizerus is less then 30 sec')
                        ind_preictal_30_0 = [];
                    end

                    %ind_during = s_st:s_end;% during seizure

                    if samples_after_seizure_0_30<length(s_end:s_next)
                        ind_postictal_0_30 = s_end:s_end + samples_after_seizure_0_30;
                    else
                        disp('seizure oneset is in 1 min after previous or in the end of recording')
                        ind_postictal_0_30 = [];
                    end

                    if samples_after_seizure_30_60<length(s_end:s_next)
                        ind_postictal_30_60 = s_end + samples_after_seizure_0_30 + 1:s_end + samples_after_seizure_30_60;
                    else
                        disp('seizure oneset is in 1 min after previous or in the end of recording')
                        ind_postictal_30_60 = [];
                    end

                    sigs(s).sigBefore_60_30 = data_eeg(ind_preictal_60_30,:);
                    sigs(s).sigBefore_30_0 = data_eeg(ind_preictal_30_0,:);
                    %sigs(s).sigDuring = data_eeg(ind_during,:);
                    sigs(s).sigAfter_0_30 = data_eeg(ind_postictal_0_30,:); 
                    sigs(s).sigAfter_30_60 = data_eeg(ind_postictal_30_60,:);
                    
    %%
                else % more than one seizure,S>1
                    if s==1% first seizure
                        s_tmp=find(d.time>=d.seizureStart(s));
                        s_st=s_tmp(1);% start of current seizure (index)
                        s_prev=1;% end of previous seizure (index) = start of recording
                        s_tmp=find(d.time>=d.seizureEnd(s));
                        s_end=s_tmp(1)-1;%  end of the current seizure
                        s_tmp=find(d.time>=d.seizureStart(s+1));
                        s_next=s_tmp(1);% start of next seizure (index)
                    elseif s==S% last seizure
                          s_tmp=find(d.time>=d.seizureStart(s));
                          s_st=s_tmp(1);% start of current seizure (index)
                          s_tmp=find(d.time>=d.seizureEnd(s-1));
                          s_prev=s_tmp(1);% end of previous seizure (index)
                          s_tmp=find(d.time>=d.seizureEnd(s));
                          s_end=s_tmp(1)-1;% end of current seizure
                          s_next=length(d.time);% start of next seizure (index) = end of signal  
                    else% middle seizure
                        s_tmp=find(d.time>=d.seizureStart(s));
                        s_st=s_tmp(1);% start of current seizure (index)
                        s_tmp=find(d.time>=d.seizureEnd(s-1));
                        s_prev=s_tmp(1);% end of previous seizure (index)
                        s_tmp=find(d.time>=d.seizureEnd(s));
                        s_end=s_tmp(1)-1;% end of current seizure
                        s_tmp=find(d.time>=d.seizureStart(s+1));
                        s_next=s_tmp(1);% start of next seizure (index)
                    end

                    if samples_before_seizure_60_30<length(s_prev:s_st)
                        ind_preictal_60_30 = s_st - samples_before_seizure_60_30:s_st - samples_before_seizure_30_0 - 1;
                    else
                        disp('seizure oneset is in the first 1 min of recording or time between seizerus is less then 1 minutes')
                        ind_preictal_60_30 = [];
                    end

                    if samples_before_seizure_30_0<length(s_prev:s_st)
                        ind_preictal_30_0 = s_st - samples_before_seizure_30_0:s_st;
                    else
                        disp('seizure oneset is in the first 30 sec of recording or time between seizerus is less then 30 sec')
                        ind_preictal_60_30 = [];
                    end

                    %ind_during = s_st:s_end;% during seizure

                    if samples_after_seizure_0_30<length(s_end:s_next)
                        ind_postictal_0_30 = s_end:s_end + samples_after_seizure_0_30;
                    else
                        disp('seizure oneset is in 1 min after previous or in the end of recording')
                        ind_postictal_0_30 = [];
                    end

                    if samples_after_seizure_30_60<length(s_end:s_next)
                        ind_postictal_30_60 = s_end + samples_after_seizure_0_30 + 1:s_end + samples_after_seizure_30_60;
                    else
                        disp('seizure oneset is in 1 min after previous or in the end of recording')
                        ind_postictal_30_60 = [];
                    end

                    %sigs(s).sigBaseline=downsample(data_baseline_filtfilted,2);
                    sigs(s).sigBefore_60_30 = data_eeg(ind_preictal_60_30,:);
                    sigs(s).sigBefore_30_0 = data_eeg(ind_preictal_30_0,:);
                    %sigs(s).sigDuring = data_eeg(ind_during,:);
                    sigs(s).sigAfter_0_30 = data_eeg(ind_postictal_0_30,:); 
                    sigs(s).sigAfter_30_60 = data_eeg(ind_postictal_30_60,:);
                end
            end
        end

        sigs_windows = struct;
        sigs_for_unmixing_matrices = struct;
        
        Fn = d.Fs/2;
        
        clear d 
        
        fields = fieldnames(sigs);
        %fields(strncmpi(fields,'sigDuring',9)) = [];
        for dilyanka = 1:numel(fields)
            %sig_bez_seiz_st = 1; sig_bez_seiz_end = 3750;

            for n_trial=1:size(sigs, 2)
                figure('units','normalized','outerposition',[0 0 1 1]);

                suptitle_name = [fields{dilyanka}, ' seizure № ', int2str(n_trial), ' of patient ', a_name];
                suptitle(suptitle_name)

                subplot(4,2,1);plot(sigs(n_trial).(fields{dilyanka})(:,4)); title(new_labels{4}); %xlim([sig_bez_seiz_st sig_bez_seiz_end])
                subplot(4,2,3);plot(sigs(n_trial).(fields{dilyanka})(:,6)); title(new_labels{6}); %xlim([sig_bez_seiz_st sig_bez_seiz_end])
                subplot(4,2,5);plot(sigs(n_trial).(fields{dilyanka})(:,14)); title(new_labels{14}); %xlim([sig_bez_seiz_st sig_bez_seiz_end])
                subplot(4,2,7);plot(sigs(n_trial).(fields{dilyanka})(:,16)); title(new_labels{16}); %xlim([sig_bez_seiz_st sig_bez_seiz_end])
                
                subplot(4,2,2);plot(sigs(n_trial).(fields{dilyanka})(:,8)); title(new_labels{8}); 
                subplot(4,2,4);plot(sigs(n_trial).(fields{dilyanka})(:,12)); title(new_labels{12}); 
                subplot(4,2,6);plot(sigs(n_trial).(fields{dilyanka})(:,13)); title(new_labels{13}); 
                subplot(4,2,8);plot(sigs(n_trial).(fields{dilyanka})(:,17)); title(new_labels{17}); 
                
                %{
                subplot(4,2,2);
                n = length(sigs(n_trial).(fields{dilyanka})(:,4));
                nfft = 512; 
                nfft2 = nfft*2;
                px=fft(sigs(n_trial).(fields{dilyanka})(:,4), nfft2);
                px=2*abs(px(1:nfft))/n; 
                f=(0:Fs/(2*(nfft-1)):Fs/2)'; 
                plot(f, px); 
                title('<-- Power spectrum'); %xlim([0 125])

                subplot(4,2,4);
                n = length(sigs(n_trial).(fields{dilyanka})(:,6));
                nfft = 512; 
                nfft2 = nfft*2;
                px=fft(sigs(n_trial).(fields{dilyanka})(:,6), nfft2);
                px=2*abs(px(1:nfft))/n; 
                f=(0:Fs/(2*(nfft-1)):Fs/2)'; 
                plot(f, px);
                title('<-- Power spectrum'); %xlim([0 125])

                subplot(4,2,6);
                n = length(sigs(n_trial).(fields{dilyanka})(:,14));
                nfft = 512; 
                nfft2 = nfft*2;
                px=fft(sigs(n_trial).(fields{dilyanka})(:,14), nfft2);
                px=2*abs(px(1:nfft))/n; 
                f=(0:Fs/(2*(nfft-1)):Fs/2)'; 
                plot(f, px);
                title('<-- Power spectrum'); %xlim([0 125])

                subplot(4,2,8);
                n = length(sigs(n_trial).(fields{dilyanka})(:,16));
                nfft = 512; 
                nfft2 = nfft*2;
                px=fft(sigs(n_trial).(fields{dilyanka})(:,16), nfft2);
                px=2*abs(px(1:nfft))/n; 
                f=(0:Fs/(2*(nfft-1)):Fs/2)'; 
                plot(f, px);
                title('<-- Power spectrum'); %xlim([0 125])
                %}
                
                name_trial = [a_name, '_seizure_', int2str(n_trial), '_', fields{dilyanka}, '.jpg'];

                prompt_name = strcat('Do you want to keep this trial? YES -- 1 / NO -- 0:\n', name_trial);
                prompt = {prompt_name};
                dlgtitle = 'Do you want to keep this trial?';
                dims = [1 40];
                definput = {'1'};
                keep = inputdlg(prompt,dlgtitle,dims,definput);
                keep = str2double(keep);

                if keep == 1
                    saveas(gcf, name_trial)
                    for n_wind=0:len_classes - win_len                                     
                        %% low pass 
                        Wp = 30/Fn;
                        Ws = 45/Fn;
                        Rp = 3;                                                        % Passband Ripple (dB)
                        Rs = 40;
                        [n,wn] = buttord(Wp, Ws, Rp, Rs);
                        [B,A] = butter(n, wn);
                        %data_baseline = filtfilt(B, A, d.data(ind_baseline,:));
                        data_filfit_low = filtfilt(B, A, sigs(n_trial).(fields{dilyanka})(n_wind*Fs + 1: (n_wind+5)*Fs, :));

                        %% high pass
                        Wp = 1/Fn;
                        Ws = 0.15/Fn;
                        Rp = 3;                                                        % Passband Ripple (dB)
                        Rs = 70;
                        [n,wn] = buttord(Wp, Ws, Rp, Rs);
                        [B,A] = butter(n, wn, 'high');
                        %data_baseline_filtfilted = filtfilt(B, A, data_baseline);
                        data_filfit_low_high = filtfilt(B, A, data_filfit_low);

                        sigs_windows(n_trial).(fields{dilyanka}){n_wind+1} = downsample(data_filfit_low_high, 2);
                        
                        if contains(fields{dilyanka}, 'Before') && n_wind == len_classes - win_len
                            sigs_for_unmixing_matrices(n_trial).(fields{dilyanka}) = sigs_windows(n_trial).(fields{dilyanka}){n_wind+1};
                        elseif  contains(fields{dilyanka}, 'After') && n_wind == 0
                            sigs_for_unmixing_matrices(n_trial).(fields{dilyanka}) = sigs_windows(n_trial).(fields{dilyanka}){n_wind+1};
                        end 
                    end
                elseif keep == 0
                    disp('Ok, deleted')
                    %sigs(n_trial).(fields{dilyanka}) = [];
                else
                    disp('Guess, it"s no. Trial was deleted')
                    %sigs(n_trial).(fields{dilyanka}) = [];
                end
                close all
            end
        end

        if length(fieldnames(sigs_windows)) ~= 0
            if do_windows == true
                struct_w_seizures_for_cspmvarica.(a_name).sigs_windows = sigs_windows;
            end
            struct_w_seizures_for_cspmvarica.(a_name).labels = new_labels;
            struct_w_seizures_for_cspmvarica.(a_name).sigs_unmixing = sigs_for_unmixing_matrices;
            struct_w_seizures_for_cspmvarica.(a_name).Fs = Fs;
        end
        clear sigs_windows sigs_for_unmixing_matrices new_labels
    end
    save('struct_w_seizures_for_cspmvarica.mat','struct_w_seizures_for_cspmvarica', '-v7.3')
end

