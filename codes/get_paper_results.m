% SCOT FUNCTIONALITY 
clc
close all 
clear variables
addpath('/Users/ivan_kotiuchyi/Documents/PyCharmProjects/epilepsy_information_theory/Functions_MVAR');
addpath('/Users/ivan_kotiuchyi/Documents/PyCharmProjects/epilepsy_information_theory/functions_');
addpath(genpath('/Users/ivan_kotiuchyi/Documents/PyCharmProjects/epilepsy_information_theory/eeglab'));

% SET UP DATA
load('/Users/ivan_kotiuchyi/Documents/PyCharmProjects/epilepsy_information_theory/datafocal_data_filtered_for_SPS.mat')

names = fieldnames(focal_data_filtered_for_SPS);
t_trial_basepre = 0; t_trial_pre = 0; t_trial_post = 0; t_trial_basepost = 0;

pmax=40; % max for model order selection
pcrit='b'; % 'a'-AIC, 'b'-BIC, 'c'-impoнаsed
pimp=10;
valore_p_value_crit = 0.05;

det_norm_eeg='y'; %detrend and normalization
det_norm_sources='n'; %detrend and normalization
pfilter=0.92; %detrend highpass filter

for subj_name_ind = 1:numel(names) %9
    disp(names{subj_name_ind})
    for seiz_ind = 1:length(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs) %N_seizures
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigBaseline') && ...
                ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBaseline)
            t_trial_basepre = t_trial_basepre + 1;
            X_sigBaseline(t_trial_basepre,:,:) = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBaseline(:,1:19)';   
        end
        
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigBefore') && ...
                ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBefore)
            t_trial_pre = t_trial_pre + 1;
            X_sigBefore(t_trial_pre,:,:) = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBefore(:,1:19)';   
        end
        
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigAfter') && ...
                ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigAfter)
            t_trial_post = t_trial_post + 1;
            X_sigAfter(t_trial_post,:,:) = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigAfter(:,1:19)';   
        end
        
%         if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'EEG_10_min_after') && ...
%                 ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).EEG_10_min_after)
%             t_trial_basepost = t_trial_basepost + 1;
%             X_eeg_10_min_after(t_trial_basepost,:,:) = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).EEG_10_min_after(:,1:19)';   
%         end
    end
end

disp(["There are" t_trial_basepre "trials in EEG 10_min_BEFORE"]) 
disp(["There are" t_trial_pre "trials in EEG 10_sec_BEFORE"])
disp(["There are" t_trial_post "trials in EEG 10_sec_AFTER"]) 
%disp(["There are" t_trial_basepost "trials in EEG 10_min_AFTER"])

X_sigBaseline = permute(X_sigBaseline, [2 3 1]);
X_sigBefore = permute(X_sigBefore, [2 3 1]);
X_sigAfter = permute(X_sigAfter, [2 3 1]);
%X_eeg_10_min_after = permute(X_eeg_10_min_after, [2 3 1]);

% base pre
[NumComponents_base_pre]=csp_numcomponents_rim_onlynum(X_sigBaseline, X_sigBefore, 0.95, 0);
[NumComponents_base_post]=csp_numcomponents_rim_onlynum(X_sigBaseline, X_sigAfter, 0.95, 0);
[NumComponents_pre_post]=csp_numcomponents_rim_onlynum(X_sigBefore, X_sigAfter, 0.95, 0);

NumComponents_All=[NumComponents_base_pre, NumComponents_base_post, NumComponents_pre_post];
NumComponents=max(NumComponents_All);

disp(['Number of components chosen: ', num2str(NumComponents)]);

UX_base_pre = funct_CSPVARICA_InfomaxICA_numComp(X_sigBaseline, X_sigBefore, NumComponents, 27, 1);
% base post
UX_base_post = funct_CSPVARICA_InfomaxICA_numComp(X_sigBaseline, X_sigAfter, NumComponents, 27, 0);
% pre post
UX_pre_post = funct_CSPVARICA_InfomaxICA_numComp(X_sigBefore, X_sigAfter, NumComponents, 27, 0);

% UX_10m_pre_10m_aft = funct_CSPVARICA_InfomaxICA(X_sigBaseline, X_eeg_10_min_after, 0.95, 27, 0);
% UX_10s_pre_10m_aft = funct_CSPVARICA_InfomaxICA(X_sigBefore, X_eeg_10_min_after, 0.95, 27, 0);
% UX_10s_aft_10m_aft = funct_CSPVARICA_InfomaxICA(X_sigAfter, X_eeg_10_min_after, 0.95, 27, 0);

%% BASE PRE INFORMATION THEORY
for subj_name_ind = 1:numel(names) %9
    results_base = struct; %base
    results_pre = struct; %pre
    
    for seiz_ind = 1:length(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs) %N_seizures
        %if contains((names{subj_name_ind}),'')
        %end
        %% sigBaseline
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigBaseline')
            if ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBaseline)
                X = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBaseline(:,1:19);
                [M,N]=size(X);

                X_n = [];
                if det_norm_eeg=='y'
                    for i=1:N
                        Xtmp=AR_filter(X,i,pfilter);
                        X_n(:,i) = (Xtmp - mean(Xtmp))/std(Xtmp);
                    end
                    X = X_n';
                else
                    X = X';
                end
                S = UX_base_pre * X;
                Q=size(S,1);
                [pAIC,pBIC,aic,mdl] = mos_idMVAR(S,pmax,0);
               %{
               % plot
                figure(1); %first nS EEGs
                for is=1:Q
                    subplot(Q,1,is); plot(X(is,:));
                    xlim([1 N]);ylabel(['X' int2str(is)]);
                end

                figure(2); % all nS sources
                for is=1:Q
                    subplot(Q,2,2*(is-1)+1); plot(S(is,:));
                    xlim([1 N]);ylabel(['S' int2str(is)])
                end
                subplot(Q,2,[2 4 6 8 10 12 14 16 18])
                plot(aic,'b.-'); hold on; plot(mdl,'r.-');
                plot(pAIC,aic(pAIC),'bo'); plot(pBIC,mdl(pBIC),'ro');
                legend(['p_{AIC}=' int2str(pAIC)],['p_{BIC}=' int2str(pBIC)]);
                title(['window N ' int2str(wind_ind)])
                %}
                % connectivity based on state space models
                % model order selection
                switch pcrit
                    case 'a', p=pAIC;
                    case 'b', p=pBIC;
                    case 'c', p=pimp;
                end
                %%% detrend and normalization
                if det_norm_sources=='y'
                    for q=1:Q
                        Stmp=AR_filter(S',q,pfilter);
                %         Sn(:,q)=Stmp;
                        Sn(q,:)=(Stmp-mean(Stmp))/std(Stmp);
                    end
                    %Sn=Sn';
                else
                    Sn=S;
                end

                %% partial variances (correlations)
                [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                qlags=10;
                Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                for jj=1:Q
                    for ii=1:Q
                        if ii~=jj
                            ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                            Hj(jj)=ret.Hy;
                            Hj_j(jj)=ret.Hy_y;
                            Hj_all(jj)=ret.Hy_yzx;
                            Hj_jk(jj,ii)=ret.Hy_yz;
                            % t-test
                            Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                            p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                        end
                    end
                end
                Sj=Hj-Hj_j;
                Tj=Hj_j-Hj_all;
                Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                psign=p_val<valore_p_value_crit;

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p = 15;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=10;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind),' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=5;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    disp('ALLO, some shit is still complex even after decreasing MVAR order to 5.')
                end

                results_base(seiz_ind).Sj = Sj;
                results_base(seiz_ind).Tj = Tj;
                results_base(seiz_ind).Tij_k = Tij_k;
                results_base(seiz_ind).p_val = p_val;
                results_base(seiz_ind).psign = psign;  
                results_base(seiz_ind).pBIC = pBIC;   
                results_base(seiz_ind).meanSj = mean(Sj);   
                results_base(seiz_ind).meanTj = mean(Tj);   
                results_base(seiz_ind).meanTij_k = nanmean(Tij_k(:));   
                results_base(seiz_ind).Results_doc = sum(sum(p_val<valore_p_value_crit));
            end
        end 
        
        %% sigBefore
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigBefore')
            if ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBefore)
                X = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBefore(:,1:19);
                [M,N]=size(X);

                X_n = [];
                if det_norm_eeg=='y'
                    for i=1:N
                        Xtmp=AR_filter(X,i,pfilter);
                        X_n(:,i) = (Xtmp - mean(Xtmp))/std(Xtmp);
                    end
                    X = X_n';
                else
                    X = X';
                end
                S = UX_base_pre * X;
                Q=size(S,1);
                [pAIC,pBIC,aic,mdl] = mos_idMVAR(S,pmax,0);
               %{
               % plot
                figure(1); %first nS EEGs
                for is=1:Q
                    subplot(Q,1,is); plot(X(is,:));
                    xlim([1 N]);ylabel(['X' int2str(is)]);
                end

                figure(2); % all nS sources
                for is=1:Q
                    subplot(Q,2,2*(is-1)+1); plot(S(is,:));
                    xlim([1 N]);ylabel(['S' int2str(is)])
                end
                subplot(Q,2,[2 4 6 8 10 12 14 16 18])
                plot(aic,'b.-'); hold on; plot(mdl,'r.-');
                plot(pAIC,aic(pAIC),'bo'); plot(pBIC,mdl(pBIC),'ro');
                legend(['p_{AIC}=' int2str(pAIC)],['p_{BIC}=' int2str(pBIC)]);
                title(['window N ' int2str(wind_ind)])
                %}
                % connectivity based on state space models
                % model order selection
                switch pcrit
                    case 'a', p=pAIC;
                    case 'b', p=pBIC;
                    case 'c', p=pimp;
                end
                %%% detrend and normalization
                if det_norm_sources=='y'
                    for q=1:Q
                        Stmp=AR_filter(S',q,pfilter);
                %         Sn(:,q)=Stmp;
                        Sn(q,:)=(Stmp-mean(Stmp))/std(Stmp);
                    end
                    %Sn=Sn';
                else
                    Sn=S;
                end

                %% partial variances (correlations)
                [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                qlags=10;
                Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                for jj=1:Q
                    for ii=1:Q
                        if ii~=jj
                            ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                            Hj(jj)=ret.Hy;
                            Hj_j(jj)=ret.Hy_y;
                            Hj_all(jj)=ret.Hy_yzx;
                            Hj_jk(jj,ii)=ret.Hy_yz;
                            % t-test
                            Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                            p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                        end
                    end
                end
                Sj=Hj-Hj_j;
                Tj=Hj_j-Hj_all;
                Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                psign=p_val<valore_p_value_crit;

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p = 15;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=10;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind),' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=5;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    disp('ALLO, some shit is still complex even after decreasing MVAR order to 5.')
                end
                
                results_pre(seiz_ind).Sj = Sj;
                results_pre(seiz_ind).Tj = Tj;
                results_pre(seiz_ind).Tij_k = Tij_k;
                results_pre(seiz_ind).p_val = p_val;
                results_pre(seiz_ind).psign = psign;  
                results_pre(seiz_ind).pBIC = pBIC;   
                results_pre(seiz_ind).meanSj = mean(Sj);   
                results_pre(seiz_ind).meanTj = mean(Tj);   
                results_pre(seiz_ind).meanTij_k = nanmean(Tij_k(:));   
                results_pre(seiz_ind).Results_doc = sum(sum(p_val<valore_p_value_crit))   ;          ;
            end
        end
        
    end
    
    results_data.(names{subj_name_ind}).results_base = results_base;
    results_data.(names{subj_name_ind}).results_pre = results_pre;
    disp([names{subj_name_ind} " is done"])
end
save('mvar_source_results_MATLAB_paper_base_pre.mat','results_data', '-v7.3')

%% BASE POST INFORMATION THEORY
for subj_name_ind = 1:numel(names) %9
    results_base = struct; %base
    results_post = struct; %post
    
    for seiz_ind = 1:length(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs) %N_seizures
        %if contains((names{subj_name_ind}),'')
        %end
        %% sigAfter
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigAfter')
            if ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigAfter)
                X = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigAfter(:,1:19);
               
                [M,N]=size(X);

                X_n = [];
                if det_norm_eeg=='y'
                    for i=1:N
                        Xtmp=AR_filter(X,i,pfilter);
                        X_n(:,i) = (Xtmp - mean(Xtmp))/std(Xtmp);
                    end
                    X = X_n';
                else
                    X = X';
                end
                S = UX_base_post * X;
                Q=size(S,1);
                [pAIC,pBIC,aic,mdl] = mos_idMVAR(S,pmax,0);
               
                % connectivity based on state space models
                % model order selection
                switch pcrit
                    case 'a', p=pAIC;
                    case 'b', p=pBIC;
                    case 'c', p=pimp;
                end
                
                %%% detrend and normalization
                if det_norm_sources=='y'
                    for q=1:Q
                        Stmp=AR_filter(S',q,pfilter);
                %         Sn(:,q)=Stmp;
                        Sn(q,:)=(Stmp-mean(Stmp))/std(Stmp);
                    end
                    %Sn=Sn';
                else
                    Sn=S;
                end

                %% partial variances (correlations)
                [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                qlags=10;
                Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                for jj=1:Q
                    for ii=1:Q
                        if ii~=jj
                            ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                            Hj(jj)=ret.Hy;
                            Hj_j(jj)=ret.Hy_y;
                            Hj_all(jj)=ret.Hy_yzx;
                            Hj_jk(jj,ii)=ret.Hy_yz;
                            % t-test
                            Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                            p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                        end
                    end
                end
                Sj=Hj-Hj_j;
                Tj=Hj_j-Hj_all;
                Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                psign=p_val<valore_p_value_crit;
                
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p = 15;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=10;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind),' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=5;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    disp('ALLO, some shit is still complex even after decreasing MVAR order to 5.')
                end

                results_post(seiz_ind).Sj = Sj;
                results_post(seiz_ind).Tj = Tj;
                results_post(seiz_ind).Tij_k = Tij_k;
                results_post(seiz_ind).p_val = p_val;
                results_post(seiz_ind).psign = psign;  
                results_post(seiz_ind).pBIC = pBIC;   
                results_post(seiz_ind).meanSj = mean(Sj);   
                results_post(seiz_ind).meanTj = mean(Tj);   
                results_post(seiz_ind).meanTij_k = nanmean(Tij_k(:));   
                results_post(seiz_ind).Results_doc = sum(sum(p_val<valore_p_value_crit))               ;
            end
        end
        
        %% sigBaseline
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigBaseline')
            if ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBaseline)
                X = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBaseline(:,1:19);
                [M,N]=size(X);

                X_n = [];
                if det_norm_eeg=='y'
                    for i=1:N
                        Xtmp=AR_filter(X,i,pfilter);
                        X_n(:,i) = (Xtmp - mean(Xtmp))/std(Xtmp);
                    end
                    X = X_n';
                else
                    X = X';
                end
                S = UX_base_post * X;
                Q=size(S,1);
                [pAIC,pBIC,aic,mdl] = mos_idMVAR(S,pmax,0);
               %{
               % plot
                figure(1); %first nS EEGs
                for is=1:Q
                    subplot(Q,1,is); plot(X(is,:));
                    xlim([1 N]);ylabel(['X' int2str(is)]);
                end

                figure(2); % all nS sources
                for is=1:Q
                    subplot(Q,2,2*(is-1)+1); plot(S(is,:));
                    xlim([1 N]);ylabel(['S' int2str(is)])
                end
                subplot(Q,2,[2 4 6 8 10 12 14 16 18])
                plot(aic,'b.-'); hold on; plot(mdl,'r.-');
                plot(pAIC,aic(pAIC),'bo'); plot(pBIC,mdl(pBIC),'ro');
                legend(['p_{AIC}=' int2str(pAIC)],['p_{BIC}=' int2str(pBIC)]);
                title(['window N ' int2str(wind_ind)])
                %}
                % connectivity based on state space models
                % model order selection
                switch pcrit
                    case 'a', p=pAIC;
                    case 'b', p=pBIC;
                    case 'c', p=pimp;
                end
                %%% detrend and normalization
                if det_norm_sources=='y'
                    for q=1:Q
                        Stmp=AR_filter(S',q,pfilter);
                %         Sn(:,q)=Stmp;
                        Sn(q,:)=(Stmp-mean(Stmp))/std(Stmp);
                    end
                    %Sn=Sn';
                else
                    Sn=S;
                end

                %% partial variances (correlations)
                [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                qlags=10;
                Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                for jj=1:Q
                    for ii=1:Q
                        if ii~=jj
                            ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                            Hj(jj)=ret.Hy;
                            Hj_j(jj)=ret.Hy_y;
                            Hj_all(jj)=ret.Hy_yzx;
                            Hj_jk(jj,ii)=ret.Hy_yz;
                            % t-test
                            Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                            p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                        end
                    end
                end
                Sj=Hj-Hj_j;
                Tj=Hj_j-Hj_all;
                Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                psign=p_val<valore_p_value_crit;

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p = 15;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=10;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind),' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=5;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    disp('ALLO, some shit is still complex even after decreasing MVAR order to 5.')
                end

                results_base(seiz_ind).Sj = Sj;
                results_base(seiz_ind).Tj = Tj;
                results_base(seiz_ind).Tij_k = Tij_k;
                results_base(seiz_ind).p_val = p_val;
                results_base(seiz_ind).psign = psign;  
                results_base(seiz_ind).pBIC = pBIC;   
                results_base(seiz_ind).meanSj = mean(Sj);   
                results_base(seiz_ind).meanTj = mean(Tj);   
                results_base(seiz_ind).meanTij_k = nanmean(Tij_k(:));   
                results_base(seiz_ind).Results_doc = sum(sum(p_val<valore_p_value_crit))             ;
            end
        end        
    end
    
    results_data.(names{subj_name_ind}).results_base = results_base;
    results_data.(names{subj_name_ind}).results_post = results_post;
    disp([names{subj_name_ind} " is done"])
end
save('mvar_source_results_MATLAB_paper_base_post.mat','results_data', '-v7.3')

%% PRE POST INFORMATION THEORY
for subj_name_ind = 1:numel(names) %9
    results_pre = struct; %pre
    results_post = struct; %post
    
    for seiz_ind = 1:length(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs) %N_seizures
        %if contains((names{subj_name_ind}),'')
        %end
        %% sigAfter
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigAfter')
            if ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigAfter)
                X = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigAfter(:,1:19);
                
                [M,N]=size(X);

                X_n = [];
                if det_norm_eeg=='y'
                    for i=1:N
                        Xtmp=AR_filter(X,i,pfilter);
                        X_n(:,i) = (Xtmp - mean(Xtmp))/std(Xtmp);
                    end
                    X = X_n';
                else
                    X = X';
                end
                S = UX_pre_post * X;
                Q=size(S,1);
                [pAIC,pBIC,aic,mdl] = mos_idMVAR(S,pmax,0);
               
                % connectivity based on state space models
                % model order selection
                switch pcrit
                    case 'a', p=pAIC;
                    case 'b', p=pBIC;
                    case 'c', p=pimp;
                end
                
                %%% detrend and normalization
                if det_norm_sources=='y'
                    for q=1:Q
                        Stmp=AR_filter(S',q,pfilter);
                %         Sn(:,q)=Stmp;
                        Sn(q,:)=(Stmp-mean(Stmp))/std(Stmp);
                    end
                    %Sn=Sn';
                else
                    Sn=S;
                end

                %% partial variances (correlations)
                [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                qlags=10;
                Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                for jj=1:Q
                    for ii=1:Q
                        if ii~=jj
                            ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                            Hj(jj)=ret.Hy;
                            Hj_j(jj)=ret.Hy_y;
                            Hj_all(jj)=ret.Hy_yzx;
                            Hj_jk(jj,ii)=ret.Hy_yz;
                            % t-test
                            Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                            p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                        end
                    end
                end
                Sj=Hj-Hj_j;
                Tj=Hj_j-Hj_all;
                Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                psign=p_val<valore_p_value_crit;
                
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p = 15;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=10;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind),' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=5;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    disp('ALLO, some shit is still complex even after decreasing MVAR order to 5.')
                end

                results_post(seiz_ind).Sj = Sj;
                results_post(seiz_ind).Tj = Tj;
                results_post(seiz_ind).Tij_k = Tij_k;
                results_post(seiz_ind).p_val = p_val;
                results_post(seiz_ind).psign = psign;  
                results_post(seiz_ind).pBIC = pBIC;   
                results_post(seiz_ind).meanSj = mean(Sj);   
                results_post(seiz_ind).meanTj = mean(Tj);   
                results_post(seiz_ind).meanTij_k = nanmean(Tij_k(:));   
                results_post(seiz_ind).Results_doc = sum(sum(p_val<valore_p_value_crit))               ;
            end
        end
        
        %% sigBefore
        if isfield(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind), 'sigBefore')
            if ~isempty(focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBefore)
                X = focal_data_filtered_for_SPS.(names{subj_name_ind}).sigs(seiz_ind).sigBefore(:,1:19);
                                
                [M,N]=size(X);

                X_n = [];
                if det_norm_eeg=='y'
                    for i=1:N
                        Xtmp=AR_filter(X,i,pfilter);
                        X_n(:,i) = (Xtmp - mean(Xtmp))/std(Xtmp);
                    end
                    X = X_n';
                else
                    X = X';
                end
                S = UX_pre_post * X;
                Q=size(S,1);
                [pAIC,pBIC,aic,mdl] = mos_idMVAR(S,pmax,0);
               %{
               % plot
                figure(1); %first nS EEGs
                for is=1:Q
                    subplot(Q,1,is); plot(X(is,:));
                    xlim([1 N]);ylabel(['X' int2str(is)]);
                end

                figure(2); % all nS sources
                for is=1:Q
                    subplot(Q,2,2*(is-1)+1); plot(S(is,:));
                    xlim([1 N]);ylabel(['S' int2str(is)])
                end
                subplot(Q,2,[2 4 6 8 10 12 14 16 18])
                plot(aic,'b.-'); hold on; plot(mdl,'r.-');
                plot(pAIC,aic(pAIC),'bo'); plot(pBIC,mdl(pBIC),'ro');
                legend(['p_{AIC}=' int2str(pAIC)],['p_{BIC}=' int2str(pBIC)]);
                title(['window N ' int2str(wind_ind)])
                %}
                % connectivity based on state space models
                % model order selection
                switch pcrit
                    case 'a', p=pAIC;
                    case 'b', p=pBIC;
                    case 'c', p=pimp;
                end
                %%% detrend and normalization
                if det_norm_sources=='y'
                    for q=1:Q
                        Stmp=AR_filter(S',q,pfilter);
                %         Sn(:,q)=Stmp;
                        Sn(q,:)=(Stmp-mean(Stmp))/std(Stmp);
                    end
                    %Sn=Sn';
                else
                    Sn=S;
                end

                %% partial variances (correlations)
                [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                qlags=10;
                Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                for jj=1:Q
                    for ii=1:Q
                        if ii~=jj
                            ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                            Hj(jj)=ret.Hy;
                            Hj_j(jj)=ret.Hy_y;
                            Hj_all(jj)=ret.Hy_yzx;
                            Hj_jk(jj,ii)=ret.Hy_yz;
                            % t-test
                            Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                            p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                        end
                    end
                end
                Sj=Hj-Hj_j;
                Tj=Hj_j-Hj_all;
                Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                psign=p_val<valore_p_value_crit;

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p = 15;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=10;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind),' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end
                    
                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    p=5;
                    disp(['for ', names{subj_name_ind}, ' seizure ', int2str(seiz_ind), ' p = ', int2str(p)])
                    [eAm,eSu,~,~]=idMVAR(Sn,p,0);
                    qlags=10;
                    Hj=nan*ones(Q,1); Hj_j=nan*ones(Q,1);
                    Hj_all=nan*ones(Q,1);Hj_jk=nan*ones(Q,Q); p_val=nan*ones(Q,Q);
                    for jj=1:Q
                        for ii=1:Q
                            if ii~=jj
                                ret = its_CElinVAR(eAm,eSu,jj,ii,qlags);
                                Hj(jj)=ret.Hy;
                                Hj_j(jj)=ret.Hy_y;
                                Hj_all(jj)=ret.Hy_yzx;
                                Hj_jk(jj,ii)=ret.Hy_yz;
                                % t-test
                                Var_u=ret.Sy_yzx; Var_r=ret.Sy_yz; Nu=p*Q; Nr=p*(Q-1);
                                p_val(jj,ii) = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,M);
                            end
                        end
                    end
                    Sj=Hj-Hj_j;
                    Tj=Hj_j-Hj_all;
                    Tij_k=Hj_jk-repmat(Hj_all,1,Q);
                    psign=p_val<valore_p_value_crit;      
                end

                if ~isreal(Sj) || ~isreal(Tj) || ~isreal(Tij_k)
                    disp('ALLO, some shit is still complex even after decreasing MVAR order to 5.')
                end
                
                results_pre(seiz_ind).Sj = Sj;
                results_pre(seiz_ind).Tj = Tj;
                results_pre(seiz_ind).Tij_k = Tij_k;
                results_pre(seiz_ind).p_val = p_val;
                results_pre(seiz_ind).psign = psign;  
                results_pre(seiz_ind).pBIC = pBIC;   
                results_pre(seiz_ind).meanSj = mean(Sj);   
                results_pre(seiz_ind).meanTj = mean(Tj);   
                results_pre(seiz_ind).meanTij_k = nanmean(Tij_k(:));   
                results_pre(seiz_ind).Results_doc = sum(sum(p_val<valore_p_value_crit))            ;
            end
        end
        
    end
    
    results_data.(names{subj_name_ind}).results_post = results_post;
    results_data.(names{subj_name_ind}).results_pre = results_pre;
    disp([names{subj_name_ind} " is done"])
end


save('mvar_source_results_MATLAB_paper_pre_post.mat','results_data', '-v7.3')
