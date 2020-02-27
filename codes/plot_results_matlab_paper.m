clc
clear all
close all
addpath(pwd)
addpath('G:\Works\Ivan\MATLAB CSPVARICA and MVARICA\Functions_MVAR');
addpath('G:\Works\Ivan\MATLAB CSPVARICA and MVARICA\Functions.2');
addpath('G:\Works\Ivan\MATLAB CSPVARICA and MVARICA\NotBoxPlot');

load('mvar_source_results_MATLAB_paper_base_pre.mat')
base_pre = results_data;
load('mvar_source_results_MATLAB_paper_base_post.mat')
base_post = results_data;
load('mvar_source_results_MATLAB_paper_pre_post.mat')
pre_post = results_data;

%load('mvar_source_results_MATLAB_paper.mat')

%addpath('/Users/ivan_kotiuchyi/Documents/Epilepsy_Italy/05_06_windows_information_theory/functions_');
%addpath('/Users/ivan_kotiuchyi/Documents/Epilepsy_Italy/raacampbell-notBoxPlot-3ce29db/code');

names = fieldnames(results_data);

meanSj_abstract_base = [];
meanTj_abstract_base = [];
meanTij_k_abstract_base = [];
mean_significant_links_base = [];

meanSj_abstract_post = [];
meanTj_abstract_post = [];
meanTij_k_abstract_post = [];
mean_significant_links_post = [];

counter_base = 0;
counter_post = 0;

for subj_name_ind = 1:numel(names)
    %results_base
    for seiz_ind = 1:length(base_pre.(names{subj_name_ind}).results_base)
        if isfield(base_pre.(names{subj_name_ind}).results_base, 'meanSj') && ...
                ~isempty(base_pre.(names{subj_name_ind}).results_base(seiz_ind).meanSj)
            counter_base = counter_base + 1;
            meanSj_abstract_base(counter_base, :) = base_pre.(names{subj_name_ind}).results_base(seiz_ind).meanSj;
            meanTj_abstract_base(counter_base, :) = base_pre.(names{subj_name_ind}).results_base(seiz_ind).meanTj;
            meanTij_k_abstract_base(counter_base, :) = base_pre.(names{subj_name_ind}).results_base(seiz_ind).meanTij_k;
            mean_significant_links_base(counter_base, :) = base_pre.(names{subj_name_ind}).results_base(seiz_ind).Results_doc;
        end
    end

    %results_pre
    for seiz_ind = 1:length(base_pre.(names{subj_name_ind}).results_pre)
        if isfield(base_pre.(names{subj_name_ind}).results_pre, 'meanSj') && ...
                ~isempty(base_pre.(names{subj_name_ind}).results_pre(seiz_ind).meanSj)
            counter_post = counter_post + 1;
            meanSj_abstract_post(counter_post, :) = base_pre.(names{subj_name_ind}).results_pre(seiz_ind).meanSj;
            meanTj_abstract_post(counter_post, :) = base_pre.(names{subj_name_ind}).results_pre(seiz_ind).meanTj;
            meanTij_k_abstract_post(counter_post, :) = base_pre.(names{subj_name_ind}).results_pre(seiz_ind).meanTij_k;
            mean_significant_links_post(counter_post, :) = base_pre.(names{subj_name_ind}).results_pre(seiz_ind).Results_doc;
        end
    end

end

figure(1)
subplot(3,4,1);
h=notBoxPlot([meanSj_abstract_base; meanSj_abstract_post], [1.* ones(length(meanSj_abstract_base),1); 2.* ones(length(meanSj_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
set(h(2).data,'markerfacecolor','c')

xticklabels({'base','pre'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')

box on; grid on; 
[pval_base_pre,~]=ranksum(meanSj_abstract_base, meanSj_abstract_post);
title(['p_{base\rightarrowpre}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(meanSj_abstract_base)+median(meanSj_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')

subplot(3,4,2);
h=notBoxPlot([meanTj_abstract_base; meanTj_abstract_post], [1.* ones(length(meanTj_abstract_base),1); 2.* ones(length(meanTj_abstract_post),1)],'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','pre'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')

box on; grid on; %
[pval_base_pre,~]=ranksum(meanTj_abstract_base, meanTj_abstract_post);
title(['p_{base\rightarrowpre}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(meanTj_abstract_base)+median(meanTj_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')

subplot(3,4,3);
h=notBoxPlot([meanTij_k_abstract_base; meanTij_k_abstract_post], [1.* ones(length(meanTij_k_abstract_base),1); 2.* ones(length(meanTij_k_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','pre'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')
box on; grid on; %title('meanTij_k new base / pre');
[pval_base_pre,~]=ranksum(meanTij_k_abstract_base, meanTij_k_abstract_post);
title(['p_{base\rightarrowpre}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(meanTij_k_abstract_base)+median(meanTij_k_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')

subplot(3,4,4);
h=notBoxPlot([mean_significant_links_base; mean_significant_links_post], [1.* ones(length(mean_significant_links_base),1); 2.* ones(length(mean_significant_links_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','pre'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')
box on; grid on; %title('meanTij_k new base / pre');
[pval_base_pre,~]=ranksum(mean_significant_links_base, mean_significant_links_post);
title(['p_{base\rightarrowpre}=' num2str(pval_base_pre)], 'fontsize',10);

%%
meanSj_abstract_base = [];
meanTj_abstract_base = [];
meanTij_k_abstract_base = [];
mean_significant_links_base = [];

meanSj_abstract_post = [];
meanTj_abstract_post = [];
meanTij_k_abstract_post = [];
mean_significant_links_post = [];


counter_base = 0;
counter_post = 0;

for subj_name_ind = 1:numel(names)
    %results_pre
    for seiz_ind = 1:length(pre_post.(names{subj_name_ind}).results_pre)
        if isfield(pre_post.(names{subj_name_ind}).results_pre, 'meanSj') && ...
                ~isempty(pre_post.(names{subj_name_ind}).results_pre(seiz_ind).meanSj)
            counter_base = counter_base + 1;
            meanSj_abstract_base(counter_base, :) = pre_post.(names{subj_name_ind}).results_pre(seiz_ind).meanSj;
            meanTj_abstract_base(counter_base, :) = pre_post.(names{subj_name_ind}).results_pre(seiz_ind).meanTj;
            meanTij_k_abstract_base(counter_base, :) = pre_post.(names{subj_name_ind}).results_pre(seiz_ind).meanTij_k;
            mean_significant_links_base(counter_base, :) = pre_post.(names{subj_name_ind}).results_pre(seiz_ind).Results_doc;

        end
    end
  
    %results_post
    for seiz_ind = 1:length(pre_post.(names{subj_name_ind}).results_post)
        if isfield(pre_post.(names{subj_name_ind}).results_post, 'meanSj') && ...
                ~isempty(pre_post.(names{subj_name_ind}).results_post(seiz_ind).meanSj)
            counter_post = counter_post + 1;
            meanSj_abstract_post(counter_post, :) = pre_post.(names{subj_name_ind}).results_post(seiz_ind).meanSj;
            meanTj_abstract_post(counter_post, :) = pre_post.(names{subj_name_ind}).results_post(seiz_ind).meanTj;
            meanTij_k_abstract_post(counter_post, :) = pre_post.(names{subj_name_ind}).results_post(seiz_ind).meanTij_k;
            mean_significant_links_post(counter_post, :) = pre_post.(names{subj_name_ind}).results_post(seiz_ind).Results_doc;

        end
    end

end

subplot(3,4,5);
h=notBoxPlot([real(meanSj_abstract_base); meanSj_abstract_post], [1.* ones(length(meanSj_abstract_base),1); 2.* ones(length(meanSj_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'pre','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
box on; grid on; %title('meanSj new pre / post');
[pval_base_pre,~]=ranksum(real(meanSj_abstract_base), meanSj_abstract_post);
title(['p_{pre\rightarrowpost}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(real(meanSj_abstract_base))+median(meanSj_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')


subplot(3,4,6);
h=notBoxPlot([meanTj_abstract_base; meanTj_abstract_post], [1.* ones(length(meanTj_abstract_base),1); 2.* ones(length(meanTj_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'pre','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')

box on; grid on; %title('meanTj new pre / post');
[pval_base_pre,~]=ranksum(meanTj_abstract_base, meanTj_abstract_post);
title(['p_{pre\rightarrowpost}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(meanTj_abstract_base)+median(meanTj_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')

subplot(3,4,7);
h=notBoxPlot([meanTij_k_abstract_base; meanTij_k_abstract_post], [1.* ones(length(meanTij_k_abstract_base),1); 2.* ones(length(meanTij_k_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'pre','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
box on; grid on; %title('meanTij_k new pre / post');
[pval_base_pre,~]=ranksum(meanTij_k_abstract_base, meanTij_k_abstract_post);
title(['p_{pre\rightarrowpost}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(meanTij_k_abstract_base)+median(meanTij_k_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')
subplot(3,4,8);
h=notBoxPlot([mean_significant_links_base; mean_significant_links_post], [1.* ones(length(mean_significant_links_base),1); 2.* ones(length(mean_significant_links_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'pre','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')
box on; grid on; %title('meanTij_k new base / pre');
[pval_base_pre,~]=ranksum(mean_significant_links_base, mean_significant_links_post);
title(['p_{pre\rightarrowpost}=' num2str(pval_base_pre)], 'fontsize',10);


%%
meanSj_abstract_base = [];
meanTj_abstract_base = [];
meanTij_k_abstract_base = [];
mean_significant_links_base = [];

meanSj_abstract_post = [];
meanTj_abstract_post = [];
meanTij_k_abstract_post = [];
mean_significant_links_post = [];

counter_base = 0;
counter_pre = 0;
counter_post = 0;

for subj_name_ind = 1:numel(names)
    %results_base
    for seiz_ind = 1:length(base_post.(names{subj_name_ind}).results_base)
        if isfield(base_post.(names{subj_name_ind}).results_base, 'meanSj') && ...
                ~isempty(base_post.(names{subj_name_ind}).results_base(seiz_ind).meanSj)
            counter_base = counter_base + 1;
            meanSj_abstract_base(counter_base, :) = base_post.(names{subj_name_ind}).results_base(seiz_ind).meanSj;
            meanTj_abstract_base(counter_base, :) = base_post.(names{subj_name_ind}).results_base(seiz_ind).meanTj;
            meanTij_k_abstract_base(counter_base, :) = base_post.(names{subj_name_ind}).results_base(seiz_ind).meanTij_k;
            mean_significant_links_base(counter_base, :) = base_post.(names{subj_name_ind}).results_base(seiz_ind).Results_doc;

        end
    end

    %results_post
    for seiz_ind = 1:length(base_post.(names{subj_name_ind}).results_post)
        if isfield(base_post.(names{subj_name_ind}).results_post, 'meanSj') && ...
                ~isempty(base_post.(names{subj_name_ind}).results_post(seiz_ind).meanSj)
            counter_post = counter_post + 1;
            meanSj_abstract_post(counter_post, :) = base_post.(names{subj_name_ind}).results_post(seiz_ind).meanSj;
            meanTj_abstract_post(counter_post, :) = base_post.(names{subj_name_ind}).results_post(seiz_ind).meanTj;
            meanTij_k_abstract_post(counter_post, :) = base_post.(names{subj_name_ind}).results_post(seiz_ind).meanTij_k;
            mean_significant_links_post(counter_post, :) = base_post.(names{subj_name_ind}).results_post(seiz_ind).Results_doc;
        end
    end

end

subplot(3,4,9);
h=notBoxPlot([meanSj_abstract_base; meanSj_abstract_post], [1.* ones(length(meanSj_abstract_base),1); 2.* ones(length(meanSj_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)

box on; grid on; %title('meanSj new base / post');
xlabel('Information storage', 'fontsize',10,'FontWeight','bold')
[pval_base_pre,~]=ranksum(meanSj_abstract_base, meanSj_abstract_post);
title(['p_{base\rightarrowpost}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%text(1.5,0.5*(median(meanSj_abstract_base)+median(meanSj_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')


subplot(3,4,10);
h=notBoxPlot([meanTj_abstract_base; meanTj_abstract_post], [1.* ones(length(meanTj_abstract_base),1); 2.* ones(length(meanTj_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
xlabel('Information transfer', 'fontsize',10,'FontWeight','bold')
box on; grid on; %title('meanTj new base / post');
[pval_base_pre,~]=ranksum(meanTj_abstract_base, meanTj_abstract_post);
title(['p_{base\rightarrowpost}=' num2str(pval_base_pre)],'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%%text(1.5,0.5*(median(meanTj_abstract_base)+median(meanTj_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')

subplot(3,4,11);
h=notBoxPlot([meanTij_k_abstract_base; meanTij_k_abstract_post], [1.* ones(length(meanTij_k_abstract_base),1); 2.* ones(length(meanTij_k_abstract_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
xlabel('Conditional information transfer', 'fontsize',10,'FontWeight','bold')
box on; grid on; %title('meanTij_k new base / post');
[pval_base_pre,~]=ranksum(meanTij_k_abstract_base, meanTij_k_abstract_post);
title(['p_{base\rightarrowpost}=' num2str(pval_base_pre)], 'fontsize',10);
%[pval_base_post,hval_base_post]=ranksum(meanTj_for_all_base, meanTj_for_all_pre(:,1));
%[pval_pre_post,hval_pre_post]=ranksum(meanTj_for_all_pre,meanTj_for_all_pre(:,1));
%%text(1.5,0.5*(median(meanTij_k_abstract_base)+median(meanTij_k_abstract_post)),['p_{base pre}=' num2str(pval_base_pre)],'HorizontalAlignment','center')
%text(2.5,0.5*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{pre post}=' num2str(pval_pre_post)],'HorizontalAlignment','center')
%text(2,0.75*(median(meanTj_for_all_pre)+median(meanTj_for_all_pre(:,1))),['p_{base post}=' num2str(pval_base_post)],'HorizontalAlignment','center')
subplot(3,4,12);
h=notBoxPlot([mean_significant_links_base; mean_significant_links_post], [1.* ones(length(mean_significant_links_base),1); 2.* ones(length(mean_significant_links_post),1)], 'markMedian',true,'jitter',0.6,'style', 'patch'); 
xticklabels({'base','post'})
a = get(gca,'xticklabels');  
set(gca,'xticklabels',a,'fontsize',10,'FontWeight','bold')
d = [h.data];
set(d(1:end), 'markerfacecolor', [0.75,0.75,0.75], 'color', [0,0,0]);
set(d, 'markersize', 6);
xtickangle(45)
%xlabel('Conditional information transfer', 'FontSize',16,'FontWeight','bold')
box on; grid on; %title('meanTij_k new base / pre');
[pval_base_post,~]=ranksum(mean_significant_links_base, mean_significant_links_post);
xlabel('Significant links', 'fontsize',10,'FontWeight','bold')
title(['p_{base\rightarrowpost}=' num2str(pval_base_post)], 'fontsize',10);