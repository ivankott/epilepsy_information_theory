function [NumComponents]=csp_numcomponents_rim_onlynum(x1, x2, dist_cov_matr, pl)
    % inputs
    % * x1: EEG signals - condition 1 (n_channels x n_samples x n_trials)
    % * x2: EEG signals - condition 2 (n_channels x n_samples x n_trials)
    % * dist_cov_matr: default is 0.95 which describe 95% of the total 
    %   distance between each of covariance matrices 
    % * pl: plot of contribution to total square Riemannian distance by 
    %   eigenvalues corresponding to each spatial filter

    % outputs
    % * NumComponents: number of CSP components to keep
    
    if nargin < 3 || isempty(dist_cov_matr)
        dist_cov_matr = 0.95;
    end
    if nargin < 4
        pl = true;
    end
    
    % CSP implementation
    % permute
    xn1 = permute(x1, [3 1 2]);
    xn2 = permute(x2, [3 1 2]);
    
    [t1, m1, n1] = size(xn1);
    sigma1 = zeros([m1,m1]);
    
    % obtain mean cov matrices for each trial in each class
    for n_trial=1:t1
        sigma1 = sigma1 + cov(permute(xn1(n_trial, :, :), [3 2 1]))/t1;
    end
    sigma1 = sigma1 / trace(sigma1);

    [t2, m2, n2] = size(xn2);
    sigma2 = zeros([m2,m2]);
    
    for n_trial=1:t2
        sigma2 = sigma2 + cov(permute(xn2(n_trial, :, :), [3 2 1]))/t2;
    end
    sigma2 = sigma2 / trace(sigma2);

    % Solve the eigenvalue problem S1�W = l�S2�W
    % W: mixing matrix (eigenvectors)
    % A: demixing matrix (inverse of W)
    % L: eigenvalues
    
    % I think we should explain in the paper, why we use sum below, 
    % in case we are going to define and understand everything
    [W, L] = eig(sigma1, sigma1 + sigma2, 'vector');
    [eigval_sorted, eigval_indices] = sort(L, 'ascend'); 

    % Riccardo's method:
    % sigma3= sigma1+sigma2;
    % sigma4 = sigma1\sigma3;
    % [W,L] = eigenshuffle(sigma4);
    % [Ws,Ls] = eig(sigma4,'vector');
    % [~,Lindex] = sort(real(Ls),1,'descend');
    % Wseq = Ws(:,Lindex);
    % W = Wseq;
    
    %2-norm normalization:
    for w=1:size(W,2)
        W(:,w)=W(:,w)./norm(W(:,w));
    end 

    %reorder in descendent order
    %[~, eigval_indices] = sort(L, 'descend');
    W = W(:,eigval_indices);
    % ? need to check it out and comment ?
    W = W.*-1;
    A = (inv(W))';  
    
    %[eigval_sorted, eigval_indices] = sort(L');
    
    % details of implementation could be found in paper https://www.dropbox.com/s/1dq7qi1mdc42vq8/Common%20Spatial%20Pattern.pdf?dl=0
    % first, we find the total Riemannian distance between the two class-related mean covariance matrices
    % Riemannian distance between the two class-related mean matrices which is directly linked with the set of CSP eigenvalues λj .
    delta_rim = [];
    for i=1:length(eigval_sorted)
        delta_rim(i) = log(eigval_sorted(i) / (1 - eigval_sorted(i))).^2;
    end
    total_delta_rim = sqrt(sum(delta_rim));
    
    function_delta_rim = abs(delta_rim);
    % then we normalize all fractional distances to get the percentage of total variance 
    % described by each spatial filter which is directly dependent on the corresponding eigenvalue.
    function_delta_rim = function_delta_rim./(total_delta_rim^2);
    
    % then let's find sorted function_delta_rim 
    [dist, ind] = sort(function_delta_rim, 'descend');
    
    % after that we sum up the fractional distances from two sides of sorted function_delta_rim
    % (as it is done in CSP) untill the sum of the fractional disctances reach it's theshold = dist_cov_matr
    dist_cov = zeros([1, length(ind)]);
    % for plot
    eigen_val_corresponding_to_dist = zeros([1,length(ind)]);
    index_of_eig_for_plt = strings(size(ind));
    for i=1:length(ind)
        if sum(dist_cov)<=dist_cov_matr
            dist_cov(i) = dist(i);
            % for plot
            eigen_val_corresponding_to_dist(i) = eigval_sorted(ind(i));
            index_of_eig_for_plt(i) = strcat('\lambda_{', int2str(i) ,'}');
        end
    end
    dist_cov = nonzeros(dist_cov);
    eigen_val_corresponding_to_dist = nonzeros(eigen_val_corresponding_to_dist);
    index_of_eig_for_plt(cellfun('isempty',index_of_eig_for_plt)) = [];
    
    NumComponents = length(dist_cov);
    disp([int2str(NumComponents), ' components (corresponding to same # of eigenvalues) were sufficient for describing 95% of the total distance between mean covariance matrices for each class.'])
    
    while (size(W,2) > NumComponents)
        i_tmp_rmv = floor(size(W,2)/2);
        W(:,i_tmp_rmv) = [];
        A(i_tmp_rmv,:) = [];
    end

    % TODO text positions
    if pl==1
        figure();
        sz = 60;
        scatter(L,function_delta_rim,sz,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
        hold on
        title({'Evolution of squared distance supported by spatial filter and corresponding eigvalue','\newline','(black dots - eigenvalues of selected filters in W)'},'FontSize',14);
        xlim([min(L)-0.01, max(L)+0.01]);
        xlabel('\lambda_{i}, eigenvalue','FontSize',16);
    	ylabel('Contribution to total squared distance','FontSize',16);

        scatter(eigen_val_corresponding_to_dist,dist_cov,sz,'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[0 0 0])
        hold on
        for nn=1:length(index_of_eig_for_plt)
            if mod(nn,2)==0
                text(eigen_val_corresponding_to_dist(nn)- 0.01,dist_cov(nn)- 0.01,index_of_eig_for_plt(nn),...
                'fontsize',14,...
                'verticalalignment','bottom');
            else
                text(eigen_val_corresponding_to_dist(nn)+ 0.01,dist_cov(nn)+ 0.01,index_of_eig_for_plt(nn),...
                'fontsize',14,...
                'verticalalignment','bottom');
            end
        end
        grid on
        hold off
    end
end