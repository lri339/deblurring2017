%%% LI
% save all features to excel

filename = 'features.xlsx';
%get LR score
load('finalscore.mat');
%get LR features
load('finaldetails.mat');
%get simple contrast feature
load('score_contrastmetric.mat');
%get simple ringing feature
load('score_ringingsimple.mat');
%get ringing feature
load('score_ringingperceptualV2.mat');
%get darkchannel feature
load('score_for_dark_channel.mat');

for i = 1:length(details)
    
    if isempty(details(i).sparsity)
        sparsity(i) = NaN;
    else
        sparsity(i) = details(i).sparsity;
    end
    
    if isempty(details(i).smallgrad)
        smallgrad(i) = NaN;
    else
        smallgrad(i) = details(i).smallgrad;
    end
    
    if isempty(details(i).metric_q)
        metric_q(i) = NaN;
    else
        metric_q(i) = details(i).metric_q;
    end
    
    if isempty(details(i).auto_corr)
        auto_corr(i) = NaN;
    else
        auto_corr(i) = details(i).auto_corr;
    end
    
    if isempty(details(i).norm_sps)
        norm_sps(i) = NaN;
    else
        norm_sps(i) = details(i).norm_sps;
    end

    if isempty(details(i).cpbd)
        cpbd(i) = NaN;
    else
        cpbd(i) = details(i).cpbd;
    end
    
    if isempty(details(i).pyr_ring)
        pyr_ring(i) = NaN;
    else
        pyr_ring(i) = details(i).pyr_ring;
    end
        
    if isempty(details(i).saturation)
        saturation(i) = NaN;
    else
        saturation(i) = details(i).saturation;
    end
        

end

humanScore = reshape(score_matrix', [1400, 1]);
M = [score', sparsity', smallgrad', metric_q', auto_corr', norm_sps', cpbd', darkchannelscore', saturation',  score_contrastmetric', score_ringingsimple', pyr_ring', score_perceptualringing', humanScore]

xlswrite(filename,M)