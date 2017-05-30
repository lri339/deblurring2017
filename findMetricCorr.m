function [ output_args ] = findMetricCorr( input_args )
%FINDMETRICCORR 
% LI
% Not used - see R code
% Code to find correlations between features metric etc

% Note: need to edit this to suppress a lot of output, unnecessary
% calculation and variables
% script from 2016 CVPR review study
demo_bt_ranking;


% This is the score of metric which we are testing
%item 12 and 13 (whioch are missing in score matrix)
%load('score1_286_12_13_missing.mat')
%score = score(1:286)
%score = reshape(score, [13,22]);
load('finalscore.mat');
score_reshaped = reshape(score, [14,100]);
%score2 = score2';
% score = score';
% bt_reshaped = reshape(score_matrix, [1400,1]);
% 
% 
% notnan_indices = find(~(isnan(bt_reshaped)|isnan(score)));
% score_corr = corr(score(notnan_indices), bt_reshaped(notnan_indices), 'type', 'Spearman');
% 
% 
% for i=1:100
%     tmp = score_matrix(i,:)';
%     notnan_indices = ~(isnan(tmp) | isnan(score_reshaped(:,i)));
%     corr_byimage(i) = corr(score_matrix(i,notnan_indices)', score_reshaped(notnan_indices,i), 'type', 'Spearman')
% end
% mean(corr_byimage)


[ meanCorr, corrs ] = getCorrs( score_matrix', score_reshaped, 'Spearman' )

load('finaldetails.mat');

[~,lengthDetails] = size(details);
for i=1:lengthDetails
    if isnan(score(i))
    sparsity(i) = NaN;
    smallgrad(i) = NaN;
    metric_q(i) = NaN;
    auto_corr(i) = NaN;
    norm_sps(i) = NaN;
    cpbd(i) = NaN;
    pyr_ring(i) = NaN;
    saturation(i) = NaN;
    else
    sparsity(i) = details(i).sparsity;
    smallgrad(i) = details(i).smallgrad;
    metric_q(i) = details(i).metric_q;
    auto_corr(i) = details(i).auto_corr;
    norm_sps(i) = details(i).norm_sps;
    cpbd(i) = details(i).cpbd;
    pyr_ring(i) = details(i).pyr_ring;
    saturation(i) = details(i).saturation;
    end
end

[ meanCorr(1), corrs ] = getCorrs( score_matrix', score_reshaped, 'Spearman' )
[ meanCorr(2), corrs ] = getCorrs( score_matrix', reshape(sparsity',[14,100]) , 'Spearman' )
[ meanCorr(3), corrs ] = getCorrs( score_matrix', reshape(smallgrad',[14,100]) , 'Spearman' )
[ meanCorr(4), corrs ] = getCorrs( score_matrix', reshape(metric_q',[14,100]) , 'Spearman' )
[ meanCorr(5), corrs ] = getCorrs( score_matrix', reshape(auto_corr',[14,100]) , 'Spearman' )
[ meanCorr(6), corrs ] = getCorrs( score_matrix', reshape(norm_sps',[14,100]) , 'Spearman' )
[ meanCorr(7), corrs ] = getCorrs( score_matrix', reshape(cpbd',[14,100]) , 'Spearman' )
[ meanCorr(8), corrs ] = getCorrs( score_matrix', reshape(pyr_ring',[14,100]) , 'Spearman' )
[ meanCorr(9), corrs ] = getCorrs( score_matrix', reshape(saturation',[14,100]) , 'Spearman' )
darkchannelscore = runDarkChannelonImages()
[ meanCorr(10), corrs ] = getCorrs( score_matrix', reshape(darkchannelscore',[14,100]) , 'Spearman' )

[ meanCorr(1), corrs ] = getCorrs( score_matrix', score_reshaped, 'Kendall' )
[ meanCorr(2), corrs ] = getCorrs( score_matrix', reshape(sparsity',[14,100]) , 'Kendall' )
[ meanCorr(3), corrs ] = getCorrs( score_matrix', reshape(smallgrad',[14,100]) , 'Kendall' )
[ meanCorr(4), corrs ] = getCorrs( score_matrix', reshape(metric_q',[14,100]) , 'Kendall' )
[ meanCorr(5), corrs ] = getCorrs( score_matrix', reshape(auto_corr',[14,100]) , 'Kendall' )
[ meanCorr(6), corrs ] = getCorrs( score_matrix', reshape(norm_sps',[14,100]) , 'Kendall' )
[ meanCorr(7), corrs ] = getCorrs( score_matrix', reshape(cpbd',[14,100]) , 'Kendall' )
[ meanCorr(8), corrs ] = getCorrs( score_matrix', reshape(pyr_ring',[14,100]) , 'Kendall' )
[ meanCorr(9), corrs ] = getCorrs( score_matrix', reshape(saturation',[14,100]) , 'Kendall' )
darkchannelscore = runDarkChannelonImages()
[ meanCorr(10), corrs ] = getCorrs( score_matrix', reshape(saturation',[14,100]) , 'Kendall' )
% this is the BT score from the human survey
%blurry original image not included (for now)
%tmp = score_matrix(1:22, 2:14)


%NOT ALL HAVE 13 IMAGES THOOUGH
%set it as default and deal with NaN need to process slightly different
%for i=1:22
%    score_corr(i) = corr(tmp(i,:)', score(:,i), 'type', 'Spearman')
%end

%mean(score_corr)

%overall corr..
tmp1 = reshape(score_matrix, [1400,1]);
tmpindices = find(~(isnan(score) | isnan(tmp1)));
corr(score(tmpindices), tmp1(tmpindices), 'type', 'Spearman')


end

