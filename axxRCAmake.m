function axxRCA = axxRCAmake(path_names, W, condsToUse, type, axxRCA)
    % Make axxRCA struct that includes the axx Wave data and the W (Ch x
    % component) matrix from fullRCA
    % INPUTS:
    %   pathnames (required): cell vector of string directory names housing Axx_c00x.mat
    %   condsToUse: vector of conditions to use
    % 
    % OUTPUTS:
    %   Struct with the following fields
    %   Wave: Wave data from Axx_c00x.mat
    %   W: linear transformation matrix to go from sensor space to RC-space
    
    if nargin < 3
        condsToUse = 0;
    else
    end
   
    if nargin < 4
        type = '';
    else
        type = sprintf('%s_', type);
    end
    
    if nargin < 5
        axxRCA.Wave = {};
    end
    for f = 1:length(path_names)
        if condsToUse == 0
            axx_matfiles = subfiles(sprintf('%s/%sExp_MATL_HCN_128_Avg/Axx_c*.mat',fileparts(path_names{f}), type), 1);
            axx_matfiles = axx_matfiles(cell2mat(cellfun(@(x) ~contains(x, 'trials'), axx_matfiles, 'uni', false)));
            condsToUse = cellfun(@(x) str2num(x(end-4)), axx_matfiles);
        else
        end
        axx_matfiles = arrayfun(@(x) ... 
                sprintf('%s/%sExp_MATL_HCN_128_Avg/Axx_c%03d.mat',fileparts(path_names{f}),type,x),condsToUse,'uni',false);
        tmpStrct = cellfun(@(x) load(x),axx_matfiles,'uni',false);
        readyStrct = cellfun(@(x) mrC.axx(x),tmpStrct,'uni',false);
        axxRCA.Wave(:,f) = cellfun(@(x) x.Wave, readyStrct,'uni',false);
    end

    nanDims = [1,2];
    axxRCA.Wave = cellfun(@(x) Zero2NaN(x,nanDims),axxRCA.Wave,'uni',false);
    axxRCA.Projected = rcaProject(axxRCA.Wave, W);
    axxRCA.Projected = cellfun(@(x) Zero2NaN(x,nanDims),axxRCA.Projected,'uni',false);
end


