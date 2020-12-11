function [histx,histy,sumyw] = hist_w1d(data,varargin)
% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: this function will give the weighted distribution given the
% inputs; if no outputs, will plot the values with optional labelx and
% labely

% Parse optional inputs
opts = inputParser;

opts.addOptional('w',ones(size(data)),@(x)all([isnumeric(x),size(x)==size(data)]));
opts.addOptional('minval', min(data), @isnumeric);  % Minimum value of data
opts.addOptional('maxval', max(data), @isnumeric);  % Maximum value of data
opts.addOptional('bins', 100, @isnumeric);          % Number of bins

opts.addOptional('smooth',false,@islogical);        % Boolean to smooth final density (using normal kernel)
opts.addOptional('kernelSigma',(4/(3*length(data)))^(1/5)*nanstd(data),@isnumeric); % 'optimal' smoothing for normal distribution (as in Matlab ksdensity and "Bowman and Azzalini")
opts.addOptional('labelx','',@ischar);              % Plotting options
opts.addOptional('labely','',@ischar);
opts.addOptional('linestyle','k-',@ischar);
opts.addOptional('plotresults',~logical(nargout),@islogical); % Plot results? will plot if no output specified

opts.parse(varargin{:});

% Process data
bins = linspace(opts.Results.minval,opts.Results.maxval,opts.Results.bins+1)'; % Get bins
histx = bins(1:end-1)+diff(bins(1:2))/2;

[histy_uw,histy_bins] = histc(data,bins);  %#ok<HISTC>
histy_uw = histy_uw(1:end-1); % Unweighted density

highb = histy_bins>length(histy_uw); % Put out of bounds values to maximum value
lowb = histy_bins<1;
histy_bins(highb) = length(histy_uw);
histy_bins(lowb) = 1;

histy_w = accumarray(histy_bins,opts.Results.w,size(histy_uw)); % Properly weight distribution

if opts.Results.smooth % Smooth if requested
    histy_w = smooth(histx,histy_w,opts.Results.kernelSigma);    
    fprintf('Smoothing with %.4f bandwidth\n',opts.Results.kernelSigma);
end

sumyw = sum(histy_w);
histy = histy_w./sumyw; % Normalize density

% Plot data if needed
if opts.Results.plotresults
    plot(histx,histy,opts.Results.linestyle)
    xlabel(opts.Results.labelx);
    ylabel(opts.Results.labely);
    curry = ylim;
    ylim([0,max(curry)]);
    xlim([min(histx),max(histx)]);
end

% From local_smooth
function y = smooth(t,x,sigma)
    y = x;
    
    if sigma == 0
        return
    end

    for i=1:length(t)
        w = normpdf(t, t(i), sigma);
        s = sum(w);       % denominator when normalizing
        s = s + (s == 0); % ensures denominator is not 0
        w = w / s;        % normalizes
        y(i,:) = (x'*w)';    
    end
