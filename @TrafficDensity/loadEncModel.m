function obj = loadEncModel(obj)
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Get encounter model speed information by altitude (for computing the event rate)

uncor = em_read(obj.filenames.uncorEncModel); 
uncor.h_cutpoints = [500,1200,3000,5000,18000]; % From uncorrelated model project report ATC-404
uncor.maxagl = 4; % Index of maximum AGL
cor = em_read(obj.filenames.corEncModel);
cor.h_cutpoints = [500,1000,3000,5000,10000,18000,29000,40000,inf]; % From correlated model project report ATC-440
cor.maxagl = 4; % Index of maximum AGL

% Process uncorrelated model (for 1200-code traffic)
[uncor.p_speed_alt,uncor.speedvar] = processEncModelSpeed(uncor,'"v"');

% Process correlated model (for discrete code traffic)
[cor.p_speed_alt,cor.speedvar] = processEncModelSpeed(cor,'"v_1"');
cor.p_speed_alt = (cor.p_speed_alt+processEncModelSpeed(cor,'"v_2"'))/2;

% Save in object
obj.uncor = uncor;
obj.cor = cor;
end

% Get the speed distribution as a function of altitude
function [p_speed_alt,speedvar] = processEncModelSpeed(model,speedvarname)
altvar = find(strcmp(model.labels_initial,'"L"'));      % Altitude index
speedvar = find(strcmp(model.labels_initial,speedvarname));    % Speed index
parents = model.G_initial(:,speedvar);                  % Parents of variable in Bayesian network
nalts = model.r_initial(altvar);                        % Number of altitude bins
nspeeds = model.r_initial(speedvar);                    % Number of speed bins
nparents = size(model.N_initial{speedvar},2);          % Number of parental instantiations
x = zeros(sum(parents),nparents); % Get parent matrix (row is each parent, column is the parental instantiation)
for pp = 1:nparents
    x(:,pp) = aind2sub(model.r_initial(parents),pp); % x is an array of parent instantiations
end
p_speed_alt = zeros(nspeeds,nalts);

% Get probability of speed given altitude layers
for aa = 1:nalts
    p_speed_alt(:,aa) = sum(model.N_initial{speedvar}(:,x(find(parents)==altvar,:) == aa),2);
    p_speed_alt(:,aa) = p_speed_alt(:,aa)./sum(p_speed_alt(:,aa));
end
end

function parameters = em_read(filename)
% EM_READ  Reads an encounter model parameters file.
% Reads an encounter model parameters file and returns a structure
% containing the parsed data.
%
%   EM_READ(FILENAME) reads the parameters contained in the specified file
%   and returns the parameters in a structure. Included in this structure
%   are the following fields:
%        labels_initial
%             n_initial
%             G_initial
%             r_initial
%             N_initial
%     labels_transition
%          n_transition
%          G_transition
%          r_transition
%          N_transition
%            boundaries
%        resample_rates
%          temporal_map
%             zero_bins

f = fopen(filename);
validate_label(f, '# labels_initial');
p.labels_initial = scanline(f, '%s', ',');
p.n_initial = numel(p.labels_initial);
validate_label(f, '# G_initial');
p.G_initial = logical(scanmatrix(f, p.n_initial));
validate_label(f, '# r_initial');
p.r_initial = scanmatrix(f);
validate_label(f, '# N_initial');
dims_initial = getdims(p.G_initial, p.r_initial, 1:p.n_initial);
p.N_initial = array2cells(scanmatrix(f), dims_initial);
validate_label(f, '# labels_transition');
p.labels_transition = scanline(f, '%s', ',');
p.n_transition = numel(p.labels_transition);
validate_label(f, '# G_transition');
p.G_transition = logical(scanmatrix(f, p.n_transition));
validate_label(f, '# r_transition');
p.r_transition = scanmatrix(f);
validate_label(f, '# N_transition');
dims_transition = getdims(p.G_transition, p.r_transition, (p.n_initial+1):p.n_transition);
p.N_transition = array2cells(scanmatrix(f), dims_transition);
validate_label(f, '# boundaries');
p.boundaries = cell(1,p.n_initial);
for i=1:p.n_initial
    p.boundaries{i} = scanmatrix(f);
end
validate_label(f, '# resample_rates');
p.resample_rates = scanmatrix(f);
fclose(f);
p.temporal_map = extract_temporal_map(p.labels_transition);
p.zero_bins = extract_zero_bins(p.boundaries);
parameters = p;
end

%Get bins where variables are set to zero
function zero_bins = extract_zero_bins(boundaries)
zero_bins = cell(1, numel(boundaries));
for i = 1:numel(boundaries)
    b = boundaries{i};
    z = [];
    if numel(b) > 2
        for j = 2:numel(b)
            if b(j - 1) < 0 && b(j) > 0
                z = j - 1;
            end
        end
    end
    zero_bins{i} = z;    
end
end

function temporal_map = extract_temporal_map(labels_transition)
% Does not assume that order of variables at t match that at
% t+1 (order assumption not true for ECEM, but true for previous models)
temporal_map = [];
for i = 1:numel(labels_transition)
    t = findstr(labels_transition{i}, '(t)');
    if ~isempty(t)
        temporal_map = [temporal_map;i,find(contains(labels_transition,[labels_transition{i}(1:t),'t+1)']))]; 
    end
end
end

function x = scanline(fid, typename, delimiter)
a = textscan(fgetl(fid), typename, 'Delimiter', delimiter);
x = a{1};
end

function x = scanmatrix(fid, num_rows)
if nargin < 2
    x = scanline(fid, '%f', ' ');
else
    x = zeros(num_rows);    
    for i = 1:num_rows
        x(i,:) = scanline(fid, '%f', ' ');
    end
end
end

function c = array2cells(x, dims)
c = cell(size(dims, 1), 1);
index = 1;
for i = 1:numel(c)
    c{i} = zeros(dims(i,1), dims(i,2));
    c{i}(:) = x(index:(index-1+numel(c{i})));
    index = index + numel(c{i});
end
end

function dims = getdims(G, r, vars)
n = size(G,1);
dims = zeros(n, 2);
for i = vars
   q = prod(r(G(:,i))); % q_i
   dims(i,:) = [r(i) q];
end
end

function validate_label(fid, s)
t = fgetl(fid);
if ~strcmp(t, s)
    error('Invalid parameters file');
end
end

function x = aind2sub(siz,ndx)
siz = siz(:);
x = zeros(size(siz));
n = length(siz);
k = [1; cumprod(siz(1:end-1))];
for i = n:-1:1
  vi = rem(ndx-1, k(i)) + 1;         
  x(i) = (ndx - vi)/k(i) + 1; 
  ndx = vi;     
end
end