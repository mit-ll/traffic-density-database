% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Script will build accumarray for @max function: this is over an order of
% magnitude faster than calling accumarray directly with @max

if verLessThan('matlab','9.6')
    warning('Can only generate mex for Matlab accumarray in Matlab Version 2019a (9.6) or newer');
    return;
end

if ~license('test','matlab_coder')
    warning('Can only generate mex for Matlab accumarray if have license for Matlab Coder');
end

%% Build
vectorType1 = coder.typeof(uint32(1), [inf 1], [true false]);
vectorType2 = coder.typeof(uint32(1), [inf 1], [true false]);

codegen accumarraymax -args {vectorType2 ,vectorType1, double(0)}