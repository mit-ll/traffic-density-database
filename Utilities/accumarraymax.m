function o = accumarraymax(a,b,c) %#codegen
% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Caller function to be compiled using Matlab Coder
    o = accumarray(a,b,[c,1],@max);
end

