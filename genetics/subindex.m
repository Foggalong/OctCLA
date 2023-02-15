% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% HACK Sometimes we will use S and D to index the full set of assets
% and sometimes we will only be indexing the free or bounded assets. A
% hack to let us do that now is this function which lets us convert
% from one to the other.

function output = subindex(X, Y)
    % SUBINDEX returns locations of elements of X appearing in Y
    %
    % Given index vectors X and Y of some other vector V, return an
    % index vector of Y giving locations of elements of X that appear
    % in Y.
    %
    % See also, ISMEMBER.

    [bool_XinY, ind_XinY] = ismember(X, Y);
    output = ind_XinY(bool_XinY);
end
