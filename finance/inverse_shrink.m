% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function newInv = inverse_shrink(invA, i)
    % INVERSE_SHRINK adjust the inverse if removing row and column i
    %
    % Takes the inverse of a symmetric matrix and an index i as input
    % and then returns the inverse of the original matrix if its
    % dimensions were both decreased by 1 through removal of row and
    % column i. That is, if after permuting the ith columns and rows
    % to the end, A = [newA, a; a, alpha] then the function takes A^-1
    % (and i) as input and returns newA^-1.
    %
    % See also, INVERSE_GROW

    % permute rows and cols of invA so that ith row and col are last
    invAperm = invA([1:(i-1) (i+1):end i], [1:(i-1) (i+1):end i]);
    % shortcut variables
    B = invAperm(1:end-1, 1:end-1);
    b = invAperm(1:end-1, end);
    beta = invAperm(end, end);
    % calculate new inverse
    newInv = B - (b*b')/beta;
end
