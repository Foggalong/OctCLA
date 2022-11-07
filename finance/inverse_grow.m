% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function newInv = inverse_grow(invA, a, alpha)
    % INVERSE_GROW adjust the inverse if gaining a row and column
    %
    % Takes the inverse of a symmetric matrix, a column vector, and a
    % scalar as inputs and then returns the inverse of the original
    % matrix if its dimensions were both increased by 1, i.e. if
    % newA = [A, a; a, alpha] it takes A^-1, a, and alpha as input and
    % returns newA^-1.
    %
    % See also, INVERSE_SHRINK

    % shortcut variabels
    c = invA*a;
    beta = 1/(alpha - c'*a);
    % preallocate the new inverse
    newInv = zeros(size(invA)+1);
    % update the top left block
    newInv(1:end-1, 1:end-1) = invA + beta*(c*c');
    % update the top right block
    newInv(1:end-1, end) = -beta*c;
    % update the bottom left block
    newInv(end, 1:end-1) = -beta*c';
    % update the bottom right block
    newInv(end, end) = beta;
end
