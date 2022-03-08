function [column] = stackcols(arrayIn)
% takes array and stacks all the columns on top of each other to make one
% vertical column
% makes it easier to run regressions
column = reshape(arrayIn,[],1);
end