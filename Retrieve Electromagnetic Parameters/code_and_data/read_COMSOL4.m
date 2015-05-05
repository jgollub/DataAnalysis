% Read frequency and S parameters in text file exported by COMSOL 4.
% 
% This function takes 1 input parameter:
%   filename     a string with the name of the text file.
% It returns 3 output parameters:
%   f            a column vector with the frequencies;
%   S11          a column vector with the reflection;
%   S21          a column vector with the transmission.
% 
% To export the data in Comsol, define derived values for S11 and S21,
% then calculate their values in a single table, and export the table
% to a text file.
% 
% This function takes the conjugate of S11 and S21 before returning
% them to use the same convention than the retrieve function.
% 
% Copyright (c) 2013,2014 Stéphane Larouche

function [f, S11, S21] = read_COMSOL4(filename)

data = textscan(fopen(filename), '%f %f %f', 'Headerlines', 5);

f = data{1};
S11 = conj(data{2});
S21 = conj(data{3});
