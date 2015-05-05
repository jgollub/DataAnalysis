% Show the result of the retrieval procedure.
% 
% This function shows the impedance, the refractive index, the
% permittivity and the permeability extracted using retrieval. It
% allows rapid evaluation and comparison of results. It does NOT
% produce publication quality figures.
% 
% This function takes 6 or 7 input arguments:
%   title      a string containing the title to give to the figure
%              window;
%   f          a vector (or a cell array of vectors) of the
%              frequencies at which the retrieval was done;
%   z          a vector (or a cell array of vectors) of the impedance
%              of the metamaterial(s);
%   n          a vector (or a cell array of vectors) of the refractive
%              index of the metamaterial(s);
%   epsilon    a vector (or a cell array of vectors) of the
%              permittivity of the metamaterial(s);
%   mu         a vector (or a cell array of vectors) of the
%              permeability of the metamaterial(s);
%   name       (optional) a string (or a cell array of strings) with
%              the names of the matematerial(s).
% 
% Examples:
%   show_retrieval('Example', f, z, n, epsilon, mu)
%   show_retrieval('Comparison of 2 designs', {f_1, f_2}, {z_1, z_2}, ...
%                  {n_1, n_2}, {epsilon_1, epsilon_2}, {mu_1, mu_2}, ...
%                  {'Design 1', 'Design 2'})
% 
% Copyright (c) 2008,2009 Stéphane Larouche

function [] = show_retrieval(window_title, f, z, n, epsilon, mu, name)



colors = {'b','r','g','c','m','y','k'};


if ~iscell(f)
	f = {f};
end
if ~iscell(z)
	z = {z};
end
if ~iscell(n)
	n = {n};
end
if ~iscell(epsilon)
	epsilon = {epsilon};
end
if ~iscell(mu)
	mu = {mu};
end
if nargin > 6
	if ~iscell(name)
		name = {name};
	end
else
	name = cell(size(z));
end

fig = figure('Name', window_title);
units = get(fig, 'units');
set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(fig, 'units', units);

subplot(2,2,1)
hold
box on
grid on
title('\it{}z')
xlabel('\it{}f')
legends = {};
for k = 1:length(f)
	plot(f{k}, real(z{k}), strcat(colors{mod(k, length(colors))}, '-'))
	plot(f{k}, imag(z{k}), strcat(colors{mod(k, length(colors))}, '--'))
	legends{end+1} = ['real({\it{}z})' ' ' name{k}];
	legends{end+1} = ['imag({\it{}z})' ' ' name{k}];
end
legend(legends)

subplot(2,2,3)
hold
box on
grid on
title('\it{}n')
xlabel('\it{}f')
legends = {};
for k = 1:length(f)
	plot(f{k}, real(n{k}), strcat(colors{mod(k, length(colors))}, '-'))
	plot(f{k}, imag(n{k}), strcat(colors{mod(k, length(colors))}, '--'))
	legends{end+1} = ['real({\it{}n})', ' ', name{k}];
	legends{end+1} = ['imag({\it{}n})', ' ', name{k}];
end
legend(legends)

subplot(2,2,2)
hold
box on
grid on
title('\it{}\epsilon')
xlabel('\it{}f')
legends = {};
for k = 1:length(f)
	plot(f{k}, real(epsilon{k}), strcat(colors{mod(k, length(colors))}, '-'))
	plot(f{k}, imag(epsilon{k}), strcat(colors{mod(k, length(colors))}, '--'))
	legends{end+1} = ['real({\it{}\epsilon})', ' ', name{k}];
	legends{end+1} = ['imag({\it{}\epsilon})', ' ', name{k}];
end
legend(legends)

subplot(2,2,4)
hold
box on
grid on
title('\it{}\mu')
xlabel('\it{}f')
legends = {};
for k = 1:length(f)
	plot(f{k}, real(mu{k}), strcat(colors{mod(k, length(colors))}, '-'))
	plot(f{k}, imag(mu{k}), strcat(colors{mod(k, length(colors))}, '--'))
	legends{end+1} = ['real({\it{}\mu})', ' ', name{k}];
	legends{end+1} = ['imag({\it{}\mu})', ' ', name{k}];
end
legend(legends)
