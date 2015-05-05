% Retrieve the effective properties of a slab of metamaterial from
% its S parameters.
% 
% This function takes 6 to 8 input parameters:
%   f                           a vector of length m with the
%                               frequencies, in increasing order, at
%                               which the S parameters where
%                               determined;
%   S11                         a vector of length m with the complex
%                               values of simulated or measured
%                               reflection;
%   S21                         a vector of length m with the complex
%                               values of simulated or measured
%                               transmission;
%   d                           the thickness of the metamaterial slab;
%   prepad                      the thickness of vacuum between the
%                               input port and the slab of
%                               metamaterial;
%   postpad                     the thickness of vacuum between the
%                               slab of metamaterial and the output
%                               port;
%   branch_selection_criterion  (optional) the criterion used to select
%                               the appropriate branch (see below);
%   ensure_continuity           (optional) a boolean indicating if the
%                               function is allowed to jump from branch
%                               to branch to ensure the continuity of
%                               the refractive index (by default true).
% It returns 2 or 4 output parameters:
%   z                           the impedance of the metamaterial;
%   n                           the refractive index of the
%                               metamaterial;
%   epsilon                     (optional) the relative permittivity of
%                               the metamaterial;
%   mu                          (optional) the relative permeability of
%                               the metamaterial.
% 
% This function uses +i*omega*t convention. Make sure you adjust the
% values of S11 and S21 accordingly.
% 
% Three branch selection criteria are implemented:
%   'abs(real(n))'                the real part of the refractive
%                                 index diverges for all branches
%                                 except one at small frequency, the
%                                 branch giving the smallest absolute
%                                 value for the smallest frequency is
%                                 selected (this should be branch 0);
%   'imag(epsilon) and imag(mu)'  in homogeneous passive materials, the
%                                 imaginary parts of epsilon and mu are
%                                 positive (this is also the case for
%                                 metamaterials outside of resonances
%                                 and when the ratio of d and the
%                                 wavelength is sufficiently small) the
%                                 branch for which the most points
%                                 respect that condition is chosen;
%   integer                       it is also possible for the user to
%                                 select the branch by providing an
%                                 integer value.
% The default selection criterion is 'abs(real(n))' which, according to
% my experience, is more reliable than 'imag(epsilon) and imag(mu)'.
% 
% For 'abs(real(n))' and 'imag(epsilon) and imag(mu)' selection
% criteria to work, the frequency domain must extend outsite of the
% resonance region.
% 
% The continuity of the real part of the refractive index as well as
% that of its derivative are controlled. To avoid spurious branch
% jumps, the frequency resolution must be sufficient.
% 
% If a global variable DEBUG exists and is true, some information is
% shown during the execution of the function. The global variable
% NB_TESTED_BRANCHES (whose default value is 10) controls the number of
% branches tested to determine which one maximizes the number of points
% respecting the condition imag(epsilon) > 0 and imag(mu) > 0.
% 
% Examples:
%   [z, n] = retrieve(f, S11, S21, 2.5e-3, 10e-3, 10e-3)
%   [z, n, epsilon, mu] = retrieve(f, S11, S21, 2.5e-3, 10e-3, 10e-3, 'abs(real(n))', false)
% 
% For more information about the formalism used here, see
%   Chen et al., "Robust method to retrieve the constitutive
%   effective parameters of metamaterials", Phys. Rev. E, vol. 70, pp.
%   016608-1-016608-7.
% To learn why the imaginary parts of effective epsilon and mu are not
% always positive in metamaterials, see
%   T. Koschny, P. Markos, D. R. Smith, and C. M. Soukoulis, "Resonant
%   and antiresonant frequency dependence of the effective parameters
%   of metamaterials", Phys. Rev. E, vol. 68, 2003, pp. 065602(R)-1-
%   065602(R)-4.
% 
% Copyright (c) 2008,2009,2014 Stéphane Larouche

function [z, n, epsilon, mu] = retrieve(f, S11, S21, d, prepad, postpad, branch_selection_criterion, ensure_continuity)


% Verify the number of input parameters.
if nargin < 6 || nargin > 8
	error('retrieve take between 6 and 8 input parameters!')
end

% Verify the number of input parameters.
if ~(nargout == 2 || nargout == 4)
	error('retrieve returns 2 or 4 output parameters!')
end

% Set default input parameter values.
if nargin < 8
	ensure_continuity = true;
end
if nargin < 7
	branch_selection_criterion = 'abs(real(n))';
end

% Check input arguments.
if ~isvector(f) || ~isreal(f) || ~all(f(2:end)-f(1:end-1) > 0)
	error('f must be a vector of real numbers is increasing order!')
end
if ~isvector(S11) || any(size(S11) ~= size(f))
	error('S11 must be vector of the same size than f!')
end
if ~isvector(S21) || any(size(S21) ~= size(f))
	error('S21 must be vector of the same size than f!')
end
if ~isscalar(d) || ~isreal(d) || d <= 0
	error('d must be a real positive scalar!')
end
if ~isscalar(prepad) || ~isreal(prepad) || prepad < 0
	error('prepad must be a real positive scalar!')
end
if ~isscalar(postpad) || ~isreal(postpad) || postpad < 0
	error('postpad must be a real positive scalar!')
end
if ~(strcmp(branch_selection_criterion, 'abs(real(n))') || strcmp(branch_selection_criterion, 'imag(epsilon) and imag(mu)') || (isscalar(branch_selection_criterion) && (round(branch_selection_criterion) == branch_selection_criterion)))
	error('branch_selection_criterion must be ''abs(real(n))'', ''imag(epsilon) and imag(mu)'', or an integer!')
end
if ~islogical(ensure_continuity)
	error('ensure_continuity must be a boolean!')
end

% Check if the global variable DEBUG is defined. If not, the default
% behavior is to don't show debug information.
if ~isempty(whos('global', 'DEBUG'))
	global DEBUG
else
	DEBUG = false;
end

% Check if the global variable NB_TESTED_BRANCHES is defined. If not,
% the default value is 10.
if ~isempty(whos('global', 'NB_TESTED_BRANCHES'))
	global NB_TESTED_BRANCHES
else
	NB_TESTED_BRANCHES = 10;
end

% The speed of light
c = 299792458.0;

% Calculate wavenumber is vacuum and 1/(k_0*d) that will be useful
% later.
k_0 = 2.0*pi*f/c;
invk_0d = (k_0*d).^-1;

% Remove padding on both side of the metamaterial slab.
exp_prepad = exp(-1.0i*k_0*prepad);
exp_postpad = exp(-1.0i*k_0*postpad);
S11 = S11 .* exp_prepad.*exp_prepad;
S21 = S21 .* exp_prepad.*exp_postpad;

% Calculate the impedance.
z = sqrt(((1+S11).^2-S21.^2)./ ...
         ((1-S11).^2-S21.^2));

% Adjust the sign of the impedance to make sure that the real part is
% positive.
z = z.*(real(z)>=0.0) - z.*(real(z)<0.0);

% Calculate the exponential of i*n*k_0*d.
expink_0d = S21./(1.0-S11.*(z-1.0)./(z+1.0));

% Calculate the refractive index according to the branch selection
% criterion.
if ischar(branch_selection_criterion) && strcmp(branch_selection_criterion, 'abs(real(n))')
	
	branch = round(-imag(log(expink_0d(1)))/2/pi);
	n = calculate_n(f, expink_0d, invk_0d, branch, ensure_continuity);

elseif ischar(branch_selection_criterion) && strcmp(branch_selection_criterion, 'imag(epsilon) and imag(mu)')
	
	% Calculate branch zero.
	branch = 0;
	n = calculate_n(f, expink_0d, invk_0d, branch, ensure_continuity);
	nb_correct_pts = sum(abs(real(n).*imag(z)) <= imag(n).*real(z));
	
	% Check if another branch gives better results.
	for tested_branch = [1:NB_TESTED_BRANCHES,-1:-1:-NB_TESTED_BRANCHES]
		n_tested_branch = calculate_n(f, expink_0d, invk_0d, tested_branch, ensure_continuity);
		nb_correct_pts_tested_branch = sum(abs(real(n_tested_branch).*imag(z)) <= imag(n_tested_branch).*real(z));
		if nb_correct_pts_tested_branch > nb_correct_pts
			branch = tested_branch;
			n = n_tested_branch;
			nb_correct_pts = nb_correct_pts_tested_branch;
		end
	end
	
else
	
	branch = branch_selection_criterion;
	n = calculate_n(f, expink_0d, invk_0d, branch, ensure_continuity);

end

% If the number of output arguments is 4, calculate the permittivity
% and the permeability
if nargout == 4
	epsilon = n./z;
	mu = n.*z;
end

% If DEBUG is true, show the real part of the refractive index of all
% tested branches.
if DEBUG
	figure
	hold
	grid on
	box on
	xlabel('{\it{}f} (Hz)')
	ylabel('real({\it{}n})')
	
	correct = abs(real(n).*imag(z)) <= imag(n).*real(z);
	correct = correct.*correct.^-1;
	plot(f, real(n), 'r-', 'Linewidth', 2)
	plot(f, correct.*real(n), 'b-', 'Linewidth', 2)
	
	n_branch = calculate_n(f, expink_0d, invk_0d, -NB_TESTED_BRANCHES, false);
	n_min = 1.3*real(n_branch(end));
	n_branch = calculate_n(f, expink_0d, invk_0d, NB_TESTED_BRANCHES, false);
	n_max = 1.3*real(n_branch(end));
	
	ftag = f(1) + 0.9*(f(end)-f(1));
	
	for branch = -NB_TESTED_BRANCHES:NB_TESTED_BRANCHES
		n_branch = calculate_n(f, expink_0d, invk_0d, branch, false);
		correct = abs(real(n_branch).*imag(z)) <= imag(n_branch).*real(z);
		correct = correct.*correct.^-1;
		plot(f, real(n_branch), 'r-')
		plot(f, correct.*real(n_branch), 'b-')
		ntag = interp1(f, real(n_branch), ftag);
		text(ftag, ntag, sprintf('m = %i', branch), 'HorizontalAlignment', 'center', 'BackgroundColor', [1,1,1])
	end
	
	YLim = get(gca, 'YLim');
	YLim(1) = max(YLim(1), n_min);
	YLim(2) = min(YLim(2), n_max);
	set(gca, 'YLim', YLim);
	
	legend({'|real({\it{}n})imag({\it{}z})| > imag({\it{}n})real({\it{}z})', '|real({\it{}n})imag({\it{}z})| \leq imag({\it{}n})real({\it{}z})'}, 'Location', 'NorthWest')
end



function [n] = calculate_n(f, expink_0d, invk_0d, branch, ensure_continuity)

% The weight given to the continuity of n and its derivative.
WEIGHT_n = 1.0/3.0;
WEIGHT_dn = 1.0-WEIGHT_n;

if ensure_continuity
	n = zeros(size(f));
	branches = zeros(size(f));
	
	imaglogexpink_0d = imag(log(expink_0d));
	
	% For the first point, use the branch provided in the function call.
	branches(1) = branch;
	n(1) = invk_0d(1)*((imaglogexpink_0d(1)+2*pi*branch) - 1.0i*real(log(expink_0d(1))));
	
	% For the second point, ensure the continuity of n only since we
	% can't calculate the derivative.
	if length(expink_0d) > 1
		branches(2) = round((imaglogexpink_0d(1)-imaglogexpink_0d(2))/2/pi) + branches(1);
		n(2) = invk_0d(2)*((imaglogexpink_0d(2)+2*pi*branches(2)) - 1.0i*real(log(expink_0d(2))));
	end
	
	% For all other points, ensure the continuity of both n and dn
	% according to the weights provided by WEIGHT_n and WEIGHT_dn.
	for i_f = 3:length(expink_0d)
		branch_n = (imaglogexpink_0d(i_f-1)-imaglogexpink_0d(i_f))/2/pi;
		delta_f = (f(i_f)-f(i_f-1))/(f(i_f-1)-f(i_f-2));
		delta_n = imaglogexpink_0d(i_f-1)+pi*branches(i_f-1) - imaglogexpink_0d(i_f-2)-pi*branches(i_f-2);
		branch_dn = (delta_f*delta_n+imaglogexpink_0d(i_f-1) - imaglogexpink_0d(i_f))/2/pi;
		branches(i_f) = round(WEIGHT_n*branch_n+WEIGHT_dn*branch_dn) + branches(i_f-1);
		n(i_f) = invk_0d(i_f)*((imaglogexpink_0d(i_f)+2*pi*branches(i_f)) - 1.0i*real(log(expink_0d(i_f))));
	end

else
	n = invk_0d.*((imag(log(expink_0d))+2*pi*branch) - 1.0i*real(log(expink_0d)));
end
