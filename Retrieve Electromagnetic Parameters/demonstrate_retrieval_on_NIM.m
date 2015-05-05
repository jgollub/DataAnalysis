% Demonstrate the retrieval for a negative refractive index
% metamaterial consisting of split ring resonators and continuous
% wires simulated using Comsol 4.4. For details about the
% structure, see
%   D. R. Smith, D. C. Vier, Th. Koschny, and C. M. Soukoulis, 
%  "Electromagnetic parameter retrieval from inhomogeneous
%  metamaterials", Phys. Rev. E, vol. 71, 2005, pp. 036617-1-036617-11.

global DEBUG
DEBUG = true;

close all

d = 0.0025;
pad = 0.0015;
criterion = 'abs(real(n))';

[f, S11, S21] = read_COMSOL4('NIM.txt');

[z, n, epsilon, mu] = retrieve(f, S11, S21, d, pad, pad, criterion, true);

show_retrieval('Demonstration of parameter retrieval on NIM', {f}, {z}, {n}, {epsilon}, {mu})
