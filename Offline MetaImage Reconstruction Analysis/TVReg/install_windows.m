%INSTALL_WINDOWS Script used to install TVReg on the Windows platform.
%
% Compiles the mex files for the TVReg package.
%
% Tested with the Lcc compiler bundled with Matlab.
%
% See readme.txt for further instructions.


%
% If you want to be able to break the execution of the programs, try to
% set CTRLCBREAK = 1, which uses a non-documented MATLAB API.
% If you do not have libut.lib in your Matlab dir, try CTRLCBREAK = 2.
% Default.
%
CTRLCBREAK=0;

if CTRLCBREAK==0
    sbreak = '';
elseif CTRLCBREAK==1
    sbreak = ['-DLIBUT -L"' matlabroot '\extern\lib\win32\lcc" -llibut'];
elseif CTRLCBREAK==2
    sbreak = ['-DLIBUT -Lexternlib -llibut'];
else
    error('Not a valid option for CTRLCBREAK')
end


any_error = false;

try
    cs = sprintf('mex %s c\\tvreg_upn_c.c c\\tools.c c\\tv_core.c ',sbreak);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use tvreg_upn because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end


try
    cs = sprintf('mex %s c\\tvreg_gpbb_c.c c\\tools.c c\\tv_core.c ',sbreak );
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use tvreg_gpbb because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end

if any_error == false
    disp('Install completed successfully.')
else
    disp('Installation did not complete successfully.')
end
