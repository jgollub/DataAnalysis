%INSTALL_LINUX Script used to install TVReg on the Linux 64bit platform.
%
% Compiles the .c-files for the TVReg package.
%
% Tested with the gcc compiler.
%
% See readme.txt for further instructions.


any_error = false;

try
    cs = sprintf(['mex -largeArrayDims c/tvreg_upn_c.c c/tools.c c/tv_core.c ' ...
                  '']);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use tvreg_upn because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end


try
    cs = sprintf(['mex -largeArrayDims c/tvreg_gpbb_c.c c/tools.c c/tv_core.c ' ...
                  '']);
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
