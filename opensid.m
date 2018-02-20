%--------------------------------------------------------------------------
%                          openSID Toolbox 
%--------------------------------------------------------------------------
% openSID toolbox loading file
% Loads the toolbox and adds the necessary paths to Matlab in order to 
% access function, classes, documentation, demos, etc.
%
% Syntax
% run('opensid')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         Keith Soal 2017
%                 questions to keithsoal at gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Initialize the opensid toolbox
clc

% Add paths to Matlab
opensid_path = fileparts(which('opensid.m'));
addpath(opensid_path);
addpath(genpath([opensid_path filesep 'CoreSafe']));

if exist([opensid_path filesep 'Lib'], 'dir')
    addpath(genpath([opensid_path filesep 'Lib']));
end

% addpath(genpath([ae_path filesep 'GUI']));
% addpath(genpath([ae_path filesep 'Doc']));
% addpath(genpath([ae_path filesep 'Demo']));
% addpath(genpath([ae_path filesep 'Interface']));
% addpath(genpath([ae_path filesep 'LibStudents']));

clear opensid_path

disp(' ___________________________________________')
disp('|                                          |')
disp('|   openSID toolbox added to Matlab path!  |')
disp('|__________________________________________|')

   