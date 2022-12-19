%set all paths of sub-directories. Using relative path-names, assuming that
%this script is in the same folder as PLDSExample.m
clc
fprintf('\nSetting up paths for code-package \n')

addpath([cd])
addpath([cd '/compiled'])
addpath([cd '/mex'])
