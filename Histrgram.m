close all
clear all
clc
% create a histogram for κmax values on the regulatory framework and its efficacy in preventing the
%buildup of atmospheric carbon within the context of thermal energy storage.
gamma = [1600 1800 2000 2200 2400 2600];
event_trigger = [876 1109 3104; 910 1364 3213;  934 1326 3321; 988 1404 3388; 1068 1472 3418; 1094 1733 3458];
self_trigger = [1109 1364 2326 231404 1472 1733 2456 2925 3110];
bar(gamma, event_trigger, 1)
% hold on
