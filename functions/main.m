% Main code 1
% P2P market, single simulation
% W. Ananduta
% 24/07/2019

clear all
close all
clc


% Generate case study
n = 4; %number of agents
sp = 1; %sparsity, in (0,1)
ty = 1; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = 0; %transaction costs: (0)none, (1)uniform, or (2)heterogenous
np = gen_case(n,sp,ty,tc);

% Distributed algorithm
[s,np] = dist_alg(np);