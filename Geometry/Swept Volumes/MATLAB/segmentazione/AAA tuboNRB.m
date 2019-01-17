clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% preprocess the volume
% definisce una sezione di ingresso tramite il piano perpendicolare
% all'inizio del primo segmento non essendo ancora stato definito
[St,v]=pre_proc_VOL(St,v);
% rimuove quello di cui non abbiamo informazione topologica
St=post_process(St,v);
figure;hold on
plot_segment(St,v)
hold on
plotc(St)
% processa tutti i segmenti per ottenerne nurbs volume
St=create_NRB_VOL(St,v);
%%
plot_vol_and_segment(St,v);
    %% Bifurcation 2 processata completamente
    load('bif_2_completa.mat')
    %% artereia coronaria ridotta processata completamente
    load('ArtCor_finita.mat')
