clear;
clc;
fprintf('Subpixel Edge Detector v1.11 ------- (2014/08/12)\n');
fprintf('Initializing path..\n');
addpath('Synthetic');
addpath('EdgeDetector');
addpath('Statistics');
addpath('GUI');
fprintf('Launching GUI...\n\n');
global Image;
Image = [];
global restoredImage;
restoredImage = [];
global edges;
edges = EdgePixel;
SubpixelGUI;

% to test function circle
% tic; i=circle(29,29,15,15,7,100,150,100); toc; 

% to test function ramp
% tic; i=ramp(30,30,15,15,30,100,150); toc;

% to test smooth image
% tic; g=smooth(i); toc

% to test basic detector
% tic; e=basicDetector(i,1,1); toc;

% to add gaussian noise to an image
% n=noise(i,0.20);

% to measure error with sinthetic images
% fullStatCircle(190,190,100,100,70,100,150,100,0.1,10,2,1);
% fullStatCircle(29,29,15,15,7,100,150,100,0.1,10,2,1);

