
%------------------------------------
% data_image.m
% function Y = data_image
%------------------------------------
% Date   : Mar. 24, 2014
% Author : Shiori Watanabe(s1190068@u-aizu.ac.jp)
% Note   : This code initialize and reading data matrices for
%           Non-negative Matrix Factorization(NMF) and this inherite
%           Hoyer's code. (http://www.cs.helsinki.fi/u/phoyer/software.html)
%
%------------------------------------

function Y = data_image
% data_image - read face image data from cbcl database
%

global imloadfunc;

%------------------------------------
% This is where the cbcl face images reside
%------------------------------------
thepath = '../data/cbcl-face-database/face/';


%------------------------------------
% Create the data matrix
%------------------------------------
Y = zeros(19*19,2429);

%------------------------------------
%Read the directory listing
%------------------------------------
D = dir(thepath);

%------------------------------------
%Step through each image, reading it into the data matrix
%Note: The (+2) is just to avoid '.' and '..' entries
%------------------------------------
fprintf('Reading in the iamges...\n');
for i = 1:2429,
    switch imloadfunc
        case 'pgma_read'
            I = utils.pgma_read([thepath D(i+2).name]);
        otherwise,
            I = imread([thepath D(i+2).name]);
    end
    Y(:,i) = reshape(I,[19*19, 1]);
    if rem(i,100) == 1, fprintf('[%d/24]',floor(i/100)); end
end

fprintf('\n');


%------------------------------------
% Same preprocessing as Lee and Seung
%------------------------------------
Y = Y - mean(Y(:)); % Y's mean to be zero
Y = Y / sqrt(mean(Y(:).^2)); % Y's power to be 1 (nomalization)
Y = Y + 0.25;
Y = Y* 0.25;
Y = min(Y,1); % setting to be max value is 1 (fix values(1<v) to be 1) 
Y = max(Y,0); % setting to be min value is 0 (fix values(0>v) to be 0)

% Additionally, this is required to avoid having any exact zeros:
% (divergence objective cannot handle them!)
Y = max(Y, 1e-4);




