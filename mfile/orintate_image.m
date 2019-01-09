% 
% UCAIR, University of Utah
% Ye Tian, phye1988@gmail.com
%
% image = orintate_image(image,orin)
%
% Change orintation of image.
%
%
% Input:
%      orin - orintation flag got from orintation_detection

function image = orintate_image(image,orin)

switch orin
    case 1
    case 2
        image = rot90(image,1);
    case 3
        image = rot90(image,2);
    case 4
        image = rot90(image,3);
    case 5
        image = flipud(image);
    case 6
        image = flipud(rot90(image,1));
    case 7
        image = flipud(rot90(image,2));
    case 8
        image = flipud(rot90(image,3));
end