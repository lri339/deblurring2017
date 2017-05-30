function [ im_skeleton ] = skeletonize( im )
%SKELETONIZE 
%   LI
% shrinks edges in image so they are only 1-2 pixels wide

im_skeleton = bwmorph(im, 'skel');

% alternatively
% nested for loop too slow,use MATLAB default function instead
% [H W] = size(im);
% im_skeleton = im;
% 
% for i = 1:(H-1)
% 	for j = 1:(W-1)
% 		if im_skeleton(i:(i+1), j:(j+1)) == [1, 1; 0, 1];
% 			im_skeleton(i:(i+1), j:(j+1)) = [1 0; 0 1];
% 		elseif im_skeleton(i:(i+1), j:(j+1)) == [1 0 ; 1 1];
% 			im_skeleton(i:(i+1), j:(j+1)) = [1 0; 0 1];
% 		elseif im_skeleton(i:(i+1), j:(j+1)) == [1 1 ; 1 0 ];
% 			im_skeleton(i:(i+1), j:(j+1)) = [0 1; 1 0 ];
% 		elseif im_skeleton(i:(i+1), j:(j+1)) == [0 1; 1 1];
% 			im_skeleton(i:(i+1), j:(j+1)) = [0 1; 1 0 ];
%         end
% 	end
% end



end


