% Plot multiple call example
figure(1); 
subplot(1,3,1)
mat1 = zeros(multiple_calls(1).gpl_match(1).cm.size(1),multiple_calls(1).gpl_match(1).cm.size(2));
mat1(multiple_calls(1).gpl_match(1).cm.index) = multiple_calls(1).gpl_match(1).cm.values;
imagesc(mat1); title('First Half'); 

subplot(1,3,2)
mat2 = zeros(multiple_calls(1).gpl_match(2).cm.size(1),multiple_calls(1).gpl_match(2).cm.size(2));
mat2(multiple_calls(1).gpl_match(2).cm.index) = multiple_calls(1).gpl_match(2).cm.values;
imagesc(mat2); title('Second Half')

subplot(1,3,3)
mat3 = zeros(multipleOut(1).gpl_match.cm.size(1),multipleOut(1).gpl_match.cm.size(2));
mat3(multipleOut(1).gpl_match.cm.index) = multipleOut(1).gpl_match.cm.values;
imagesc(mat3); title('Combined')

% Plot window overlap example
figure(2);
subplot(1,3,1)
mat1 = zeros(overlap_calls(3).gpl_match.cm.size(1),overlap_calls(3).gpl_match.cm.size(2));
mat1(overlap_calls(3).gpl_match.cm.index) = overlap_calls(3).gpl_match.cm.values;
imagesc(mat1); title('First Half'); 

subplot(1,3,2)
mat2 = zeros(overlap_calls(4).gpl_match.cm.size(1),overlap_calls(4).gpl_match.cm.size(2));
mat2(overlap_calls(4).gpl_match.cm.index) = overlap_calls(4).gpl_match.cm.values;
imagesc(mat2); title('Second Half')

subplot(1,3,3)
mat3 = zeros(overlapOut(2).gpl_match.cm.size(1),overlapOut(2).gpl_match.cm.size(2));
mat3(overlapOut(2).gpl_match.cm.index) = overlapOut(2).gpl_match.cm.values;
imagesc(mat3); title('Combined')