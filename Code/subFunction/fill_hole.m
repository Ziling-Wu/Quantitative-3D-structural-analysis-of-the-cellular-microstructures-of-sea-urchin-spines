function new_img = fill_hole(img)
BW = bwlabeln(img);
N = tabulate(BW(:));
%%
for i = 1:size(N,1)-1
    stas(i,1)=N(i+1,1);
    stas(i,2)=N(i+1,2);
end
img(BW ~= stas(stas(:,2)==max(stas(:,2)),1)) =0;
%%
inv_img = ~img;
BW = bwlabeln(inv_img);
N = tabulate(BW(:));
%%
for i = 1:size(N,1)-1
    stas(i,1)=N(i+1,1);
    stas(i,2)=N(i+1,2);
end
inv_img(BW ~= stas(stas(:,2)==max(stas(:,2)),1)) =0;
%
new_img =~inv_img;