function im1 = segementation_edge(im)
im = im*10^3;
%% binary
divide_loc = 0.5;
im(im<divide_loc)=0;
im(im>=divide_loc)=1;
%%
f = im;
%%
f1=bwareaopen(f,100);
%%
f2 = ~bwareaopen(~f1,100);
%% median filter
f3 = medfilt2(f2,[8,8]);
%%
f4=bwareaopen(f3,50);
%%
im1 = ~bwareaopen(~f4,50);