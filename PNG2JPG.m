clc,clear;
file_path = 'F:\Meeting\202310MGH\tem\test\';                 
img_path_list = dir(strcat(file_path,'*.png'));
img_num = length(img_path_list);               

%%
for j = 1:img_num
    image_name = img_path_list(j).name;          
    photo = imread(strcat(file_path,image_name));
    name_png = strsplit(image_name,'.');
    name = string(name_png(1));
    filename = strcat(name,'.jpg');  
    jpg_file = fullfile('F:\Meeting\202310MGH\tem\topic1\',filename); 
    imwrite(photo,jpg_file,'jpg');     
end
%%
