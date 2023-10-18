clc,clear;
file_path = 'F:\Meeting\202310MGH\tem\test\';                 %待转化图片的文件路径
img_path_list = dir(strcat(file_path,'*.png'));%读取该路径下文件中所有png格式图片（可改）
img_num = length(img_path_list);               %获取图片数量

%%
for j = 1:img_num
    image_name = img_path_list(j).name;          %获取图片名字
    photo = imread(strcat(file_path,image_name));%读图片，strcat为字符串拼接
    name_png = strsplit(image_name,'.');
    name = string(name_png(1));
    filename = strcat(name,'.jpg');  %根据需求改名字(根据需求变)
    jpg_file = fullfile('F:\Meeting\202310MGH\tem\topic1\',filename); %新建文件
    imwrite(photo,jpg_file,'jpg');              %将png转换为bmp
end
%%