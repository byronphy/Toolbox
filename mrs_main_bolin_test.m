close all;clear;clc

%% 导入数据和参数
FIDfile = ['D:\Data\哈医大多核20231023\'];
fileList = [14111];
%fileList = [14169];
% FIDfile = ['D:\Data\哈医大多核20231023\兔肝成像1\'];
% fileList = [18918:18928];
% fileList = [19034:19043];
avg =20;%10  
mreadP = 256;%FIDpiont3
BW_P = 5000;%Hz 
freq_ref_P = 51.8959; % MHz
ppmRange = 25; % 
isAdjustP = 1;  %
IsRegis=1;
ByteOrder = 2;
gamma = [42.576,40.053,11.262,17.235];% H F Na P [42.576,40.053,11.262,17.235]
num_nucleus = 4;
% H F Na paramter
mreadH = 128;    %点数  
spoke = 128;    %相位编码步数
% P spectrum parameter
X = 16;Y = 16; Z=1;
imgSize = 32;
padfactor = 1;

%% GetFIdFromBidary 这段代码没看懂，读取FID的，需要确认
% P sepctrum recon
ns = 20;
choosens = 20;
fidpoints = mreadP*2; % 考虑了实部和虚部

FID = zeros([mreadP,X,Y,Z,4]);
for fileIdx_P = fileList
    id=fopen([FIDfile,num2str(fileIdx_P),'\prefid'], 'r', 'l');			%use little endian format if asked for
    for k = 1:Z
        for j = 1:Y
            for i = 1:X
                spokeidx = i+(j-1)*X+(k-1)*X*Y;
                fseek(id,(((32+mreadH*2*4*4)*(num_nucleus-1)+(32+fidpoints*4*4)*1)*ns*(spokeidx-1)+((32+mreadH*2*4*4)*(num_nucleus-1)+(32+fidpoints*4*4)*1)*(choosens-1)), 'bof');
                for num_U = num_nucleus : -1 : 1
                    %%读fid1_1文件时，将下面这段包头信息注释掉
                    dataindex     = fread(id, 1,'*ubit16');  %数据索引
                    rxChannelNum  = fread(id, 1,'*ubit8');  %接收通道数
                    rxID          = fread(id, 1,'*ubit8');  %接收机ID
                    dataLength    = fread(id, 1,'*ubit32');  %长度，包括fid头和fid数据,以字节为单位
                    phaseIndex    = fread(id, 1,'*ubit32');  %相位编码索引
                    freqIndex     = fread(id, 1,'*ubit16');  %频率编码索引
                    isAdding      = fread(id, 1,'*ubit8');  %是否累加
                    isDummyScan   = fread(id, 1,'*ubit8');  %是否空扫
                    channelEnable = fread(id, 1,'*ubit64');  %通道可用标志，仅低八位
                    reserve       = fread(id, 1,'*ubit64');  %保留位
                    if rxID == 4  % 表示31P的接收线圈（共4个通道）
                        for cn=1:rxChannelNum
                            fidall = fread(id, fidpoints, 'int32');
                            fidreal = fidall(1:2:end,1);
                            fidimag = fidall(2:2:end,1);
                            FID_one(1:(fidpoints / 2),i,j,k,cn) = fidreal.'+1i*fidimag.';
                        end
                    end
                end
            end
        end
    end
    FID = FID+FID_one;
end
FID = FID/avg/length(fileList)/(mreadP*(1+padfactor))*1000;
FID = squeeze(FID); % 256x16x16x4
FID = permute(FID,[1,3,2,4]);
FID = FID(:,end:-1:1,:,:);

%% k空间填零，根据H和P的gamma来在k空间中继续填零
FID_kspace = padarray(FID,[0 (imgSize-X)/2 (imgSize-Y)/2 0],0,'both'); % 图像空间XY两侧填零,256x32x32x4

gamma = [42.576,40.053,11.262,17.235];% H F Na P
if isAdjustP ==1
    n = round(gamma(1)/gamma(4)*imgSize);
else 
    n = imgSize;
end
padding = n-imgSize;
FID_kspace = padarray(FID_kspace,[0 round(padding/2) round(padding/2) 0],0,'both'); % isAdjustP: 256x80x80x4

%% 2D-FFT 得到图像空间
FID_rspace = zeros([size(FID_kspace,1),imgSize,imgSize,size(FID_kspace,4)]);
for c = 1 : size(FID_rspace,4) % 四个通道
    for t = 1:mreadP % 对每个FID时刻点的k空间
        tmp = ifftshift(ifft2(ifftshift(squeeze(FID_kspace(t,:,:,c)))));
        FID_rspace(t,:,:,c) = tmp(round(padding/2)+1:round(padding/2)+imgSize,round(padding/2)+1:round(padding/2)+imgSize);
    end
end

%% 对图像空间的FID进行预处理（具体每个体素的FID）
% 时域填零
FID_rspace_zerofill = cat(1,FID_rspace,zeros(round(padfactor*mreadP),imgSize,imgSize,4)); % 第一个维度填零，时间维度，信号后面补零，512x32x32x4
N = size(FID_rspace_zerofill,1);
time_range = linspace(0,N/BW_P*1000,N);

%% 对所有FID做FFT得到波谱
% 多体素波谱
spec = fftshift(fft(FID_rspace_zerofill,[],1),1); % 512x32x32x4
%spec_sos = sqrt(sum(abs(spec).^2, 4)); % 四个通道模平方和的开方,512x32x32
spec_sos = sqrt(sum(abs(spec(:,:,:,1:3)).^2, 4)); % 前3个通道模平方和的开方,512x32x32

%% 对31P幅度谱值求和，找最大的位置，并对附近的体素都进行频谱平移校正
img_all_peak = squeeze(sum(spec_sos, 1));
[img_MIP, index] = max(spec_sos,[],1);
img_MIP = squeeze(img_MIP); index = squeeze(index);
mask = img_all_peak > 1;
mask = imfill(mask,'holes');   
img_index = index .* mask;
img_all_peak = img_all_peak-min(img_all_peak(:));
img_all_peak = img_all_peak/max(abs(img_all_peak(:)));
img_MIP = img_MIP/max(abs(img_MIP(:)));

% 选择图像空间的某一个体素（选MIP图中最大的地方，理论上是PCr含量最大的地方）
[~,tmpidx] = sort(img_MIP(:),'descend');
Y_loc = ceil(tmpidx(1)/imgSize);
X_loc = tmpidx(1)-(Y_loc-1)*imgSize;

% X_loc = ceil(tmpidx(1)/imgSize);
% Y_loc = tmpidx(1)-(X_loc-1)*imgSize;

% X_loc = 24;
% Y_loc = 21;

figure(1); 
subplot(131),imshow(img_all_peak,[]); title('all spectrum'); hold on, plot(X_loc,Y_loc,'or', 'MarkerSize',5, 'MarkerFaceColor', 'r');
subplot(132),imshow(img_MIP,[]); title('MIP');
subplot(133),imshow(img_index,[]); colorbar; title('MIP spectrum index');

%% 单体素FID：(X_loc,Y_loc)
FID_rspace_loc = FID_rspace_zerofill(:,X_loc,Y_loc,:); % 第1个通道图像某体素的FID, 512x4

% Plot 单通道
figure(2);
subplot(2,2,1);
plot(time_range, real(FID_rspace_loc(:,1)),'b', 'LineWidth', 2);hold on;
plot(time_range, imag(FID_rspace_loc(:,1)),'r', 'LineWidth', 2);hold off;
legend('Real', 'Imaginary');
xlabel('Time (ms)');ylabel('Signal');title('Zero-filled FID (Channel 1)');
subplot(2,2,2);
plot(time_range, real(FID_rspace_loc(:,2)),'b', 'LineWidth', 2);hold on;
plot(time_range, imag(FID_rspace_loc(:,2)),'r', 'LineWidth', 2);hold off;
legend('Real', 'Imaginary');
xlabel('Time (ms)');ylabel('Signal');title('Zero-filled FID (Channel 2)');
subplot(2,2,3);
plot(time_range, real(FID_rspace_loc(:,3)),'b', 'LineWidth', 2);hold on;
plot(time_range, imag(FID_rspace_loc(:,3)),'r', 'LineWidth', 2);hold off;
legend('Real', 'Imaginary');
xlabel('Time (ms)');ylabel('Signal');title('Zero-filled FID (Channel 3)');
subplot(2,2,4);
plot(time_range, real(FID_rspace_loc(:,4)),'b', 'LineWidth', 2);hold on;
plot(time_range, imag(FID_rspace_loc(:,4)),'r', 'LineWidth', 2);hold off;
legend('Real', 'Imaginary');
xlabel('Time (ms)');ylabel('Signal');title('Zero-filled FID (Channel 4)');

%% Apodization，消除填零导致的振铃伪影
% FID_rspace_apd = mrs_apod(FID_rspace_zerofill, 5000, 5);
% FID_rspace_zerofill = mrs_apod(FID_rspace_zerofill, 5000, 5);

%% 多体素波谱频移校正：更新spec, 512x32x32x4
% 先从PCr含量最大的几个体素的频移程度估计整体频移
voxel_num = 10;
shift_values = zeros(voxel_num,1);
for n = 1:voxel_num
    voxel_y = ceil(tmpidx(n)/imgSize);
    voxel_x = tmpidx(n)-(voxel_y-1)*imgSize;
    spec_sos_voxel = spec_sos(:,voxel_x,voxel_y);
    spec_sos_voxel = smooth(spec_sos_voxel,8,'lowess'); % 平滑幅度谱，窗宽8
    [~,Idx] = sort(spec_sos_voxel,'descend'); % 降序排列波谱，幅值最大的为中心频率
    shift_values(n) = round(N/2)-Idx(1)+1; % 对第n个体素进行频移校正：移动距离为shift_values(i),10
    % tmp_ =circshift(spec_sos_voxel,shift_values(n));
end
shift_value = round(mean(shift_values));

% 对所有通道/所有体素进行校正
for channel = 1:size(spec,4)
    for j = 1:size(spec,3)
        for i = 1:size(spec,2)
            spec(:,i,j,channel) = spec_shift(spec(:,i,j,channel),shift_value);
        end
    end
end

%% 单体素波谱spec_loc, 512x4
spec_loc = squeeze(spec(:,X_loc,Y_loc,:)); % 图像某体素的波谱，包括4个通道 512x4
% 联合通道
smooth_width = 4;
% spec_loc_sos = sqrt(sum(abs(spec_loc).^2, 2)); % 四通道
spec_loc_sos = sqrt(sum(abs(spec_loc(:,1:4)).^2, 2));% 前4个通道模平方和的开方,512x32x32
%spec_loc_sos = smooth(spec_loc_sos,smooth_width,'lowess'); % 平滑处理，窗宽8

spec_loc_sos13 = sqrt(sum(abs(spec_loc(:,1:3)).^2, 2));% 前3个通道模平方和的开方,512x32x32
%spec_loc_sos13 = smooth(spec_loc_sos13,smooth_width,'lowess'); % 平滑处理，窗宽8

spec_loc_sos12 = sqrt(sum(abs(spec_loc(:,1:2)).^2, 2));% 前3个通道模平方和的开方,512x32x32
%spec_loc_sos12 = smooth(spec_loc_sos12,smooth_width,'lowess'); % 平滑处理，窗宽8

ppm_range = linspace(-BW_P/2,BW_P/2,N)/freq_ref_P;
% begin_ppm = -round(max(ppm_range)); % 全频段
% end_ppm = round(max(ppm_range));
begin_ppm = -25; % 关注的频段
end_ppm = 25;
begin_index = round(begin_ppm *N*freq_ref_P/BW_P+N/2);
end_index = round(end_ppm *N*freq_ref_P/BW_P+N/2);

% Plot 单通道
figure(3);
subplot(2,2,1);
plot(ppm_range(begin_index:end_index), real(spec_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on; 
plot(ppm_range(begin_index:end_index), real(spec_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on; 
plot(ppm_range(begin_index:end_index), real(spec_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1); hold on;
plot(ppm_range(begin_index:end_index), real(spec_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1); hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('Real Spectrum of ^{31}P');
subplot(2,2,2);
plot(ppm_range(begin_index:end_index), imag(spec_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('Imaginary Spectrum of ^{31}P');
subplot(2,2,3);
plot(ppm_range(begin_index:end_index), abs(spec_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('Magnitude Spectrum of ^{31}P');
subplot(2,2,4);
plot(ppm_range(begin_index:end_index), angle(spec_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold on;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Phase (rad)');
title('Phase Spectrum of ^{31}P');


% Plot 联合通道
figure(4);
plot(ppm_range(begin_index:end_index), spec_loc_sos(begin_index:end_index),'b', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), spec_loc_sos13(begin_index:end_index),'r', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), spec_loc_sos12(begin_index:end_index),'g', 'LineWidth', 2);hold off;
legend('1-4','1-3','1-2');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('Magnitude Spectrum of ^{31}P (Combined Channel)');
%%
% % 单通道
% figure(4);
% plot(ppm_range(begin_index:end_index), abs(spec_loc(begin_index:end_index,2)),'b', 'LineWidth', 2);
% legend('Magnitude');set(gca,'XDir','reverse');
% xlabel('ppm');ylabel('Signal');title('Magnitude Spectrum of ^{31}P (Channel 2)');

%% 多体素：零阶相位校正，spec(512,32,32,4) --> spec_pc(512,32,32,4) 
ppm_pc_range = [-20,20];  %感兴趣物质所在的区间，只针对这部分信号，然后全局校正
index_pc_range = round(ppm_pc_range *N*freq_ref_P/BW_P+N/2);

deg_step = 1; % deg
phs_range = (-180:deg_step:180)/180*pi;
phs_len = length(phs_range);
phs_corr = zeros(1,phs_len);
spec_pc = zeros(N,imgSize,imgSize, size(spec,4));
for channel = 1:size(spec,4)
    for j = 1:size(spec,3)
        for i = 1:size(spec,2)
            % 找某一范围内的最大实部谱对应的零相位
            for phs_index = 1:phs_len
                phs = phs_range(phs_index);
                phs_corr(phs_index) = sum(real(mrs_rephase(spec(index_pc_range,i,j,channel),phs)));
            end
            % 零阶相位校正
            [~,op_index] = max(phs_corr);
            op_phs = phs_range(op_index); 
            spec_pc(:,i,j,channel) = mrs_rephase(spec(:,i,j,channel), op_phs); % 512x32x32x4
            %spec_pc(:,i,j)=mrs_zeroPHC(spec(:,i,j));
        end
    end
end
size(spec_pc)

% 选单体素看看效果
spec_pc_loc = squeeze(spec_pc(:,X_loc,Y_loc,:)); % 图像某体素的波谱，包括4个通道 512x4
% 联合通道
smooth_width = 4;
spec_pc_loc_sos = sqrt(sum(abs(spec_pc_loc).^2, 2)); % 四通道
%spec_pc_loc_sos = smooth(spec_pc_loc_sos,smooth_width,'lowess'); % 平滑处理，窗宽8
spec_pc_loc_sos13 = sqrt(sum(abs(spec_pc_loc(:,1:3)).^2, 2));% 前3个通道模平方和的开方,512x32x32
%spec_pc_loc_sos13 = smooth(spec_pc_loc_sos13,smooth_width,'lowess'); % 平滑处理，窗宽8
spec_pc_loc_sos12 = sqrt(sum(abs(spec_pc_loc(:,1:2)).^2, 2));% 前2个通道模平方和的开方,512x32x32
%spec_pc_loc_sos12 = smooth(spec_pc_loc_sos12,smooth_width,'lowess'); % 平滑处理，窗宽8

ppm_range = linspace(-BW_P/2,BW_P/2,N)/freq_ref_P;
% begin_ppm = -round(max(ppm_range)); % 全频段
% end_ppm = round(max(ppm_range));
begin_ppm = -25; % 关注的频段
end_ppm = 25;
begin_index = round(begin_ppm *N*freq_ref_P/BW_P+N/2);
end_index = round(end_ppm *N*freq_ref_P/BW_P+N/2);

% Plot 单通道
figure(5);
subplot(2,2,1);
plot(ppm_range(begin_index:end_index), real(spec_pc_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on; 
plot(ppm_range(begin_index:end_index), real(spec_pc_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on; 
plot(ppm_range(begin_index:end_index), real(spec_pc_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1); hold on;
plot(ppm_range(begin_index:end_index), real(spec_pc_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1); hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('0th-order Phase-Corrected Real spec_pctrum of ^{31}P');
subplot(2,2,2);
plot(ppm_range(begin_index:end_index), imag(spec_pc_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_pc_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_pc_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_pc_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('0th-order Phase-Corrected Imaginary spec_pctrum of ^{31}P');
subplot(2,2,3);
plot(ppm_range(begin_index:end_index), abs(spec_pc_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_pc_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_pc_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_pc_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('0th-order Phase-Corrected Magnitude spec_pctrum of ^{31}P');
subplot(2,2,4);
plot(ppm_range(begin_index:end_index), angle(spec_pc_loc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_pc_loc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_pc_loc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_pc_loc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold on;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Phase (rad)');
title('0th-order Phase-Corrected Phase spec_pctrum of ^{31}P');

% Plot 联合通道
figure(6);
plot(ppm_range(begin_index:end_index), spec_pc_loc_sos(begin_index:end_index),'b', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), spec_pc_loc_sos13(begin_index:end_index),'r', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), spec_pc_loc_sos12(begin_index:end_index),'g', 'LineWidth', 2);hold off;
legend('1-4','1-3','1-2');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('0th-order Phase-Corrected Magnitude Spectrum of ^{31}P (Combined Channel)');


%% 单体素：一阶相位校正，spec_loc(512,4) --> spec_loc_linearpc(512,4) 
spec_1oc_linearpc = zeros(size(spec_loc));
options=optimset('TolX',1e-8,'MaxFunEvals',1e8, 'MaxIter',1e8);
for channel = 1:size(spec,4)
    phc=fminsearch(@(x) mrs_entropy(x,spec_loc(:,channel)), [0 0], options);
    phc0 = phc(1);
    phc1 = phc(2);
    ph=(phc0+phc1.*(1:N)/N);
    spec_1oc_linearpc(:,channel)=mrs_rephase(spec_loc(:,channel), ph');
end

% 联合通道 square root of sum of squared (SOS)
smooth_width = 4;
spec_1oc_linearpc_sos = sqrt(sum(abs(spec_1oc_linearpc).^2, 2)); % 四通道
spec_1oc_linearpc_sos = smooth(spec_1oc_linearpc_sos,smooth_width,'lowess'); % 平滑处理，窗宽8
spec_1oc_linearpc_sos13 = sqrt(sum(abs(spec_1oc_linearpc(:,1:3)).^2, 2));% 前3个通道模平方和的开方,512x32x32
spec_1oc_linearpc_sos13 = smooth(spec_1oc_linearpc_sos13,smooth_width,'lowess'); % 平滑处理，窗宽8
spec_1oc_linearpc_sos12 = sqrt(sum(abs(spec_1oc_linearpc(:,1:2)).^2, 2));% 前3个通道模平方和的开方,512x32x32
spec_1oc_linearpc_sos12 = smooth(spec_1oc_linearpc_sos12,smooth_width,'lowess'); % 平滑处理，窗宽8

ppm_range = linspace(-BW_P/2,BW_P/2,N)/freq_ref_P;
% begin_ppm = -round(max(ppm_range)); % 全频段
% end_ppm = round(max(ppm_range));
begin_ppm = -25; % 关注的频段
end_ppm = 25;
begin_index = round(begin_ppm *N*freq_ref_P/BW_P+N/2);
end_index = round(end_ppm *N*freq_ref_P/BW_P+N/2);

% Plot 单通道
figure(7);
subplot(2,2,1);
plot(ppm_range(begin_index:end_index), real(spec_1oc_linearpc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on; 
plot(ppm_range(begin_index:end_index), real(spec_1oc_linearpc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on; 
plot(ppm_range(begin_index:end_index), real(spec_1oc_linearpc(begin_index:end_index,3)),'g-', 'LineWidth', 1); hold on;
plot(ppm_range(begin_index:end_index), real(spec_1oc_linearpc(begin_index:end_index,4)),'m-', 'LineWidth', 1); hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('1st-order Phase-Corrected Real spec_pctrum of ^{31}P');
subplot(2,2,2);
plot(ppm_range(begin_index:end_index), imag(spec_1oc_linearpc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_1oc_linearpc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_1oc_linearpc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), imag(spec_1oc_linearpc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('1st-order Phase-Corrected Imaginary spec_pctrum of ^{31}P');
subplot(2,2,3);
plot(ppm_range(begin_index:end_index), abs(spec_1oc_linearpc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_1oc_linearpc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_1oc_linearpc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), abs(spec_1oc_linearpc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold off;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('1st-order Phase-Corrected Magnitude spec_pctrum of ^{31}P');
subplot(2,2,4);
plot(ppm_range(begin_index:end_index), angle(spec_1oc_linearpc(begin_index:end_index,1)),'b-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_1oc_linearpc(begin_index:end_index,2)),'r-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_1oc_linearpc(begin_index:end_index,3)),'g-', 'LineWidth', 1);hold on;
plot(ppm_range(begin_index:end_index), angle(spec_1oc_linearpc(begin_index:end_index,4)),'m-', 'LineWidth', 1);hold on;
legend('Channel 1','Channel 2','Channel 3','Channel 4');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Phase (rad)');
title('1st-order Phase-Corrected Phase spec_pctrum of ^{31}P');

% Plot 联合通道
figure(8);
plot(ppm_range(begin_index:end_index), spec_1oc_linearpc_sos(begin_index:end_index),'b', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), spec_1oc_linearpc_sos13(begin_index:end_index),'r', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), spec_1oc_linearpc_sos12(begin_index:end_index),'g', 'LineWidth', 2);hold off;
legend('1-4','1-3','1-2');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('1st-order Phase-Corrected Magnitude Spectrum of ^{31}P (Combined Channel)');

%% 基线校正，多项式函数效果不太好，似乎并不需要校正（原因是太多峰了，实际上可能是振铃伪影）
y = spec_1oc_linearpc_sos;
degree = 3;  % 三次多项式
coeff = polyfit(ppm_range, y, degree);
y_trend = polyval(coeff, ppm_range');
y_detrend = y - y_trend;

[y_min,min_index] = min(y_detrend);
y_min_list = ones(size(y))*y_min;
if y_min <0
    y_trend = y_trend + y_min;
    y_detrend = y_detrend - y_min;
end

spec_loc_detrend = y_detrend;

figure(13);
plot(ppm_range(begin_index:end_index), spec_1oc_linearpc_sos(begin_index:end_index),'b', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), y_trend(begin_index:end_index),'r', 'LineWidth', 2);hold on;
plot(ppm_range(begin_index:end_index), y_detrend(begin_index:end_index),'black', 'LineWidth', 2);hold off;
legend('Magnitude','Trend','Detrend');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('Detrended 1st Phase Corrected Spectrum');

%% 分峰拟合
spectrum = spec_loc_detrend; 
% PCr
PCr_ppm_range = [-1.41,2.16]; % [-3,3]
[PCr_fitted, PCr_pars] = FitLorentzPeak(spectrum,PCr_ppm_range);

% Pi
Pi_ppm_range = [4.8,6.12]; % [3.29,4.24]
[Pi_fitted, Pi_pars] = FitLorentzPeak(spectrum,Pi_ppm_range);

% alphaATP
alphaATP_ppm_range = [-8.2,-6.3]; % [-6,-9]
[alphaATP_fitted, alphaATP_pars] = FitLorentzPeak(spectrum,alphaATP_ppm_range);

% betaATP
betaATP_ppm_range = [-17,-14.8]; % [-17.5,-14.5]
[betaATP_fitted, betaATP_pars] = FitLorentzPeak(spectrum,betaATP_ppm_range);

% gammaATP
gammaATP_ppm_range = [-2.54,-1.22]; % [-2.3,-1.79]
[gammaATP_fitted, gammaATP_pars] = FitLorentzPeak(spectrum,gammaATP_ppm_range);

figure(14);
plot(ppm_range, spectrum,'b', 'LineWidth', 2);hold on;
plot(ppm_range,PCr_fitted,'r','LineWidth', 2);hold on;
plot(ppm_range,Pi_fitted,'g','LineWidth', 2);hold on;
plot(ppm_range,alphaATP_fitted,'c','LineWidth', 2);hold on;
plot(ppm_range,betaATP_fitted,'m','LineWidth', 2);hold on;
plot(ppm_range,gammaATP_fitted,'y','LineWidth', 2);hold off;
legend('Magnitude','PCr','Pi','alphaATP','betaATP','gammaATP');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('Individual Peak Fitted');

%% 计算化学物质峰面积 S=A*FWHM
        % x0 = location of the peak 
        % fwhm = full width at half maximum
        % A = height of the peak 
index_range = 1:N;
%PCr_Area = PCr_pars(2)*PCr_pars(3)
PCr_Area = trapz(index_range,PCr_fitted);
Pi_Area = trapz(index_range,Pi_fitted);
alphaATP_Area = trapz(index_range,alphaATP_fitted);
betaATP_Area = trapz(index_range,betaATP_fitted);
gammaATP_Area = trapz(index_range,gammaATP_fitted);












%% 对sos幅度谱的频率轴平移校正
function spec_sos_shifted = spec_shift(spec_sos,shift_value)
    spec_sos_shifted = zeros(size(spec_sos,1),size(spec_sos,2),size(spec_sos,3));
    for voxel_y = 1:size(spec_sos,3)
        for voxel_x = 1:size(spec_sos,2)
            spec_sos_voxel = spec_sos(:,voxel_x,voxel_y);
            spec_sos_shifted(:,voxel_x,voxel_y) = circshift(spec_sos_voxel,shift_value);
        end
    end
end

%% 一阶相位校正部分的函数
function f = mrs_entropy(x,spectrum)
% Entropy is defined as the normalized derivative of the NMR spectral data
% ARGS :
% x = phc0 and phc1
% spectrum = a spectrum before automatic first-order phase correction 
% RETURNS : 
% f = entropy value (Using the first derivative)

    %initial parameters
    stepsize=1; 
    func_type=1;

    %dephase
    L=length(spectrum);
    phc0=x(1);
    phc1=x(2);
    
    % linear phase
    n=length(spectrum);
    ph=(phc0+phc1.*(1:n)/n);
    
    spectrum=real(mrs_rephase(spectrum, ph));

    % Calculation of first derivatives 
    if (func_type == 1)
        ds1 = abs((spectrum(3:L)-spectrum(1:L-2))/(stepsize*2));
    else
        ds1 = ((spectrum(3:L)-spectrum(1:L-2))/(stepsize*2)).^2;
    end  
    p1 = ds1./sum(ds1);

    %Calculation of Entropy
    [M,K]=size(p1);
    for i=1:M
        for j=1:K
            if (p1(i,j)==0)%in case of ln(0)
               p1(i,j)=1; 
            end
        end
    end
    h1  = -p1.*log(p1);
    H1  = sum(h1);
    %Calculation of penalty
    Pfun	= 0.0;
    as      = spectrum - abs(spectrum);
    sumas   = sum(sum(as));
    if (sumas < 0)
       Pfun = Pfun + sum(sum((as./2).^2));
    end
    P = 1000*Pfun;

    % The value of objective function
    f = H1+P;

end

%% 分峰拟合函数
function [y_fitted, pars_fitted] = FitLorentzPeak(spectrum, ppm_range)
    % pars_fitted: 
        % y0 = baseline amplitude
        % x0 = location of the peak 
        % fwhm = full width at half maximum
        % A = height of the peak 
    N = length(spectrum);
    freq_ref_P = 51.8959; BW_P = 5000;
    index_range = round(ppm_range *N*freq_ref_P/BW_P+N/2);
    x = index_range(1):index_range(2);
    data = spectrum(x);
    [A_ini,index] = max(data);
    %par_initials=[(data(1)+data(end))/2, x(index), 2.5, A_ini]; 
    par_initials=[x(index), 2.5, A_ini]; 
    [f, pars_fitted] = lorentzFit(par_initials, data, x');
    %y_fitted = mrs_lorentzFun(1:N, pars_fitted(1), pars_fitted(2), pars_fitted(3), pars_fitted(4));
    y_fitted = lorentzFun(1:N, pars_fitted(1), pars_fitted(2), pars_fitted(3));
end

function  [y_fitted, pars_fitted]= lorentzFit(pars0, data, x)
    %options = optimset('TolX',1e-8,'MaxFunEvals',1e8, 'MaxIter',1e8);
    options = optimset('TolX',1e-8,'MaxFunEvals',1e8, 'MaxIter',10000);
    %[pars_fitted, ~] = fminsearch(@(pars) se_fun(pars, data, x), pars0,options);
    %y_fitted = mrs_lorentzFun(x, pars_fitted(1), pars_fitted(2), pars_fitted(3), pars_fitted(4));
    [pars_fitted, ~] = fminsearch(@(pars) error_fun(pars, data, x), pars0,options);
    y_fitted = lorentzFun(x, pars_fitted(1), pars_fitted(2), pars_fitted(3));
end

function se = error_fun(pars,data, x)
% SE_FUN defines the objective function to be minimised 
% ARGS :
% pars = parameters to be estimated ([y0 x0 fwhm A])
% data = data to be fitted with a lorentzian function
% x = input vectors 
    x0=pars(1);   % location of the peak 
    fwhm=pars(2); % full width at half maximum
    A=pars(3);    % height of the peak 
    % lorentzian model
    est_peak= lorentzFun(x, x0, fwhm, A );
    % squared error
    se = sum((est_peak-data).^2);
end

function y = lorentzFun(x, x0, fwhm, A)
     m=(x-x0)*2./fwhm;
     y= A.*(1./(1+m.^2));
end

