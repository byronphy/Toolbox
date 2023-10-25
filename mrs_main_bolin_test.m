close all;clear;clc

%% 导入数据和参数
FIDfile = ['D:\Data\哈医大多核20231023\'];
fileList = [14111];
% FIDfile = ['D:\Data\哈医大多核20231023\兔肝成像3\'];
% fileList = [18918:18928];
% fileList = [19085:19093];
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
spec = fftshift(fft(FID_rspace_zerofill,[],1),1); % 512x32x32x4
spec_sos = sqrt(sum(abs(spec).^2, 4)); % 四个通道模平方和的开方,512x32x32

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

% 选择图像空间的某一个体素
[~,tmpidx] = sort(img_MIP(:),'descend');
% Y_loc = ceil(tmpidx(1)/imgSize);
% X_loc = tmpidx(1)-(Y_loc-1)*imgSize;

X_loc = ceil(tmpidx(1)/imgSize);
Y_loc = tmpidx(1)-(X_loc-1)*imgSize;

%X_loc = 21;
%Y_loc = 27;

figure(1); 
subplot(131),imshow(img_all_peak,[]); title('all spectrum'); hold on, plot(X_loc,Y_loc,'or', 'MarkerSize',5, 'MarkerFaceColor', 'r');
subplot(132),imshow(img_MIP,[]); title('MIP');
subplot(133),imshow(img_index,[]); colorbar; title('MIP spectrum index');

%% Plot FID

FID_rspace_plot = FID_rspace_zerofill(:,X_loc,Y_loc,1); % 第1个通道图像某体素的FID
figure(2);
plot(time_range, real(FID_rspace_plot),'b', 'LineWidth', 2);hold on;
plot(time_range, imag(FID_rspace_plot),'r', 'LineWidth', 2);hold off;
legend('Real', 'Imaginary');
xlabel('Time (ms)');ylabel('Signal');title('FID (zerofilled)');

%% Apodization，消除填零导致的振铃伪影
% FID_rspace_apd = mrs_apod(FID_rspace_zerofill, 5000, 5);
% FID_rspace_zerofill = mrs_apod(FID_rspace_zerofill, 5000, 5);

%% 对波谱所有频谱校正(不宜过大，保证是磷谱才行)
% 先从最大值的几个体素的频移程度估计整体频移
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


%% 画出波谱
% ppm_range = linspace(BW_P/2,-BW_P/2,N)/freq_ref_P;
% begin_ppm = 25;
% end_ppm = -25;
% begin_index = round(-begin_ppm *N*freq_ref_P/BW_P+N/2);
% end_index = round(-end_ppm *N*freq_ref_P/BW_P+N/2);

ppm_range = linspace(-BW_P/2,BW_P/2,N)/freq_ref_P;
begin_ppm = -round(max(ppm_range));
end_ppm = round(max(ppm_range));
begin_index = round(begin_ppm *N*freq_ref_P/BW_P+N/2);
end_index = round(end_ppm *N*freq_ref_P/BW_P+N/2);

spec_loc = spec(:,X_loc,Y_loc,:); % 图像某体素的波谱，包括4个通道 512x4

figure(3);
subplot(2,2,1);
plot(ppm_range(begin_index:end_index), real(spec_loc(begin_index:end_index,1)),'b', 'LineWidth', 2); % 画第一个通道
legend('Real');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('Real Spectrum of ^{31}P, f_{ref}=51.8959 MHz (Channel 1)');
subplot(2,2,2);
plot(ppm_range(begin_index:end_index), imag(spec_loc(begin_index:end_index,1)),'r', 'LineWidth', 2);% 画第一个通道 
legend('Imaginary');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('Imaginary Spectrum of ^{31}P, f_{ref}=51.8959 MHz (Channel 1)');
subplot(2,2,3);
plot(ppm_range(begin_index:end_index), abs(spec_loc(begin_index:end_index,1)),'b', 'LineWidth', 2);
legend('Magnitude');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Signal');
title('Magnitude Spectrum of ^{31}P, f_{ref}=51.8959 MHz (Channel 1)');
subplot(2,2,4);
plot(ppm_range(begin_index:end_index), angle(spec_loc(begin_index:end_index,1)),'r', 'LineWidth', 2);
legend('Phase');set(gca,'XDir','reverse');xlabel('ppm');ylabel('Phase (rad)');
title('Phase Spectrum of ^{31}P, f_{ref}=51.8959 MHz (Channel 1)');

% 四通道几何平均
spec_sos = spec_shift(spec_sos,shift_value);
spec_loc_sos = spec_sos(:,X_loc,Y_loc); % 图像某体素的波谱 512
%spec_loc_sos = smooth(spec_loc_sos,8,'lowess'); % 平滑处理，窗宽8

figure(4);
plot(ppm_range(begin_index:end_index), spec_loc_sos(begin_index:end_index),'b', 'LineWidth', 2);
legend('Magnitude');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('Spectrum of ^{31}P, f_{ref}=51.8959 MHz (Combined Channel)');

%% 自动零阶相位校正(四个通道)
ppm_pc_range = [-20,20];
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
            [~,op_index] = max(phs_corr);
            op_phs = phs_range(op_index); 
            spec_pc(:,i,j,channel) = mrs_rephase(spec(:,i,j,channel), op_phs); % 512x32x32x4
            %spec_pc(:,i,j)=mrs_zeroPHC(spec(:,i,j));
        end
    end
end
size(spec_pc)

% 零阶相位校正后
spec_pc_loc = spec_pc(:,X_loc,Y_loc,:); % 图像某体素的波谱
figure(5);
subplot(2,2,1);
plot(ppm_range(begin_index:end_index), real(spec_pc_loc(begin_index:end_index,1)),'b', 'LineWidth', 2);
legend('Real');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('0th Phase Corrected Real Spectrum (Channel 1)');
subplot(2,2,2);
plot(ppm_range(begin_index:end_index), imag(spec_pc_loc(begin_index:end_index,1)),'r', 'LineWidth', 2);
legend('Imaginary');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('0th Phase Corrected Imaginary Spectrum (Channel 1)');
subplot(2,2,3);
plot(ppm_range(begin_index:end_index), abs(spec_pc_loc(begin_index:end_index,1)),'b', 'LineWidth', 2);
legend('Magnitude');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('0th Phase Corrected Magnitude Spectrum (Channel 1)');
subplot(2,2,4);
plot(ppm_range(begin_index:end_index), angle(spec_pc_loc(begin_index:end_index,1)),'r', 'LineWidth', 2);
legend('Phase');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Phase (rad)');title('0th Phase Corrected Imaginary Spectrum (Channel 1)');

% 四通道平均
spec_pc_sos = sqrt(sum(abs(spec_pc).^2, 4)); % 四个通道模平方和的开方,512x32x32
spec_pc_sos = spec_shift(spec_pc_sos,shift_value);
spec_pc_loc_sos = spec_pc_sos(:,X_loc,Y_loc); % 图像某体素的波谱
%spec_pc_loc_sos = smooth(spec_pc_loc_sos,8,'lowess'); % 平滑处理，窗宽8
figure(6);
plot(ppm_range(begin_index:end_index), spec_pc_loc_sos(begin_index:end_index),'b', 'LineWidth', 2);
legend('Magnitude');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('0th Phase Corrected Magnitude Spectrum (Combined Channel)');

%% 对特定位置体素的波谱作一阶相位校正(全部体素都做的话时间太长了)
spec_loc_1pc = zeros(N,size(spec,4));
options=optimset('TolX',1e-8,'MaxFunEvals',1e8, 'MaxIter',1e8);
for channel = 1:size(spec,4)
    phc=fminsearch(@(x) mrs_entropy(x,spec_loc(:,channel)), [0 0], options);
    phc0 = phc(1);
    phc1 = phc(2);
    ph=(phc0+phc1.*(1:N)/N);
    spec_loc_1pc(:,channel)=mrs_rephase(spec_loc(:,channel), ph');
end

figure(7);
subplot(2,2,1);
plot(ppm_range(begin_index:end_index), real(spec_loc_1pc(begin_index:end_index,1)),'b', 'LineWidth', 2);
legend('Real');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('1st Phase Corrected Real Spectrum (Channel 1)');
subplot(2,2,2);
plot(ppm_range(begin_index:end_index), imag(spec_loc_1pc(begin_index:end_index,1)),'r', 'LineWidth', 2);
legend('Real', 'Imaginary');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('1st Phase Corrected Imaginary Spectrum (Channel 1)');
subplot(2,2,3);
plot(ppm_range(begin_index:end_index), abs(spec_loc_1pc(begin_index:end_index,1)),'b', 'LineWidth', 2);
legend('Magnitude');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('1st Phase Corrected Magnitude Spectrum (Channel 1)');
subplot(2,2,4);
plot(ppm_range(begin_index:end_index), angle(spec_loc_1pc(begin_index:end_index,1)),'r', 'LineWidth', 2);
legend('Phase');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('1st Phase Corrected Phase Spectrum (Channel 1)');

% 四通道平均
spec_loc_1pc_sos = sqrt(sum(abs(spec_loc_1pc).^2, 2)); % 图像中心的波谱
spec_loc_1pc_sos = spec_shift(spec_loc_1pc_sos,shift_value);
spec_loc_1pc_sos = smooth(spec_loc_1pc_sos,8,'lowess'); % 平滑处理，窗宽8
figure(8);
plot(ppm_range(begin_index:end_index), spec_loc_1pc_sos(begin_index:end_index),'b', 'LineWidth', 2);
legend('Magnitude');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('1st Phase Corrected Magnitude Spectrum (Combined Channel)');

%% 基线校正，多项式函数效果不太好，似乎并不需要校正（原因是太多峰了，实际上可能是振铃伪影）
% y = spec_loc_1pc_sos;
% degree = 3;  % 三次多项式
% coeff = polyfit(ppm_range, y, degree);
% y_trend = polyval(coeff, ppm_range);
% y_detrend = y - y_trend;
% 
% % van_matrix = [time_range'.^3, time_range'.^2, time_range', ones(N,1)]; %三阶多项式
% % coeff = van_matrix \ y'; % 线性拟合
% % y_detrend = y - van_matrix * coeff;
% figure(12);
% plot(ppm_range(begin_index:end_index), spec_loc_1pc_sos(begin_index:end_index),'b', 'LineWidth', 2);hold on;
% plot(ppm_range(begin_index:end_index), y_trend(begin_index:end_index),'r', 'LineWidth', 2);hold on;
% plot(ppm_range(begin_index:end_index), y_detrend(begin_index:end_index),'black', 'LineWidth', 2);hold off;
% legend('Magnitude','Trend','Detrend');set(gca,'XDir','reverse');
% xlabel('ppm');ylabel('Signal');title('Detrended 1st Phase Corrected Spectrum');

%% 最小值法基线校正
% y = spec_loc_1pc_sos;
% y_trend = ones(N,1)*min(y);
% y_detrend = y - y_trend;
% spec_loc_1pc_sos_detrend = y_detrend;
% figure(12);
% plot(ppm_range(begin_index:end_index), spec_loc_1pc_sos(begin_index:end_index),'b', 'LineWidth', 2);hold on;
% plot(ppm_range(begin_index:end_index), y_trend(begin_index:end_index),'r', 'LineWidth', 2);hold on;
% plot(ppm_range(begin_index:end_index), y_detrend(begin_index:end_index),'black', 'LineWidth', 2);hold off;
% legend('Magnitude','Trend','Detrend');set(gca,'XDir','reverse');
% xlabel('ppm');ylabel('Signal');title('Detrended 1st Phase Corrected Spectrum');

%% 分峰拟合
spectrum = spec_loc_1pc_sos;
% PCr
PCr_ppm_range = [-0.8,2.7]; % [-3,3]
[PCr_fitted, PCr_pars] = FitLorentzPeak(spectrum,PCr_ppm_range);

% Pi
Pi_ppm_range = [3.28,5.62]; % [3.28,5.62]
[Pi_fitted, Pi_pars] = FitLorentzPeak(spectrum,Pi_ppm_range);

% alphaATP
alphaATP_ppm_range = [-9,-6]; % [-6,-9]
[alphaATP_fitted, alphaATP_pars] = FitLorentzPeak(spectrum,alphaATP_ppm_range);

% betaATP
betaATP_ppm_range = [-17.5,-14.5]; % [-17.5,-14.5]
[betaATP_fitted, betaATP_pars] = FitLorentzPeak(spectrum,betaATP_ppm_range);

% gammaATP
gammaATP_ppm_range = [-5,-2]; % [-5,-2]
[gammaATP_fitted, gammaATP_pars] = FitLorentzPeak(spectrum,gammaATP_ppm_range);

figure(13);
plot(ppm_range, spec_loc_1pc_sos,'b', 'LineWidth', 2);hold on;
plot(ppm_range,PCr_fitted,'r','LineWidth', 2);hold on;
plot(ppm_range,Pi_fitted,'g','LineWidth', 2);hold on;
plot(ppm_range,alphaATP_fitted,'c','LineWidth', 2);hold on;
plot(ppm_range,betaATP_fitted,'m','LineWidth', 2);hold on;
plot(ppm_range,gammaATP_fitted,'y','LineWidth', 2);hold off;
legend('Magnitude','PCr','Pi','alphaATP','betaATP','gammaATP');set(gca,'XDir','reverse');
xlabel('ppm');ylabel('Signal');title('Individual Peak Fitted');


%% 计算化学物质峰面积 S=A*FWHM
        % y0 = baseline amplitude
        % x0 = location of the peak 
        % fwhm = full width at half maximum
        % A = height of the peak 
PCr_Area = (PCr_pars(4)-PCr_pars(1))*PCr_pars(3);

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
    par_initials=[(data(1)+data(end))/2, x(index), 2.5, A_ini]; 
    [f, pars_fitted] = mrs_lorentzFit(par_initials, data, x');
    y_fitted = mrs_lorentzFun(1:N, pars_fitted(1), pars_fitted(2), pars_fitted(3), pars_fitted(4));
end

% %% 求面积函数
% function area = cal_Area(pars)
%     baseline = pars(1);
%     peak_loc = pars(2);
%     FWHM = pars(3);
%     peak_hight = pars(4);
%     if baseline < 0
%         area = peak_hight*FWHM;
%     else
%         area = (peak_hight-baseline)
% 
% end