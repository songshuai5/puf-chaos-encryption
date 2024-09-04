clear all
close all
clc

% 1. 读取彩色图像
img = imread('your_image.png');  % 加载彩色图像
img = double(img);  % 将图像转换为 double 类型

% 2. 定义图像分块大小
block_size = [32, 32];  % 分块大小为 32x32
[rows, cols, channels] = size(img);

% 3. 计算填充后的图像尺寸以确保整除
padded_rows = ceil(rows / block_size(1)) * block_size(1);
padded_cols = ceil(cols / block_size(2)) * block_size(2);

% 对图像进行填充
padded_img = zeros(padded_rows, padded_cols, channels);
padded_img(1:rows, 1:cols, :) = img;

% 更新分块数
row_blocks = padded_rows / block_size(1);
col_blocks = padded_cols / block_size(2);

% 3. 生成 PUF 密钥
rng(123); % 固定随机种子以确保实验可重复
puf_key_blocks = cell(row_blocks, col_blocks, channels);

for c = 1:channels
    for i = 1:row_blocks
        for j = 1:col_blocks
            puf_key_blocks{i, j, c} = uint8(randi([0 255], block_size)); % 转换为 uint8 类型
        end
    end
end

% 4. 生成 Lorenz 混沌序列
sigma = 10;
rho = 28;
beta = 8/3;
dt = 0.01;  % 时间步长
num_steps = block_size(1) * block_size(2);  % 每个块中的像素数

% 初始条件
x = 0.1;
y = 0.1;
z = 0.1;

lorenz_sequence_blocks = cell(row_blocks, col_blocks, channels);

for c = 1:channels
    for i = 1:row_blocks
        for j = 1:col_blocks
            x = 0.1; y = 0.1; z = 0.1; % 重置初始条件
            lorenz_sequence = zeros(1, num_steps);
            for k = 1:num_steps
                dx = sigma * (y - x) * dt;
                dy = (x * (rho - z) - y) * dt;
                dz = (x * y - beta * z) * dt;
                x = x + dx;
                y = y + dy;
                z = z + dz;
                lorenz_sequence(k) = mod(abs(x), 1);  % 取 x 作为混沌序列的一部分
            end
            lorenz_sequence = uint8(255 * lorenz_sequence);  % 转换为 uint8 类型
            lorenz_sequence_blocks{i, j, c} = reshape(lorenz_sequence, block_size);
        end
    end
end
%%
% 5. 图像加密
encrypted_img = zeros(padded_rows, padded_cols, channels, 'uint8');  % 使用 uint8 类型存储加密图像

for c = 1:channels
    for i = 1:row_blocks
        for j = 1:col_blocks
            % 提取当前块
            block = uint8(padded_img((i-1)*block_size(1)+1:i*block_size(1), ...
                                     (j-1)*block_size(2)+1:j*block_size(2), c)); % 转换为 uint8 类型
            
            % 进行双重加密：PUF 和 Lorenz 混沌
            encrypted_block = bitxor(block, puf_key_blocks{i, j, c});
            encrypted_block = bitxor(encrypted_block, lorenz_sequence_blocks{i, j, c});
            
            % 存储加密后的块
            encrypted_img((i-1)*block_size(1)+1:i*block_size(1), ...
                          (j-1)*block_size(2)+1:j*block_size(2), c) = encrypted_block;
        end
    end
end

% 6. 显示并保存加密后的图像
figure;
subplot(1, 2, 1);
imshow(uint8(img)); 
title('Original Image');

subplot(1, 2, 2);
imshow(encrypted_img(1:rows, 1:cols, :));  % 显示裁剪回原始大小的加密图像
title('Encrypted Image');
imwrite(encrypted_img(1:rows, 1:cols, :), 'encrypted_image.png');
%%
% 7. 图像解密
decrypted_img = zeros(padded_rows, padded_cols, channels, 'uint8');  % 使用 uint8 类型存储解密图像

for c = 1:channels
    for i = 1:row_blocks
        for j = 1:col_blocks
            % 提取加密块
            encrypted_block = encrypted_img((i-1)*block_size(1)+1:i*block_size(1), ...
                                            (j-1)*block_size(2)+1:j*block_size(2), c);
            
            % 进行双重解密：PUF 和 Lorenz 混沌
            decrypted_block = bitxor(encrypted_block, lorenz_sequence_blocks{i, j, c});
            decrypted_block = bitxor(decrypted_block, puf_key_blocks{i, j, c});
            
            % 存储解密后的块
            decrypted_img((i-1)*block_size(1)+1:i*block_size(1), ...
                          (j-1)*block_size(2)+1:j*block_size(2), c) = decrypted_block;
        end
    end
end

% 8. 显示并保存解密后的图像
figure;
subplot(1, 2, 1);
imshow(encrypted_img(1:rows, 1:cols, :));
title('Encrypted Image');

subplot(1, 2, 2);
imshow(decrypted_img(1:rows, 1:cols, :));  % 显示裁剪回原始大小的解密图像
title('Decrypted Image');
%imwrite(decrypted_img(1:rows, 1:cols, :), 'decrypted_image.png');

% 9. 验证加解密是否成功
if isequal(uint8(img), decrypted_img(1:rows, 1:cols, :))
    disp('Decryption successful, the image is restored correctly.');
else
    disp('Decryption failed, the image restoration is incorrect.');
end
