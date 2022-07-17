function B = imgaussfilt_local(A,varargin)
%IMGAUSSFILT_LOCAL 2-D Gaussian filtering of images
%
%   This function was adapted within the fastACI toolbox to be 'forward 
%   compatible': It is used when your MATLAB version does not contain the 
%   imgaussfilt.m (introduced starting from R2015). In this script, some of
%   the input parameters are fixed.
%
%   B = imgaussfilt(A) filters image A with a 2-D Gaussian smoothing kernel
%   with standard deviation of 0.5. A can have any number of dimensions.
%
%   B = imgaussfilt_local(A,SIGMA) filters image A with a 2-D Gaussian smoothing
%   kernel with standard deviation specified by SIGMA. SIGMA can be a
%   scalar or a 2-element vector with positive values. If sigma is a
%   scalar, a square Gaussian kernel is used.
%
%   B = imgaussfilt_local(___,Name,Value,...) filters image A with a 2-D
%   Gaussian smoothing kernel with Name-Value pairs used to control aspects
%   of the filtering. 
%
%   Parameters include:
%
%   'FilterSize'    -   Scalar or 2-element vector, of positive, odd 
%                       integers that specifies the size of the Gaussian
%                       filter. If a scalar Q is specified, then a square
%                       Gaussian filter of size [Q Q] is used.
%
%                       Default value is 2*ceil(2*SIGMA)+1.
%
%   'Padding'       -   String or character vector or numeric scalar that
%                       specifies padding to be used on image before
%                       filtering. 
%
%                       If a scalar (X) is specified, input image values
%                       outside the bounds of the image are implicitly
%                       assumed to have the value X. 
%
%                       If a string is specified, it can be 'replicate',
%                       'circular' or 'symmetric'. These options are
%                       analogous to the padding options provided by
%                       imfilter. 
%
%                       'replicate'
%                       Input image values outside the bounds of the image
%                       are assumed equal to the nearest image border
%                       value.
%
%                       'circular'
%                       Input image values outside the bounds of the image
%                       are computed by implicitly assuming the input image
%                       is periodic.
%
%                       'symmetric'
%                       Input image values outside the bounds of the image
%                       are computed by mirror-reflecting the array across
%                       the array border. 
%
%                       Default value is 'replicate'.
%
%   'FilterDomain'  -   String or character vector that specifies domain in
%                       which filtering is performed. It can be 'spatial',
%                       'frequency' or 'auto'. For 'spatial', convolution
%                       is performed in the spatial domain, for 'frequency',
%                       convolution is performed in the frequency domain
%                       and for 'auto', convolution may be performed in
%                       spatial or frequency domain based on internal
%                       heuristics.
%
%                       Default value is 'auto'.
%
%   Example
%   ---------
%   % Smooth an image with Gaussian filters of increasing standard deviations
%       I = imread('cameraman.tif');
%
%       subplot(2,2,1), imshow(I), title('Original Image');
%
%       Iblur = imgaussfilt_local(I, 2);
%       subplot(2,2,2), imshow(Iblur)
%       title('Gaussian filtered image, \sigma = 2')
%
%       Iblur = imgaussfilt_local(I, 4);
%       subplot(2,2,3), imshow(Iblur)
%       title('Gaussian filtered image, \sigma = 4')
%
%       Iblur = imgaussfilt_local(I, 6);
%       subplot(2,2,4), imshow(Iblur)
%       title('Gaussian filtered image, \sigma = 6')
%

narginchk(1, Inf);

sigma       = varargin{1}*[1 1]; % options.Sigma;
hsize       = 2*ceil(2*sigma)+1; % from computeFilterSizeFromSigma

if strcmp(varargin{2},'Padding')
    padding = varargin{3};
end
domain      = 'auto';

[domain, separableFlag] = il_chooseFilterImplementation(A, hsize, domain);

switch domain
    case 'spatial'
        B = il_spatialGaussianFilter(A, sigma, hsize, padding, separableFlag);
        
    case 'frequency'
        B = il_frequencyGaussianFilter(A, sigma, hsize, padding);
        
    otherwise
        assert(false, 'Internal Error: Unknown filter domain');
end

end

%--------------------------------------------------------------------------
% Spatial Domain Filtering
%--------------------------------------------------------------------------
function A = il_spatialGaussianFilter(A, sigma, hsize, padding, separableFlag)

dtype = class(A);

if separableFlag

    [hCol,hRow] = il_createSeparableGaussianKernel(sigma, hsize);
    
    switch class(A)
        case {'int32','uint32'}
            A = double(A);
        case {'uint8','int8','uint16','int16'}
            A = single(A);
        case {'single','double'}
            % No-op
        otherwise
            assert(false,'Unexpected datatype');
    end

    [~, padSize] = il_computeSizes(A, hsize);

    A = il_filterDoubleSeparableWithConv(A, hCol,hRow, hsize, padSize, padding);
    
    if ~isa(A,dtype)
        A = cast(A,dtype);
    end
    
else
    h = createGaussianKernel(sigma, hsize);
    
    A = imfilter(A, h, padding, 'conv', 'same');
end
                    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [finalSize, pad] = il_computeSizes(a, hSize)

rank_a = ndims(a);
rank_h = numel(hSize);

% Pad dimensions with ones if filter and image rank are different
size_h = [hSize ones(1,rank_a-rank_h)];
size_a = [size(a) ones(1,rank_h-rank_a)];

%Same output
finalSize = size_a;

%Calculate the number of pad pixels
filter_center = floor((size_h + 1)/2);
pad = size_h - filter_center;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = il_filterDoubleSeparableWithConv(a, hcol,hrow, hSize, padSize, padding)

sameSize = 1;

if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end

imageSize = size(a);

nonSymmetricPadShift = 1-mod(hSize,2);
ndimsH = numel(hSize);
prePadSize = padSize;
prePadSize(1:ndimsH) = padSize(1:ndimsH)-nonSymmetricPadShift;

if sameSize && any(nonSymmetricPadShift == 1)
    a = padarray_algo(a, prePadSize, method, padVal, 'pre');
    a = padarray_algo(a, padSize, method, padVal, 'post');
else
    a = padarray_algo(a,padSize,method,padVal,'both');
end


if ismatrix(a)
    result = conv2(hcol, hrow, a,'valid');
else % Stack behavior
    result = zeros(imageSize, 'like',a);
    for i = 1:size(a,3)
       result(:,:,i) = conv2(hcol, hrow, a(:,:,i),'valid');
    end
end

end

function TF = useSeparableFiltering(A, hsize)

isKernel1D = any(hsize==1);

minKernelElems = getSeparableFilterThreshold(class(A));

TF = ~isKernel1D && prod(hsize) >= minKernelElems;

end

% function TF = useIPPL(A, outSize)
% 
% prefFlag = images.internal.useIPPLibrary();
% 
% if ~isImageIPPFilterType(class(A))
%     TF = false;
%     return;
% end
% 
% tooBig = isImageTooBigForIPPFilter(A, outSize);
% 
% TF = prefFlag && ~tooBig;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hcol,hrow] = il_createSeparableGaussianKernel(sigma, hsize)

isIsotropic = sigma(1)==sigma(2) && hsize(1)==hsize(2);

hcol = images.internal.createGaussianKernel(sigma(1), hsize(1));

if isIsotropic
    hrow = hcol;
else
    hrow = images.internal.createGaussianKernel(sigma(2), hsize(2));
end

hrow = reshape(hrow, 1, hsize(2));

end

%--------------------------------------------------------------------------
% Frequency Domain Filtering
%--------------------------------------------------------------------------
function A = il_frequencyGaussianFilter(A, sigma, hsize, padding)

sizeA = size(A);
dtype = class(A);
outSize = sizeA;

A = il_padImage(A, hsize, padding);

h = images.internal.createGaussianKernel(sigma, hsize);

% cast to double to preserve precision unless single
if ~isfloat(A)
    A = double(A);
end

fftSize = size(A);
if ismatrix(A)
    A = ifft2( fft2(A) .* fft2(h, fftSize(1), fftSize(2)), 'symmetric' );
else
    fftH = fft2(h, fftSize(1), fftSize(2));
    
    dims3toEnd = prod(fftSize(3):fftSize(end));
    
    %Stack behavior
    for n = 1 : dims3toEnd
        A(:,:,n) = ifft2( fft2(A(:,:,n), fftSize(1), fftSize(2)) .* fftH, 'symmetric' );
    end
end

% cast back to input type
if ~strcmp(dtype,class(A))
    A = cast(A, dtype);
end

A = il_unpadImage(A, outSize);

end

%--------------------------------------------------------------------------
% Common Functions
%--------------------------------------------------------------------------
function [domain, separableFlag] = il_chooseFilterImplementation(A, hsize, domain)

ippFlag = 1; % useIPPL(A, size(A));

separableFlag = 0; % useSeparableFiltering(A, hsize);

if strcmp(domain, 'auto')
    domain = il_chooseFilterDomain(A, hsize, ippFlag);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, padSize] = il_padImage(A, hsize, padding)

padSize = il_computePadSize(size(A), hsize);

if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end

A = padarray_algo(A, padSize, method, padVal, 'both');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function padSize = il_computePadSize(sizeA, sizeH)

rankA = numel(sizeA);
rankH = numel(sizeH);

sizeH = [sizeH ones(1,rankA-rankH)];

padSize = floor(sizeH/2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = il_unpadImage(A, outSize)

start = 1 + size(A) - outSize;
stop  = start + outSize - 1;

subCrop.type = '()';
subCrop.subs = {start(1):stop(1), start(2):stop(2)};

for dims = 3 : ndims(A)
    subCrop.subs{dims} = start(dims):stop(dims);
end

A = subsref(A, subCrop);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function domain = il_chooseFilterDomain(A, hsize, ippFlag)
%IL_CHOOSEFILTERDOMAIN chooses appropriate domain for 2-D filtering
%   domain = chooseFilterDomain(A, hsize, ippFlag) determines whether
%   filtering image A with 2-D kernel of size hsize is faster in the
%   frequency domain or spatial domain with IPP availability specified by
%   ippFlag.
%
%   Note that this is only for internal use by imgaussfilt.

if ippFlag
    bigImageThreshold = 2.25e6;         %1.5k x 1.5k
    if isa(A,'single')
        bigKernelThreshold = 1.44e4;    %120 x 120
    else
        bigKernelThreshold = 4e4;       %200 x 200    
    end
else
    bigImageThreshold = 1e6;            %1k x 1k
    bigKernelThreshold = 1.6e3;         %40 x 40
end

Asize = [size(A,1) size(A,2)];
imageIsBig = prod(Asize)>=bigImageThreshold;
kernelIsBig = prod(hsize)>=bigKernelThreshold;

if imageIsBig && kernelIsBig
    domain = 'frequency';
else
    domain = 'spatial';
end

end