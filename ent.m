function [jhist,jent,MI] = ent(varargin)
%takes images of equal size and pixel value ranges and returns the
%joint histogram, joint entropy, and mutual information. No error checking
%used, as this function is built for speed. For additional speed,
%instructions on how to remove what you do not need are in the comments
%below
%made my Jason Agne 10MAY13

%%assumes 256x256x... (pixel values range from 0-255)
%%this assumption can be changed here if desired
dimen = 256;

x = numel(varargin{1});
t = 1:x;
J = varargin{1};
yy = zeros(nargin,x);

yy(1,:) = J(:)+1;

for i=2:nargin
    J = varargin{i};
    yy(i,:) = (dimen^(i-1))*J(:);
end

xx = sum(yy,1);
xx = sort(xx);

yy(1,1:x-1) = xx(2:x);

zz = yy(1,:) - xx;
zz(x) = 1;
zz = t(zz ~=0);

yy = xx(zz);

t = numel(zz);

zz(2:t) = zz(2:t)-zz(1:t-1);

%normalizes joint historgam output. Remove for counts.
xx = zz/x;

%if you do not need the joint histogram or mutual information, comment this
%out
dimen = repmat(dimen,[1 nargin]);

%if you do not need the joint histogram but do need the mutual information,
%remove jhist from the returned variables, but do NOT comment this out.
%if you need neither the joint histogram nor the mutual information,
%comment out the following two lines
jhist = zeros(dimen);
jhist(yy) = xx;

%if you do not need joint entropy, comment out the following line and
%remove jent from the list of returned variables
jent = -sum(xx.*log2(xx));

%if you do not need mutual information, comment out all of the below and
%remove MI from the list of returned variables
xx = zeros(nargin,dimen(1));
for i=1:nargin
    %notsum is a separate function
    xx(i,:) = squeeze(notsum(jhist,i));
end

xx(xx>1e-12) = (xx(xx>1e-12)).^-1;

MI = jhist;
for i=1:nargin
    x = 1:nargin;
    dimen(i) = 1;
    x([2 i]) = x([i 2]);
    MI = MI.*repmat(permute(xx(i,:),x),dimen);
    dimen(i) = max(dimen(:));
end
    

MI = sum(jhist(:).*log2(MI(:)+(MI(:)<=1e-12)));

%additional speed may be possible by allocating some of these variables
%outside of the function. 
end