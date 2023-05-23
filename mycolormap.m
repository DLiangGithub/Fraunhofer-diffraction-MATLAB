function map = mycolormap(wavelength)
%MYCOLORMAP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

if nargin < 2
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

I = wave2rgb(wavelength);
r = linspace(0,I(1),m);
g = linspace(0,I(2),m);
b = linspace(0,I(3),m);
map = [r' g' b'];
end

