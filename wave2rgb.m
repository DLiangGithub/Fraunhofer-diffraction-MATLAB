% 输入指定波长，转换成RGB数值，并输出
% 输入波长范围：380~780 nm
 
function I = wave2rgb(lam)
 

 lambda = lam*1e9;
if  (lambda >= 380.0) && (lambda < 440.0)
    r = -1.0 * (lambda - 440.0) / (440.0 - 380.0);
    g = 0.0;
    b = 1.0;
elseif (lambda >= 440.0) && (lambda < 490.0)
    r = 0.0;
    g = (lambda - 440.0) / (490.0 - 440.0);
    b = 1.0;
elseif (lambda >= 490.0) && (lambda < 510.0)
    r = 0.0;
    g = 1.0;
    b = -1.0 * (lambda - 510.0) / (510.0 - 490.0);
elseif (lambda >= 510.0) && (lambda < 580.0)
    r = (lambda - 510.0) / (580.0 - 510.0);
    g = 1.0;
    b = 0.0;
elseif (lambda >= 580.0) && (lambda < 645.0)
    r = 1.0;
    g = -1.0 * (lambda - 645.0) / (645.0 - 580.0);
    b = 0.0;
elseif (lambda >= 645.0) && (lambda <= 780.0)
    r = 1.0;
    g = 0.0;
    b = 0.0;
else
    r = 0.0;
    g = 0.0;
    b = 0.0;
end
 
% 在可见光谱的边缘处强度较低。
if  (lambda >= 380.0) && (lambda < 420.0)
    attenuation = 0.3 + 0.7 * (lambda - 380) / (420 - 380);
    r = r * attenuation;
    g = 0.0;
    b = 1.0 * attenuation;
elseif(lambda >= 701.0) && (lambda < 780.0)
    attenuation = 0.30 + 0.70 * (780.0 - lambda) / (780.0 - 700.0);
    r = r * attenuation;
    g = 0.0;
    b = 0.0;
end
r = round(r*255);
g = round(g*255);
b = round(b*255);
 
I = [r,g,b];

end