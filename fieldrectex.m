path(path,'/Users/Livvy/Documents/BME 8730/Field_II_ver_3_24_windows');
%path(path,'/home/user/field_II/m_files');
field_init(-1)
rect=[1 0/1000 0/1000 0 2/1000 0/1000 0 2/1000 5/1000 0 0/1000 5/1000 0 ...
1 2/1000 5/1000 1/1000 2.5/1000 0
1 -2/1000 0/1000 0 0/1000 0/1000 0 0/1000 5/1000 0 -2/1000 5/1000 0 ...
1 2/1000 5/1000 -1/1000 2.5/1000 0
1 -2/1000 -5/1000 0 0/1000 -5/1000 0 0/1000 0/1000 0 -2/1000 0/1000 0 ...
1 2/1000 5/1000 -1/1000 -2.5/1000 0
1 0/1000 -5/1000 0 2/1000 -5/1000 0 2/1000 0/1000 0 0/1000 0/1000 0 ...
1 2/1000 5/1000 1/1000 -2.5/1000 0];
center=[0 0 0];
focus=[0 0 70]/1000;
Th = xdc_rectangles (rect, center, focus);
xdc_show(Th)