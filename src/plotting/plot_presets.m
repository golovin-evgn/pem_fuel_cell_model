function ps = plot_presets()
%
% plot_presets - gives predifined plots presets
%

ps.fs = 14;  % font size
ps.ms = 14;  % marker size
ps.lw = 2;   % line width
ps.lw_th = 1.5;   % line width thin
ps.lw_ti = 1;   % line width tiny


% colors
ps.color_black = [0 0 0];
ps.color_gray = [0.6 0.6 0.6];
ps.color_grey = [0.25 0.25 0.25];
ps.color_green = [50/255,205/255,50/255];
ps.color_tomato = [255/255, 99/255, 71/255];
ps.color_orange = [255/255, 140/255, 0];
ps.color_red = [220/255,20/255,60/255];
ps.color_blue = [0 71/255 171/255];

ps.color_gr = [ps.color_black; ps.color_gray; ps.color_grey];

% color map shades of red 
ps.color_red1 = [150/255, 0, 24/255];
ps.color_red_bright = [255/255, 193/255, 204/255];

% color map blue-green (schemecolor.com)
ps.color_bgmap = [24/255 32/255 96/255;
                  32/255 74/255 158/255;
                  43/255 152/255 211/255;
                  161/255 214/255 232/255;
                  170/255 236/255 217/255;
                  2/255 195/255 154/255;
                  0/255 168/255 150/255;
                  3/255 63/255 99/255];

ps.color_brmap = [24/255 32/255 96/255;
                  220/255 20/255 60/255;
                  43/255 152/255 211/255;
                  0/255 168/255 150/255;
                  150/255, 0, 24/255;
                  161/255 214/255 232/255;
                  170/255 236/255 217/255;
                  2/255 195/255 154/255;
                  3/255 63/255 99/255];

ps.color_rbmap =   [69, 73, 140; % blue
                    206, 40, 54; % red
                    69, 172, 89; % green
                    235, 191, 0; % yellow
                    238, 107, 49; % orange
                    79, 170, 209] /255; % light blue

ps.color_rblmap =  [158, 171, 214; % blue
                    224, 147, 145; % red
                    195, 232, 181; % green
                    255, 238, 184; % yellow
                    225, 181, 125; % orange
                    146, 206, 238] /255; % light blue

% line styles
ps.style1 = '--';
ps.style2 = '-';
ps.style3 = '-.';


