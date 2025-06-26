function info (str, lvl, MAX)
% INFO - Display progress information
%
% Inputs:
%   str - formatted string of information to displayed
%   lvl - importance level. If lvl > MAX, it will not be displayed.
%   MAX - Maximum importance level. Default: 1.

if nargin < 3
  MAX = 1;
  if nargin < 2
    lvl = 0;
  end
end

if lvl <= MAX
  fprintf('%s: %s%s\n', datestr(now), repmat('  ',1,lvl), str);
end

end